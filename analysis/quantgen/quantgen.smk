#!python3

# =================
# author: Sara A. Carioscia, Biology Dept., Johns Hopkins University
# email: scarios1@jhu.edu
# last updated: January 26, 2025
# aim: Compute heritability and genetic correlation for recombination, aneuploidy, and 
#       published fertility-related traits. Process input files as necessary. 
# =================

# Usage: snakemake --snakefile quantgen.smk --use-conda --profile ~/code/rockfish_smk_profile -p

configfile: "config.yaml"

chromosomes = [str(i) for i in range(1, 24)]

# Create all heritability and genetic correlation results 
rule all:
    input:
        "results/heritability_published_merged.txt",
    	"results/genetic_correlation_merged.txt",
    	expand("results/pheWAS_results_{rsid}.tsv", rsid=[config["rsid"]]),
    	"results/queried_snps_across_traits.tsv",


# -------- Step 1: Steps to standardize Natera summary stats and supporting files for use in LDSC ------- #

rule process_dbsnp: 
	"""Generate cpra and rsid table for dbsnp."""
	input: 
		dbsnp=config["dbsnp"],
	output:
		cpra2rsid_info="results/intermediate_files/dbsnp151_hg38_info.txt"
	resources:
		time="6:00:00",
		mem_mb=1000
	threads: 16
	shell: 
		'bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t%ID\n" {input.dbsnp} | bgzip -@ {threads} > {output.cpra2rsid_info}'


rule rename_summary_stats:
    """
    For aneuploidy and recombination summary stats, downselect and reorder columns.
    For public datasets, copy existing summary stats.
    """
    input:
        summary_stats=lambda wildcards: config["summary_stats"][wildcards.name]["file"]
    output:
        summary_stats_renamed="results/intermediate_files/{name}_renamed_summary_stats.tsv"
    resources:
        time="1:00:00",
        mem_mb=500
    params:
        filetype=lambda wildcards: config["summary_stats"][wildcards.name]["type"]
    shell:
        """
        if [[ "{params.filetype}" == "recombination" ]]; then
            awk 'BEGIN {{ OFS="\\t"; print "SNP", "A2", "A2", "BETA", "SE", "P" }} NR > 1 {{ print $3, $6, $6, $10, $11, $13 }}' {input.summary_stats} > {output.summary_stats_renamed};
        elif [[ "{params.filetype}" == "aneuploidy" ]]; then
            zcat {input.summary_stats} | awk 'BEGIN {{ OFS="\\t"; print "SNP", "A1", "A2", "BETA", "SE", "P" }} {{ print $8, $13, $12, $3, $4, $6 }}' > {output.summary_stats_renamed};
        else
            cp {input.summary_stats} {output.summary_stats_renamed};
        fi
        """


rule cpra2rsid:
    """
    For aneuploidy and recombination summary stats, create RSID column.
    For public datasets, copy existing summary stats.
    """
    input:
        cpra2rsid_exec=config["cpra2rsid_exec"],
        summary_stats="results/intermediate_files/{name}_renamed_summary_stats.tsv", 
        dbsnp=config["dbsnp"],
        dictionary=rules.process_dbsnp.output.cpra2rsid_info
    output:
        summary_stats_cpra="results/intermediate_files/{name}_summary_stats_cpra.tsv"
    params:
        filetype=lambda wildcards: config["summary_stats"][wildcards.name]["type"],
    threads: 1
    resources:
        time="3:00:00",
        mem_mb=128000,
        disk_mb=200000
    shell:
        """
        if [[ "{params.filetype}" == "recombination" || "{params.filetype}" == "aneuploidy" ]]; then
            python3 {input.cpra2rsid_exec} --sumstats {input.summary_stats} --dbsnp {input.dbsnp} --dictionary {input.dictionary} --output {output.summary_stats_cpra};
        else
            cp {input.summary_stats} {output.summary_stats_cpra};
        fi
        """


# -------- Step 2: Process all summary stats for use in LDSC ------- #

rule munge_summary_stats: 
    """Apply the LDSC data-cleaning script to each summary statistic."""
    input:
        munge_exec=config["munge_exec"],
        summary_stats=rules.cpra2rsid.output.summary_stats_cpra,
        allele_list=config["allele_list"]
    output:
        summary_stats_munged="results/intermediate_files/{name}_munged.sumstats.gz",
        log="results/intermediate_files/{name}_munged.log"
    threads: 1
    resources:
        time="0:20:00",
        mem_mb=128000
    params:
        num_individuals=lambda wildcards: config["summary_stats"][wildcards.name]["N"],
        chunksize=50000,
        outfix="results/intermediate_files/{name}_munged"
    conda:
        "ldsc_env.yaml"
    shell: 
        """
        python2 {input.munge_exec} --sumstats {input.summary_stats} \
        --merge-alleles {input.allele_list} \
        --out {params.outfix} \
        --a1-inc --N {params.num_individuals} --chunksize {params.chunksize}
        """



# -------- Step 3: Calculate heritability on all summary stats ------- #

rule heritability:
    """Calculate heritability of each trait."""
    input:
        ldsc_exec=config["ldsc_exec"],
        summary_stats=expand(
            "results/intermediate_files/{name}_summary_stats_cpra.tsv",
            name=[
                key for key, value in config["summary_stats"].items() 
                if value["type"] not in {"aneuploidy", "recombination"}
            ]
        )
    output:
        heritability="results/heritability/{name}_heritability.log"
    resources:
        time="0:05:00",
        mem_mb=128000
    params:
        ld_scores=lambda wildcards: config["summary_stats"][wildcards.name]["ld_scores"],
        outfix="results/heritability/{name}_heritability"
    conda:
        "ldsc_env.yaml"
    shell:
        """
        python2 {input.ldsc_exec} \
        --h2 {input.summary_stats} \
        --ref-ld-chr {params.ld_scores} \
        --w-ld-chr {params.ld_scores} \
        --out {params.outfix}
        """


rule merge_heritability_results:
    """Merge heritability results into a single file."""
    input:
        logs=expand(
            "results/heritability/{name}_heritability.log", 
            name=[
                key for key, value in config["summary_stats"].items() 
                if value["type"] not in {"aneuploidy", "recombination"}
            ]
        )
    output:
        merged="results/heritability_published_merged.txt"
    resources:
        time="0:05:00",
        mem_mb=4000,
        disk_mb=5000
    run:
        import re

        # Define the regex patterns to extract relevant values
        patterns = {
            "total_h2": r"Total Observed scale h2:\s+([\d.]+)\s+\(([\d.]+)\)",
            "lambda_gc": r"Lambda GC:\s+([\d.]+)",
            "mean_chi2": r"Mean Chi\^2:\s+([\d.]+)",
            "intercept": r"Intercept:\s+([\d.]+)\s+\(([\d.]+)\)",
            "snps": r"After merging with regression SNP LD, (\d+) SNPs remain"
        }

        # Prepare to write the merged output
        with open(output.merged, "w") as outfile:
            # Write the header row
            outfile.write("Trait\tTotal_h2\tTotal_h2_SE\tLambda_GC\tMean_Chi2\tIntercept\tIntercept_SE\tSNPs\n")

            # Iterate through each input log file
            for log_file in input.logs:
                with open(log_file, "r") as infile:
                    content = infile.read()

                # Extract values using regex
                trait = log_file.split("/")[-1].replace("_heritability.log", "")
                total_h2, total_h2_se = re.search(patterns["total_h2"], content).groups()
                lambda_gc = re.search(patterns["lambda_gc"], content).group(1)
                mean_chi2 = re.search(patterns["mean_chi2"], content).group(1)
                intercept, intercept_se = re.search(patterns["intercept"], content).groups()
                snps = re.search(patterns["snps"], content).group(1)

                # Write extracted values to the output file
                outfile.write(f"{trait}\t{total_h2}\t{total_h2_se}\t{lambda_gc}\t{mean_chi2}\t{intercept}\t{intercept_se}\t{snps}\n")


# -------- Step 4: Calculate genetic correlation on relevant pairs of summary stats ------- #

from itertools import combinations

# Create pairings for each set of summary stats based on population 
def pairwise_comparisons(config, population_filter):
    # Filter traits based on the specified population
    traits = [
        trait for trait, details in config["summary_stats"].items()
        if details["population"] == population_filter
    ]
    
    # Generate all pairwise combinations
    pairwise = list(combinations(traits, 2))
    
    return pairwise

# Pairings between published summary stats and European-specific subsets of Natera
pairwise_european = pairwise_comparisons(config, "European")


rule pairwise_genetic_correlation:
	"""Calculate genetic correlation between each pairing of traits."""
	input:
		ldsc_exec=config["ldsc_exec"],
		trait1_file=lambda wildcards: f"results/intermediate_files/{wildcards.trait1}_munged.sumstats.gz",
		trait2_file=lambda wildcards: f"results/intermediate_files/{wildcards.trait2}_munged.sumstats.gz",
	output:
		genetic_correlation="results/genetic_correlation/{trait1}-{trait2}.log"
	resources:
		time="0:05:00",
		mem_mb=128000
	params:
		ld_scores=lambda wildcards: config["summary_stats"][wildcards.trait1]["ld_scores"],
		outfix="results/genetic_correlation/{trait1}-{trait2}"
	conda:
		"ldsc_env.yaml"
	shell:
		"""
		python2 {input.ldsc_exec} \
		--rg {input.trait1_file},{input.trait2_file} \
		--ref-ld-chr {params.ld_scores} \
		--w-ld-chr {params.ld_scores} \
		--out {params.outfix}
		"""


rule merge_genetic_correlation:
	"""Merge genetic correlation results for the relevant pairings."""
	input:
		expand("results/genetic_correlation/{trait1}-{trait2}.log",
			trait1=[t1 for t1, t2 in pairwise_european],
			trait2=[t2 for t1, t2 in pairwise_european])
	output:
		intermediate=temp("results/genetic_correlation_merged_space.txt"),
		merged="results/genetic_correlation_merged.txt"
	resources:
		time="0:05:00",
		mem_mb=4000,
		disk_mb=5000
	shell:
		"""
		echo "p1 p2 rg se z p h2_obs h2_obs_se h2_int h2_int_se gcov_int gcov_int_se" > {output.intermediate}
		awk '/^p1/{{flag=1; next}} flag && !/^Analysis finished|^Total time elapsed/{{print; flag=0}}' {input} >> {output.intermediate}
		awk -v OFS="\t" '$1=$1' {output.intermediate} > {output.merged}
		"""


# -------- Step 5: Conduct pheWAS on traits in this analysis.------- #
rule pheWAS:
    """Extract lines matching a given RSID from summary stats files and merge them into a single table."""
    input:
        summary_stats_cpra=expand("results/intermediate_files/{name}_summary_stats_cpra.tsv", name=config["summary_stats"].keys())
    output:
        merged_results="results/pheWAS_results_{rsid}.tsv"
    params:
        rsid=config["rsid"]
    shell:
        """
        # Create an empty file for the merged results
        echo -e "File\tExtracted_Line" > {output.merged_results}

        # Loop through all input files and extract the line containing the RSID
        for file in {input.summary_stats_cpra}; do
            # Extract the line containing the RSID
            result=$(grep "{params.rsid}" ${{file}})
            # If a result is found, append it to the output with the filename
            if [[ -n "$result" ]]; then
                echo -e "$(basename ${{file}})\t$result" >> {output.merged_results}
            fi
        done
        """


rule extract_snps:
    """Extract SNPs that were genome-wide significant from both aneuploidy and recombination traits."""
    input:
        filter_recomb_exec=config["filter_recomb_exec"],
        extractsnps_exec=config["extractsnps_exec"],
        lead_variants_recombination="/scratch16/rmccoy22/abiddan1/natera_recomb/analysis/gwas/results/gwas_output/regenie/finalized/natera_recombination_gwas.sumstats.replication.rsids.tsv",
        aneuploidy_summary_stats="/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/quantgen/results/intermediate_files/maternal_meiotic_aneuploidy_by_mother_summary_stats_cpra.tsv",
    output:
        filtered_recomb=temp("results/intermediate_files/recomb_hits_filtered.tsv"),
        significant_snps="results/intermediate_files/gw_significant_snps.txt"
    shell:
        """
        ml r/4.3.0
        Rscript {input.filter_recomb_exec} {input.lead_variants_recombination} {output.filtered_recomb}
        bash {input.extractsnps_exec} {input.aneuploidy_summary_stats} {output.filtered_recomb} {output.significant_snps}
        """


rule merge_summary_stats: 
    """Query SNPs that were significant in aneuploidy and recombination phenotypes across all traits."""
    input:
        merge_summary_stats=config["merge_summary_stats_exec"],
        significant_snps=rules.extract_snps.output.significant_snps,
        summary_stats_cpra=expand("results/intermediate_files/{name}_summary_stats_cpra.tsv", name=config["summary_stats"].keys()),
    output:
        snps_across_traits="results/queried_snps_across_traits.tsv"
    resources:
    	mem_mb=128000,
        disk_mb=128000,
    shell:
        """
        ml r/4.3.0
        Rscript --vanilla {input.merge_summary_stats} {input.significant_snps} {output.snps_across_traits} {input.summary_stats_cpra}
        """

