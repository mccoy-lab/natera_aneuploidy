#!python3

# =================
# author: Sara A. Carioscia, Biology Dept., Johns Hopkins University
# email: scarios1@jhu.edu
# last updated: February 10, 2025
# aim: Compute heritability and genetic correlation for recombination, aneuploidy, and
#       published fertility-related traits. Process input files as necessary.
# =================

# Usage: snakemake --snakefile quantgen.smk --use-conda --profile ~/code/rockfish_smk_profile -p


configfile: "config.yaml"


chromosomes = [str(i) for i in range(1, 24)]


# Create all heritability and genetic correlation results
rule all:
    input:
        #"results/heritability_published_merged.txt",
        #"results/intermediate_files/CentromereDist_Female_eur_summary_stats_cpra.tsv"
        #"results/genetic_correlation_merged.txt",
        #expand("results/pheWAS_results_{rsid}.tsv", rsid=config["rsid"]),
        #"results/queried_snps_across_traits.tsv",
        #"results/intermediate_files/MeanCO_Female_munged.sumstats.gz",
        #"results/intermediate_files/maternal_meiotic_aneuploidy_by_mother_munged.sumstats.gz",
        #"results/intermediate_files/age_at_menarche_reproGen_munged.sumstats.gz",
        #expand("/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/quantgen/results/ld_scores/filtered_natera_vcf/plink_files/spectrum_imputed_chr{chrom}_rehead_filterDR29_plink.bed", chrom=chromosomes),
        #"results/intermediate_files/ld_scores_EUR/hgdp1kgp_chr21.filtered.SNV_INDEL.phased.shapeit5.european_only.bed",
        expand("results/ld_scores_EUR_rsid/LDscore.19.l2.ldscore.gz", chrom=chromosomes)


# -------- Step 1: Steps to standardize Natera summary stats and supporting files for use in LDSC ------- #

rule process_dbsnp:
    """Generate cpra and rsid table for dbsnp."""
    input:
        dbsnp=config["dbsnp"],
    output:
        cpra2rsid_info="results/intermediate_files/dbsnp151_hg38_info.txt",
    resources:
        time="6:00:00",
        mem_mb=1000,
    threads: 16
    shell:
        'bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t%ID\n" {input.dbsnp} | bgzip -@ {threads} > {output.cpra2rsid_info}'


rule rename_summary_stats:
    """
    For aneuploidy and recombination summary stats, downselect and reorder columns.
    For public datasets, copy existing summary stats.
    """
    input:
        summary_stats=lambda wildcards: config["summary_stats"][wildcards.name]["file"],
    output:
        summary_stats_renamed="results/intermediate_files/{name}_renamed_summary_stats.tsv",
    resources:
        time="10:00",
        mem_mb=500,
    params:
        filetype=lambda wildcards: config["summary_stats"][wildcards.name]["type"],
    shell:
        """
        if [[ "{params.filetype}" == "recombination" ]]; then
            awk 'BEGIN {{ OFS="\\t"; print "SNP", "A1", "A2", "BETA", "SE", "P" }} NR > 1 {{ print $3, $6, $4, $10, $11, $13 }}' {input.summary_stats} > {output.summary_stats_renamed};
        elif [[ "{params.filetype}" == "aneuploidy" ]]; then
            tmp1=$(mktemp)
            zcat {input.summary_stats} | awk -F'\t' '{{ split($8, arr, ":"); \
                split(arr[3], alleles, ":"); new_col = (alleles[1] == $9) ? alleles[2] : alleles[1]; 
                print $0 "\t" new_col; }}' | gzip  > "$tmp1"
            zcat "$tmp1" | awk 'BEGIN {{ OFS="\\t"; print "SNP", "A1", "A2", "BETA", "SE", "P" }} {{ print $8, $9, $14, $3, $4, $6 }}' > {output.summary_stats_renamed};
        else
            cp {input.summary_stats} {output.summary_stats_renamed};
        fi
        """


rule cpra2rsid:
    """
    For European subsets of aneuploidy and recombination summary stats, create RSID column.
    For public datasets, and for full datasets of aneuploidy/recombination,
    copy existing summary stats.
    """
    input:
        cpra2rsid_exec=config["cpra2rsid_exec"],
        dbsnp=rules.process_dbsnp.output.cpra2rsid_info,
        summary_stats="results/intermediate_files/{name}_renamed_summary_stats.tsv",
    output:
        summary_stats_cpra="results/intermediate_files/{name}_summary_stats_cpra.tsv",
    params:
        traittype=lambda wildcards: config["summary_stats"][wildcards.name]["type"],
        population=lambda wildcards: config["summary_stats"][wildcards.name]["population"],
        filetype="summary_stats"
    threads: 1
    resources:
        time="1:30:00",
        mem_mb="132G",
    shell:
        """
        if [[ "{params.traittype}" == "recombination" || "{params.traittype}" == "aneuploidy" ]]; then
            tmp1=$(mktemp)
            tmp2=$(mktemp)

            python3 {input.cpra2rsid_exec} {input.dbsnp} {input.summary_stats} "$tmp1" {params.filetype};
            cut --complement -f 2 "$tmp1" > "$tmp2"
            sed -e '1s/NA/SNP/' "$tmp2" > {output.summary_stats_cpra}

            rm "$tmp1" "$tmp2"
        else
            cp {input.summary_stats} {output.summary_stats_cpra};
        fi
        """


# -------- Step 2: Process all summary stats for use in LDSC ------- #

rule munge_summary_stats:
    """Apply the LDSC data-cleaning script to each summary statistic."""
    input:
        munge_exec=config["munge_exec"],
        #summary_stats=rules.cpra2rsid.output.summary_stats_cpra,
        summary_stats="results/intermediate_files/{name}_summary_stats_cpra.tsv",
    output:
        summary_stats_munged="results/intermediate_files/{name}_munged.sumstats.gz",
        log="results/intermediate_files/{name}_munged.log",
    threads: 1
    resources:
        time="0:20:00",
        mem_mb=128000,
    params:
        num_individuals=lambda wildcards: config["summary_stats"][wildcards.name]["N"],
        chunksize=50000,
        outfix="results/intermediate_files/{name}_munged",
    conda:
        "ldsc_env.yaml"
    shell:
        """
        python2 {input.munge_exec} --sumstats {input.summary_stats} \
        --out {params.outfix} --N {params.num_individuals} \
        --chunksize {params.chunksize}
        """


# -------- Step 3: Calculate LD scores for Natera data (CPRA) ------- #

rule filter_Natera_vcf: 
	"""Filter imputed VCF to keep only sites with DR2>0.9 for use in LD scores."""
	input:
        vcf=config["imputed_parents_vcf"],
        vcfidx=config["imputed_parents_vcfidx"],
	output:
		filtered_vcf="results/ld_scores/filtered_natera_vcf/spectrum_imputed_chr{chrom}_rehead_filterDR29.cpra.vcf.gz",
		filtered_vcf_tabix="results/ld_scores/filtered_natera_vcf/spectrum_imputed_chr{chrom}_rehead_filterDR29.cpra.vcf.gz.tbi",
	resources:
		mem_mb="10G",
		time="2:00:00",
	threads: 16
	wildcard_constraints:
		chrom="|".join([str(i) for i in range(1, 24)]),
	params:
		dosage_r2_thresh=0.9
	conda:
		"geno.yaml"
	shell:
		"""
		bcftools view -i \'DR2>{params.dosage_r2_thresh}\' --threads {threads} {input.vcf} -Oz -o {output.filtered_vcf}; tabix -f {output.filtered_vcf}
		"""

rule vcf2bed: 
	"""Create plink output files for use in LD score calculation."""
	input:
		filtered_vcf="results/ld_scores/filtered_natera_vcf/spectrum_imputed_chr{chrom}_rehead_filterDR29.cpra.vcf.gz",
	output:
		bed="results/ld_scores/filtered_natera_vcf/plink_files/spectrum_imputed_chr{chrom}_rehead_filterDR29_plink.bed",
		bim="results/ld_scores/filtered_natera_vcf/plink_files/spectrum_imputed_chr{chrom}_rehead_filterDR29_plink.bim",
		fam="results/ld_scores/filtered_natera_vcf/plink_files/spectrum_imputed_chr{chrom}_rehead_filterDR29_plink.fam",
		log="results/ld_scores/filtered_natera_vcf/plink_files/spectrum_imputed_chr{chrom}_rehead_filterDR29_plink.log",
	resources:
		mem_mb="10G",
		time="45:00"
	params:
		outfix="results/ld_scores/filtered_natera_vcf/plink_files/spectrum_imputed_chr{chrom}_rehead_filterDR29_plink"
	conda:
		"geno.yaml"
	shell:
		"""
		plink --vcf {input.filtered_vcf} --memory 9000 --double-id --make-bed --out {params.outfix}
		"""


rule create_ldscores: 
    """Calculate LD Scores for the Natera dataset for each chromosome."""
    input:
        ldsc_exec=config["ldsc_exec"],
    output:
        ld_score_gz="results/ld_scores/LDscore.{chrom}.l2.ldscore.gz",
    resources:
        time="6:00:00",
        mem_mb=128000,
        disk_mb=128000
    params:
        outfix="results/ld_scores/LDscore.{chrom}",
        imputed_parents_prefix=lambda wildcards: config["imputed_parents_template_filtered"].format(chrom=wildcards.chrom),
        window=300,
        maf=0.005
    conda:
        "ldsc_env.yaml"
    shell:
        """
        python2 {input.ldsc_exec} --out {params.outfix} --bfile {params.imputed_parents_prefix} --l2 --ld-wind-kb {params.window} --maf {params.maf}
        """


# -------- Step 5: Calculate LD Scores for European individuals (RSID) ------- #

rule filter_EUR_individuals:
	"""Create vcf file for just EUR individuals for each chr."""
	input:
		metadata="/scratch4/rmccoy22/sharedData/populationDatasets/GnomAD_Genomes_HGDP_TGP/gnomad_meta_updated.tsv",
	output:
		samples="results/intermediate_files/ld_scores_EUR/eur.samples.txt",
	resources:
		mem_mb=500,
		time="30:00"
	shell:
		"""
		grep "EUR" {input.metadata} | awk '{{print $2}}' > {output.samples}
		"""

rule bcf2bed_hgdp1kgp:
	"""Create plink output files for use in calculating LD scores."""
	input:
		samples=rules.filter_EUR_individuals.output.samples,
		#bcf="/scratch4/rmccoy22/sharedData/populationDatasets/GnomAD_Genomes_HGDP_TGP/hgdp1kgp_chr{chrom}.filtered.SNV_INDEL.phased.shapeit5.bcf",
		bcf=lambda wildcards: f"/scratch4/rmccoy22/sharedData/populationDatasets/GnomAD_Genomes_HGDP_TGP/hgdp1kgp_chr{wildcards.chrom}.filtered.SNV_INDEL.phased.shapeit5.bcf"
			if wildcards.chrom in [str(c) for c in range(1, 23)] else
			"/scratch4/rmccoy22/sharedData/populationDatasets/1KGP_NYGC/GRCh38_phased_vcfs/1kGP_high_coverage_Illumina.chr23.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
	output:
		vcf="results/intermediate_files/ld_scores_EUR/hgdp1kgp_chr{chrom}.filtered.SNV_INDEL.phased.shapeit5.european_only.vcf.gz",
		bed="results/intermediate_files/ld_scores_EUR/hgdp1kgp_chr{chrom}.filtered.SNV_INDEL.phased.shapeit5.european_only.bed",
		bim="results/intermediate_files/ld_scores_EUR/hgdp1kgp_chr{chrom}.filtered.SNV_INDEL.phased.shapeit5.european_only.bim",
		fam="results/intermediate_files/ld_scores_EUR/hgdp1kgp_chr{chrom}.filtered.SNV_INDEL.phased.shapeit5.european_only.fam",
		log="results/intermediate_files/ld_scores_EUR/hgdp1kgp_chr{chrom}.filtered.SNV_INDEL.phased.shapeit5.european_only.log",
	threads: 8
	resources:
		mem_mb="10G",
		time="30:00"
	params:
		outfix="results/intermediate_files/ld_scores_EUR/hgdp1kgp_chr{chrom}.filtered.SNV_INDEL.phased.shapeit5.european_only",
	shell:
		"""
		# If chromosome is 23, replace chrX with chr23 in the VCF file 
		if [ "{wildcards.chrom}" == "23" ]; then
            bcftools view -S {input.samples} --force-samples -m2 -M2 -c 1 -q 0.005:minor {input.bcf} | \
            sed 's/^chrX$/chr23/' | \
            bgzip -@ {threads} > {output.vcf}
        else
            bcftools view -S {input.samples} --force-samples -m2 -M2 -c 1 -q 0.005:minor {input.bcf} | \
            bgzip -@ {threads} > {output.vcf}
        fi
		
		plink --vcf {output.vcf} --memory 9000 --double-id --make-bed --out {params.outfix}
		"""

rule create_ldscores_EUR: 
    """Calculate LD Scores for European subset of the 1kgphgp dataset for each chromosome."""
    input:
        ldsc_exec=config["ldsc_exec"],
        bed="results/intermediate_files/ld_scores_EUR/hgdp1kgp_chr{chrom}.filtered.SNV_INDEL.phased.shapeit5.european_only.bed",
    output:
        ld_score="results/ld_scores_EUR/LDscore.{chrom}.l2.ldscore.gz",
    resources:
        time="30:00",
        mem_mb=128000,
        disk_mb=128000
    params:
        outfix="results/ld_scores_EUR/LDscore.{chrom}",
        hgdp1kgp_prefix=lambda wildcards: config["hgdp1kgp_template"].format(chrom=wildcards.chrom),
        window=300,
        maf=0.005
    conda:
        "ldsc_env.yaml"
    shell:
        """
        python2 {input.ldsc_exec} --out {params.outfix} --bfile {params.hgdp1kgp_prefix} --l2 --ld-wind-kb {params.window} --maf {params.maf}
        """

rule cpra2rsid_ldscores_EUR:
    """Convert EUR ld scores from CPRA to RSID."""
    input:
        cpra2rsid_exec=config["cpra2rsid_exec"],
        dbsnp=rules.process_dbsnp.output.cpra2rsid_info,
        ld_score_in="results/ld_scores_EUR/LDscore.{chrom}.l2.ldscore.gz",
        ld_score_M_in="results/ld_scores_EUR/LDscore.{chrom}.l2.M",
        ld_score_M_5_50_in="results/ld_scores_EUR/LDscore.{chrom}.l2.M_5_50",
    output:
        ld_score_temp=temp("results/intermediate_files/LDscore.{chrom}.temp.l2.ldscore.gz"),
        ld_score_gz="results/ld_scores_EUR_rsid/LDscore.{chrom}.l2.ldscore.gz",
        ld_score_M="results/ld_scores_EUR_rsid/LDscore.{chrom}.l2.M",
        ld_score_M_5_50="results/ld_scores_EUR_rsid/LDscore.{chrom}.l2.M_5_50",
    resources:
        time="1:30:00",
        mem_mb="132G",
    params:
        filetype="summary_stats"
    shell:
        """
        python3 {input.cpra2rsid_exec} {input.dbsnp} {input.ld_score_in} {output.ld_score_temp} {params.filetype}
        awk 'BEGIN{{OFS=FS="\t"}} NR==1 {{$1="SNP"; $3="CPRA"}} 1' {output.ld_score_temp}| bgzip > {output.ld_score_gz}
        # Copy SNP count files to renamed directory for use with RSID ld scores
        cp {input.ld_score_M_in} {output.ld_score_M}
        cp {input.ld_score_M_5_50_in} {output.ld_score_M_5_50}
        """


# -------- Step 6: Calculate heritability on all summary stats ------- #

rule heritability:
    """Calculate heritability of each trait."""
    input:
        ldsc_exec=config["ldsc_exec"],
        summary_stats="results/intermediate_files/{name}_munged.sumstats.gz",
    output:
        heritability="results/heritability/{name}_heritability.log",
    resources:
        time="0:05:00",
        mem_mb=128000,
    params:
        ld_scores=lambda wildcards: config["summary_stats"][wildcards.name]["ld_scores"],
        outfix="results/heritability/{name}_heritability",
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
				if (
                    value["type"] not in {"recombination", "aneuploidy"}  # include all published 
                    # include recombination and aneuploidy only for their full datsets
                    or (value["type"] in {"recombination", "aneuploidy"} and value["population"] == "All") 
                )
			]
		),
	output:
		merged="results/heritability_published_merged.txt",
	resources:
		time="0:05:00",
		mem_mb=4000,
		disk_mb=5000,
	run:
		import re

		# Define the regex patterns to extract relevant values
		patterns = {
			"total_h2": r"Total Observed scale h2:\s+([\d.]+)\s+\(([\d.]+)\)",
			"lambda_gc": r"Lambda GC:\s+([\d.]+)",
			"mean_chi2": r"Mean Chi\^2:\s+([\d.]+)",
			"intercept": r"Intercept:\s+([\d.]+)\s+\(([\d.]+)\)",
			"snps": r"After merging with regression SNP LD, (\d+) SNPs remain",
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


# -------- Step 7: Calculate genetic correlation on relevant pairs of summary stats ------- #

from itertools import combinations


# Create pairings for each set of summary stats based on population
def pairwise_comparisons(config, population_filter):
    # Filter traits based on the specified population
    traits = [
        trait
        for trait, details in config["summary_stats"].items()
        if details["population"] == population_filter
    ]

    # Generate all pairwise combinations
    pairwise = list(combinations(traits, 2))

    return pairwise


# Pairings between full summary stats for within recombination/aneuploidy comparisions
pairwise_all = pairwise_comparisons(config, "All")
# Pairings between published summary stats and European-specific subsets of Natera
pairwise_european = pairwise_comparisons(config, "European")


rule pairwise_genetic_correlation:
    """Calculate genetic correlation between each pairing of traits."""
    input:
        ldsc_exec=config["ldsc_exec"],
        trait1_file=lambda wildcards: f"results/intermediate_files/{wildcards.trait1}_munged.sumstats.gz",
        trait2_file=lambda wildcards: f"results/intermediate_files/{wildcards.trait2}_munged.sumstats.gz",
    output:
        genetic_correlation="results/genetic_correlation/{trait1}-{trait2}.log",
    resources:
        time="0:05:00",
        mem_mb=128000,
    params:
        ld_scores=lambda wildcards: config["summary_stats"][wildcards.trait1][
            "ld_scores"
        ],
        outfix="results/genetic_correlation/{trait1}-{trait2}",
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


pairwise_inputs = sorted(
    set(
        f"results/genetic_correlation/{min(t1, t2)}-{max(t1, t2)}.log"
        for subset in [pairwise_all, pairwise_european]
        for t1, t2 in subset
    )
)
# Split into groups
n_tranches = 20
tranche_size = (len(pairwise_inputs) + n_tranches - 1) // n_tranches
tranches = [pairwise_inputs[i:i + tranche_size] for i in range(0, len(pairwise_inputs), tranche_size)]

rule merge_genetic_correlation_tranche:
    """Merge genetic correlation results for the relevant pairings."""
    input:
        lambda wildcards: tranches[int(wildcards.tranche)]
    output:
        intermediate=temp("results/genetic_correlation_tranche_{tranche}_space.txt"),
        merged=temp("results/genetic_correlation_tranche_{tranche}.txt"),
    resources:
        time="0:05:00",
        mem_mb=4000,
        disk_mb=5000,
    shell:
        """
        echo "p1 p2 rg se z p h2_obs h2_obs_se h2_int h2_int_se gcov_int gcov_int_se" > {output.intermediate}
        awk '/^p1/{{flag=1; next}} flag && !/^Analysis finished|^Total time elapsed/{{print; flag=0}}' {input} >> {output.intermediate}
        awk -v OFS="\t" '$1=$1' {output.intermediate} > {output.merged}
        """


rule merge_genetic_correlation_final:
    """Concatenate and finalize genetic correlation results."""
    input:
        expand("results/genetic_correlation_tranche_{tranche}.txt", tranche=range(n_tranches))
    output:
        "results/genetic_correlation_merged.txt",
    shell:
        """
        # Extract header from the first tranche
        head -n 1 {input[0]} > {output}
        # Concatenate all files, skipping subsequent headers
        tail -n +2 -q {input} >> {output}
        """


# -------- Step 8: Conduct pheWAS on traits in this analysis.------- #

rule pheWAS:
    """Extract lines matching a given RSID from summary stats files and merge them into a single table."""
    input:
        summary_stats_cpra=expand(
            "results/intermediate_files/{name}_summary_stats_cpra.tsv",
            name=config["summary_stats"].keys(),
        ),
    output:
        merged_results=expand("results/pheWAS_results_{rsid}.tsv", rsid=config["rsid"]),
    params:
        rsid=config["rsid"],  # List of RSIDs from the config
    shell:
        """
        # Loop over each RSID and create the file for each
        for rsid in {params.rsid}; do
            # Create an empty file for the merged results
            echo -e "File\tSNP\tA1\tA2\tBeta\tSE\tP" > results/pheWAS_results_${{rsid}}.tsv

            # Loop through all input files and extract the line containing the RSID
            for file in {input.summary_stats_cpra}; do
                # Extract the line containing the RSID
                result=$(grep "${{rsid}}" ${{file}})
                # If a result is found, append it to the output with the filename
                if [[ -n "$result" ]]; then
                    # Remove the suffix "_summary_stats_cpra.tsv" from the filename
                    modified_name=$(basename ${{file}} | sed 's/_summary_stats_cpra.tsv$//')
                    echo -e "$modified_name\t$result" >> results/pheWAS_results_${{rsid}}.tsv
                fi
            done

            # Replace all white space with tabs
            sed -E -i 's/[[:space:]]+/\\t/g' results/pheWAS_results_${{rsid}}.tsv
        done
        """


rule extract_snps:
    """Extract SNPs that were genome-wide significant from both aneuploidy and recombination traits."""
    input:
        filter_recomb_exec=config["filter_recomb_exec"],
        extractsnps_exec=config["extractsnps_exec"],
        lead_variants_recombination=config["lead_variants_recombination"],
        aneuploidy_summary_stats=config["aneuploidy_summary_stats"]
    output:
        filtered_recomb=temp("results/intermediate_files/recomb_hits_filtered.tsv"),
        significant_snps="results/intermediate_files/gw_significant_snps.txt",
    shell:
        """
        ml r/4.3.0
        Rscript {input.filter_recomb_exec} {input.lead_variants_recombination} {output.filtered_recomb}
        bash {input.extractsnps_exec} {input.aneuploidy_summary_stats} {output.filtered_recomb} {output.significant_snps}
        """


# Filter to query only the aneuploidy and recombination summary stats based on the full dataset,
# not the European subsets
def filter_summary_stats(config):
    filtered_files = []
    for name, data in config["summary_stats"].items():
        if data["type"] in ["aneuploidy", "recombination"]:
            if data["population"].lower() != "european":  # Include only non-European
                filtered_files.append(
                    f"results/intermediate_files/{name}_summary_stats_cpra.tsv"
                )
        else:
            # Include any entries for all other phenotype types
            filtered_files.append(
                f"results/intermediate_files/{name}_summary_stats_cpra.tsv"
            )
    return filtered_files


rule merge_summary_stats:
    """Query SNPs that were significant in aneuploidy and recombination phenotypes across all traits."""
    input:
        merge_summary_stats=config["merge_summary_stats_exec"],
        significant_snps=rules.extract_snps.output.significant_snps,
        summary_stats_cpra=filter_summary_stats(config),
    output:
        snps_across_traits="results/queried_snps_across_traits.tsv",
    resources:
        mem_mb=128000,
        disk_mb=128000,
    shell:
        """
        ml r/4.3.0
        Rscript --vanilla {input.merge_summary_stats} {input.significant_snps} {output.snps_across_traits} {input.summary_stats_cpra}
        """
