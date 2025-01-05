#!python3


configfile: "config.yaml"

# Create all heritability and genetic correlation results 
rule all:
    input:
        "results/test/genetic_correlation_merged.txt"


# -------- Step 1: Steps to standardize Natera summary stats and supporting files for use in LDSC ------- #

rule process_dbsnp: 
	"""Generate cpra and rsid table for dbsnp."""
	input: 
		dbsnp=config["dbsnp"],
	output:
		cpra2rsid_info="/results/intermediate_files/dbsnp151_hg38_info.txt"
	threads: 16
	shell: 
		'bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t%ID\n" {input.dbsnp} | bgzip -@ {params.threads} > {output.cpra2rsid_info}'


rule rename_summary_stats: 
	"""For aneuploidy and recombination summary stats, downselect and reorder columns. For public datasets, copy existing summary stats."""
	input:
		summary_stats=lambda wildcards: config["summary_stats"][wildcards.name]["file"]
	outputs:
		summary_stats_renamed=temp(lambda wildcards: f"/results/intermediate_files/{wildcards.name}_renamed_summary_stats.tsv")
	params:
		filetype=lambda wildcards: config["summary_stats"][wildcards.name]["type"],
	run:
        if params.filetype == "recombination":
            shell(
                awk 'BEGIN { OFS="\t"; print "SNP", "CHR", "BP", "A1", "BETA", "SE", "P" } NR > 1 { print $3, $1, $2, $6, $10, $11, $13 }' {input.summary_stats} > {output.summary_stats_renamed}
            )
        elif params.filetype == "aneuploidy":
            shell(
                zcat {input.summary_stats} | awk 'BEGIN { OFS="\t"; print "SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P" } { print $8, $10, $1, $13, $12, $3, $4, $6 }' > {output.summary_stats_renamed}
            )
        else:
            shell('cp {input.summary_stats} {output.summary_stats_renamed}')


rule cpra2rsid: 
	"""For aneuploidy and recombination summary stats, add RSID column. For public datasets, copy existing summary stats."""
	input:
		cpra2rsid_exec=config["cpra2rsid_exec"],
		summary_stats=rules.rename_summary_stats.output.summary_stats_renamed,
		dbsnp=config["dbsnp"],
		dictionary=rules.process_dbsnp.output.cpra2rsid_info,
	output:
		summary_stats_cpra_intermediate=temp(lambda wildcards: f"/results/intermediate_files/{wildcards.name}_summary_stats_cpra_intermediate.tsv")
		summary_stats_cpra=lambda wildcards: f"/results/intermediate_files/{wildcards.name}_summary_stats_cpra.tsv"
	run:
		if params.filetype in {"recombination", "aneuploidy"}:
			shell('python3 {input.cpra2rsid_exec} --sumstats {input.summary_stats} --dbsnp {input.dbsnp} --dictionary {input.dictionary} --output {output.summary_stats_cpra_intermediate}')
			shell('awk 'NR==1 {gsub(/^SNP/, "CPRA"); gsub(/RSID$/, "SNP")} {print}' {output.summary_stats_cpra_intermediate} > {output.summary_stats_cpra}')
			if params.filetype == "recombination":
				shell('awk 'BEGIN {OFS="\t"} NR==1 {print $1, $2, $3, $4, "A2", $5, $6, $7, $8; next} {split($1, a, ":"); print $1, $2, $3, $4, a[3], $5, $6, $7, $8}' {output.summary_stats_cpra} > {output.summary_stats_cpra_intermediate}') 
				shell('cp {output.summary_stats_cpra_intermediate} {output.summary_stats_cpra}')
        else:
            shell('cp {input.summary_stats} {output.summary_stats_cpra}')


rule create_ldscores: 
	"""Calculate LD Scores for the Natera dataset for each chromosome."""
	input:
		ldsc_exec=config["ldsc_exec"],
		imputed_parents=config["imputed_parents"],
	output:
		ld_score=lambda wildcards: f"/results/ld_scores/LDscore.{wildcards.chrom}"
	# output: LD scores for all chromosomes 
	params:
		window=1000
	shell:
		"""
		conda activate ldsc
		python2 {input.ldsc_exec} --bfile {input.imputed_parents} --l2 --ld-wind-kb {params.window} --out {output.ld_score}
		conda deactivate
		"""


# -------- Step 2: Process all summary stats for use in LDSC ------- #

# uses conda 
rule munge_summary_stats: 
	"""Apply the LDSC data-cleaning script to each summary statistic."""
	input:
		munge_exec=config["munge_exec"],
		summary_stats=rules.cpra2rsid.output.summary_stats_cpra,
		allele_list=config["allele_list"],
	output:
		summary_stats_munged="/results/intermediate_files/"
	resources:
		time="0:20:00",
		mem_mb=
	params:
		num_individuals=lambda wildcards: config["summary_stats"][wildcards.name]["N"],
		chunksize=50000,
	shell: 
		"""
		conda activate ldsc
		python2 {input.munge_exec} --sumstats {input.summary_stats} --merge-alleles {input.allele_list} --out {output.summary_stats_munged} --a1-inc --N {params.num_individuals} --chunksize {params.chunksize}
		conda deactivate
		"""


# -------- Step 3: Calculate heritability on all summary stats ------- #

rule heritability:
	"""Calculate heritability of each trait."""
	input:
		ldsc_exec=config["ldsc_exec"],
		summary_stats=rules.munge_summary_stats.output.summary_stats_munged,
	output:
		heritability="/results/"
	params:
		ld_scores=lambda wildcards: config["summary_stats"][wildcards.name]["ld_scores"],
	shell:
		"""
		python2 {input.ldsc_exec} \
		--h2 munged_sumstats.gz \
		--ref-ld-chr {params.ld_scores}
		--w-ld-chr {params.ld_scores}
		--out {output.heritability}
		"""


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
    
    # Skip pairings where both traits have population "European"
    # and both types are "aneuploidy" or "recombination"
    if population_filter == "European":
        pairwise = [
            (t1, t2) for t1, t2 in pairwise
            if not (
                config["summary_stats"][t1]["type"] in {"aneuploidy", "recombination"} and
                config["summary_stats"][t2]["type"] in {"aneuploidy", "recombination"}
            )
        ]
    return pairwise


# Generate pairwise comparisons for each population
# Natera summary stats (population all)
pairwise_all = pairwise_comparisons(config, "All")
# Pairings between published and with European-specific subsets of Natera 
pairwise_european = pairwise_comparisons(config, "European")


rule pairwise_genetic_correlation:
    input:
        trait1_file=lambda wildcards: config["summary_stats"][wildcards.trait1]["file"],
        trait2_file=lambda wildcards: config["summary_stats"][wildcards.trait2]["file"],
        ld_scores=lambda wildcards: config["summary_stats"][wildcards.trait1]["ld_scores"]
    output:
        "results/test/{trait1}-{trait2}.txt"
    shell:
        "cat {input.trait1_file} {input.trait2_file} {input.ld_scores} > {output}"


rule merge_genetic_correlation:
    input:
        # Combine outputs from the pairwise comparisons
        expand("results/test/{trait1}-{trait2}.txt", 
               trait1=[t1 for t1, t2 in pairwise_all], 
               trait2=[t2 for t1, t2 in pairwise_all]) +
        expand("results/test/{trait1}-{trait2}.txt", 
               trait1=[t1 for t1, t2 in pairwise_european], 
               trait2=[t2 for t1, t2 in pairwise_european])
    output:
        "results/test/genetic_correlation_merged.txt"
    shell:
        "cat {input} > {output}"









