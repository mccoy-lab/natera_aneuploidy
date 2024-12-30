#!python3


configfile: "config.yaml"


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
		summary_stats_renamed="/results/intermediate_files/" # can be made temp 
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
		summary_stats_cpra_intermediate=
		summary_stats_cpra="/results/intermediate_files/"
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
		ld_score="/results/intermediate_files/" # once for each chromosome (maybe use as an outfix)
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


# -------- Step 3: Calculate heritability and genetic correlation on all summary stats ------- #

# run heritability on each summary stat 
# could make two different rules, one to use the calculated and one to use the downloaded 
# but i think i can define the ld scores in the config
rule heritability:
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


# run genetic correlation for aneuploidy and recombination 
rule gc_aneu_recomb: 
	# input: summary stats (output from munged)
	# input: calculated LD scores
	# output: one file for each run 
	shell:
		"""
		# aneu; rec1; rec2; rec3
		# rec1; rec2; rec3
		# rec2; rec3
		"""


# run genetic correlation for other traits 
# this will use the GWAS on European only 
rule gc_published: 
	# input: summary stats (output from munged)
	# input: downloaded LD scores 
	# output: one file for each run
	shell:
		"""
		# aneu; PL1; PL2; infert1
		# recomb; PL1; PL2; infert1
		# PL1; PL2; infert1
		# PL2; infert1
		"""








