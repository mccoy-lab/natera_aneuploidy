#!python3

# Usage: nohup snakemake -p --cores 48 --snakefile gwas.smk > nohup_date.out 2>&1 &
# Usage on rockfish: nohup snakemake -p --snakefile gwas.smk -j 200 --profile ~/code/rockfish_smk_profile/ &
# Optional: add -j 12 to submit as 12 jobs, etc.
# Optional: add -n for a dry run
# Executed from /scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/

# -------- Variables and filepaths for GWAS ------- #
#king_exec = "~/code/king"
general_outputs_fp = "results/"
vcf_fp = "/data/rmccoy22/natera_spectrum/genotypes/opticall_parents_031423/genotypes/"
metadata = (
    "/data/rmccoy22/natera_spectrum/data/summary_metadata/spectrum_metadata_merged.csv"
)
pcs_out = "results/parental_genotypes_pcs/"
ploidy_calls = "/data/rmccoy22/natera_spectrum/karyohmm_outputs/compiled_output/natera_embryos.karyohmm_v30a.bph_sph_trisomy.full_annotation.031624.tsv.gz"
gwas_results = "results/gwas/"
imputed_vcf_fp = "/data/rmccoy22/natera_spectrum/genotypes/imputed_parents_101823_cpra/"


# Dictionary of number of files to chunk each vcf into in `split`
# executed by /scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/scripts/count_and_split.sh
chunks_dict = {
    "chr1": 60,
    "chr2": 60,
    "chr3": 60,
    "chr4": 60,
    "chr5": 60,
    "chr6": 60,
    "chr7": 60,
    "chr8": 60,
    "chr9": 60,
    "chr10": 60,
    "chr11": 60,
    "chr12": 60,
    "chr13": 60,
    "chr14": 60,
    "chr15": 50,
    "chr16": 40,
    "chr17": 40,
    "chr18": 40,
    "chr19": 30,
    "chr20": 30,
    "chr21": 20,
    "chr22": 20,
    "chr23": 60,
}

# Parameters pipeline will run on 
phenotypes = ["embryo_count", "embryo_count_euploid", "maternal_age", "maternal_meiotic_aneuploidy", "haploidy", "triploidy"]
parents = ["mother", "father"]
dataset_type = ["discovery", "test"]

# shell.prefix("set -o pipefail; ")


# -------- Rules section -------- #
rule all:
    input:
        expand(
            gwas_results + "gwas_{phenotype}_by_{parent}_{dataset_type}_total.tsv.gz",
            phenotype=phenotypes,
            parent=parents,
            dataset_type=dataset_type,
        ),


# -------- 0. Preprocess Genetic Data -------- #
rule vcf2pgen: 
	"""Convert imputed VCF files to PGEN format for use in KING and PCA"""
	input: 
		input_vcf=imputed_vcf_fp + "spectrum_imputed_chr{chrom}_rehead_filter_cpra.vcf.gz",
	output: 
		pgen=temp("results/gwas/spectrum_imputed_chr{chrom}.pgen"),
		psam=temp("results/gwas/spectrum_imputed_chr{chrom}.psam"),
		pvar=temp("results/gwas/spectrum_imputed_chr{chrom}.pvar"),
	threads: 24
	params:
		outfix=lambda wildcards: gwas_results + f"spectrum_imputed_chr{wildcards.chrom}",
	shell:
		"plink2 --vcf {input.input_vcf} dosage=DS --double-id --maf 0.005 --threads {threads} --make-pgen --out {params.outfix}"


rule run_king:
    """Reformat parental genotypes vcf and run king to identify related individuals"""
    input:
        vcf_input=vcf_fp + "opticall_concat_total.norm.b38.vcf.gz",
    output:
        king_bed=temp(
            general_outputs_fp + "opticall_concat_total.norm.b38.alleleorder.bed"
        ),
        king_bim=temp(
            general_outputs_fp + "opticall_concat_total.norm.b38.alleleorder.bim"
        ),
        king_fam=temp(
            general_outputs_fp + "opticall_concat_total.norm.b38.alleleorder.fam"
        ),
        king_log=temp(
            general_outputs_fp + "opticall_concat_total.norm.b38.alleleorder.log"
        ),
    params:
        plink_outfix=general_outputs_fp + "opticall_concat_total.norm.b38.alleleorder",
        king_outfix=general_outputs_fp,
    shell:
        """
        plink --vcf {input.vcf_input} --double-id --allow-extra-chr --make-bed --out {params.plink_outfix} 
        {king_exec} -b {output.king_bed} --unrelated --degree 2 --prefix {params.king_outfix}
        """


rule run_plink_pca:
    """Run plink to get principal components for parental genotypes"""
    input:
        concat_vcf=vcf_fp + "opticall_concat_total.norm.b38.vcf.gz",
        concat_vcf_tbi=vcf_fp + "opticall_concat_total.norm.b38.vcf.gz.tbi",
    output:
        eigenvec=pcs_out + "parental_genotypes.eigenvec",
        eigenval=pcs_out + "parental_genotypes.eigenval",
        log=pcs_out + "parental_genotypes.log",
    resources:
        time="1:00:00",
        mem_mb="4G",
    params:
        outfix=pcs_out + "parental_genotypes",
        pcs=20,
    threads: 24
    shell:
        "plink2 --vcf {input.concat_vcf} --pca {params.pcs} approx --threads {threads} --out {params.outfix}"


rule discovery_validate_split:
    """Split families into discovery/validation sets for use in GWAS"""
    input:
        discovery_validate_R="scripts/discovery_test_split.R",
        metadata=metadata,
        fam_file="/data/rmccoy22/natera_spectrum/genotypes/opticall_parents_031423/genotypes/opticall_concat_total.norm.b38.fam",
        king_to_remove=general_outputs_fp + "unrelated_toberemoved.txt",
    output:
        discovery_validate_maternal=general_outputs_fp
        + "discover_validate_split_mother.txt",
        discovery_validate_paternal=general_outputs_fp
        + "discover_validate_split_father.txt",
    shell:
        "Rscript --vanilla {input.discovery_validate_R} {input.metadata} {input.fam_file} {input.king_to_remove} {output.discovery_validate_maternal} {output.discovery_validate_paternal}"


# -------- 1. Subset Genetic Data -------- #
rule get_chrom_pos:
    input:
        input_vcf=imputed_vcf_fp
        + "spectrum_imputed_chr{chrom}_rehead_filter_cpra.vcf.gz",
    output:
        chrom_mapfile=gwas_results
        + "subsets/spectrum_imputed_chr{chrom}_chrom_pos.txt",
    resources:
        mem_mb="2G",
    threads: 1
    shell:
        """
        bcftools query -f'%CHROM\t%POS\n' {input.input_vcf} > {output.chrom_mapfile}
        """


rule make_vcf_regions:
    input:
        chrom_mapfile=rules.get_chrom_pos.output.chrom_mapfile,
        get_regions="scripts/gwas/get_regions.py",
    output:
        regions_file=gwas_results
        + "subsets/spectrum_imputed_chr{chrom}_regions.txt",
    resources:
        mem_mb="2G",
    params:
        nchunks=lambda wildcards: chunks_dict[f"chr{wildcards.chrom}"],
    threads: 1
    shell:
        """
        python {input.get_regions} {params.nchunks} {input.chrom_mapfile} {output.regions_file}
        """


rule bed_split_vcf:
    input:
        regions_file=rules.make_vcf_regions.output.regions_file,
        input_vcf=imputed_vcf_fp
        + "spectrum_imputed_chr{chrom}_rehead_filter_cpra.vcf.gz",
    output:
        bcf=gwas_results
        + "subsets/spectrum_imputed_chr{chrom}_rehead_filter_cpra_{chunk}.bcf",
        bed=gwas_results
        + "subsets/spectrum_imputed_chr{chrom}_rehead_filter_cpra_{chunk}.bed",
        bim=gwas_results
        + "subsets/spectrum_imputed_chr{chrom}_rehead_filter_cpra_{chunk}.bim",
        fam=gwas_results
        + "subsets/spectrum_imputed_chr{chrom}_rehead_filter_cpra_{chunk}.fam",
        log=gwas_results
        + "subsets/spectrum_imputed_chr{chrom}_rehead_filter_cpra_{chunk}.log",
        nosex=gwas_results
        + "subsets/spectrum_imputed_chr{chrom}_rehead_filter_cpra_{chunk}.nosex",
    resources:
        mem_mb="6G",
    params:
        nchunks=lambda wildcards: chunks_dict[f"chr{wildcards.chrom}"],
        outfix=gwas_results
        + "subsets/spectrum_imputed_chr{chrom}_rehead_filter_cpra_{chunk}",
    threads: 16
    shell:
        """
        region=$(awk -v n={wildcards.chunk} "NR==n+1 {{print}}" {input.regions_file})
        bcftools view -r $region -Ob {input.input_vcf} > {output.bcf}
        plink --bcf {output.bcf} --double-id --allow-extra-chr --make-bed --out {params.outfix}
        """

# -------- 2. Create aneuploidy phenotypes -------- #
rule generate_aneuploidy_phenotypes:
    """Make file for each aneuploidy phenotype"""
    input:
        rscript="scripts/phenotypes/generate_phenotype_files.R",
        ploidy_calls=ploidy_calls,
        metadata=metadata,
    output:
        phenotype_file="results/phenotypes/{phenotype}_by_{parent}.csv",
    wildcard_constraints:
        parent="mother|father",
        phenotype="embryo_count|embryo_count_euploid|maternal_age|maternal_meiotic_aneuploidy|haploidy|triploidy",
    resources:
        time="0:30:00",
        mem_mb="10G",
    params:
        filter_day_5="TRUE",
        bayes_factor_cutoff=2,
        nullisomy_threshold=5,
        min_prob=0.9,
        max_meiotic=3,
        min_ploidy=15,
    shell:
        """
        ml gcc r/4.0.2
        Rscript --vanilla {input.rscript} {input.ploidy_calls} {wildcards.parent} {input.metadata} {wildcards.phenotype} {params.filter_day_5} {params.bayes_factor_cutoff} {params.nullisomy_threshold} {params.min_prob} {params.max_meiotic} {params.min_ploidy} {output.phenotype_file}
        """


# -------- 3. Execute GWAS and concatenate files -------- #
rule run_gwas_subset:
    """Run GWAS for each set of parameters, using the subsetted bed files"""
    input:
        gwas_rscript="scripts/gwas/gwas_all.R",
        metadata=metadata,
        bed=rules.bed_split_vcf.output.bed,
        discovery_test=general_outputs_fp + "discover_validate_split_{parent}.txt",
        parental_pcs=rules.run_plink_pca.output.eigenvec,
        phenotype_file=rules.generate_aneuploidy_phenotypes.output.phenotype_file,
        bim=rules.bed_split_vcf.output.bim,
    output:
        gwas_output=gwas_results
        + "subset_gwas_{phenotype}_by_{parent}_{dataset_type}_{chrom}_{chunk}.tsv",
    threads: 16
    resources:
        time="0:30:00",
        mem_mb="10G",
    wildcard_constraints:
        dataset_type="discovery|test",
        phenotype="embryo_count|embryo_count_euploid|maternal_age|maternal_meiotic_aneuploidy|haploidy|triploidy",
        parent="mother|father",
    shell:
        """
        ml gcc r/4.0.2
        Rscript --vanilla {input.gwas_rscript} {input.metadata} {input.bed} {input.discovery_test} {input.parental_pcs} {input.phenotype_file} {input.bim} {wildcards.dataset_type} {wildcards.phenotype} {wildcards.parent} {threads} {output.gwas_output}
        """


rule merge_subsets:
    """Create single file for GWAS for each chromosome, merging all subsets"""
    input:
        lambda wildcards: expand(
            gwas_results
            + "subset_gwas_{{phenotype}}_by_{{parent}}_{{dataset_type}}_{{chrom}}_{chunk}.tsv",
            phenotype=wildcards.phenotype,
            parent=wildcards.parent,
            dataset_type=wildcards.dataset_type,
            chrom=wildcards.chrom,
            chunk=range(chunks_dict.get(f"chr{wildcards.chrom}", 0)),
        ),
    output:
        gwas_output=gwas_results
        + "gwas_{phenotype}_by_{parent}_{dataset_type}_{chrom}.tsv",
    shell:
        "cat {input} > {output.gwas_output}"


rule merge_chroms:
    """Create single file for each phenotype/parent/dataset, merging all chromosomes"""
    input:
        expand(
            gwas_results
            + "gwas_{{phenotype}}_by_{{parent}}_{{dataset_type}}_{chrom}.tsv",
            chrom=range(1, 23),
        ),
    output:
        merged_file=gwas_results
        + "gwas_{phenotype}_by_{parent}_{dataset_type}_total.tsv.gz",
    shell:
        "cat {input} | gzip > {output.merged_file}"
