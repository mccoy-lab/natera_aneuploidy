#!python3

# Usage: nohup snakemake -p --cores 48 --snakefile gwas.smk > nohup_date.out 2>&1 &
# Optional: add -j 12 to submit as 12 jobs, etc.
# Optional: add -n to do a dry run
# Executed from /scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/

# -------- Setting variables and paths for pre-GWAS processing steps and GWAS outputs ------- #
king_exec = "~/code/king"
general_outputs_fp = "results/"
vcf_fp = "/data/rmccoy22/natera_spectrum/genotypes/opticall_parents_031423/genotypes/"
metadata = (
    "/data/rmccoy22/natera_spectrum/data/summary_metadata/spectrum_metadata_merged.csv"
)
pcs_out = "results/parental_genotypes_pcs/"
gwas_results = "results/gwas/"


# Define the parameters that the pipeline will run on 
#chroms = range(1, 24)
phenotypes = ["embryo_count", "haploidy", "maternal_meiotic_aneuploidy", "triploidy", "parental_triploidy"]
parents = ["mother", "father"]
dataset_type = ["discovery", "test"]

shell.prefix("set -o pipefail; ")

# -------- Rule all to run whole pipeline -------- #
rule all:
    input:
        expand(
            gwas_results + "gwas_{phenotype}_by_{parent}_{dataset_type}_total.txt.gz",
            #phenotype=["embryo_count", "haploidy", "maternal_meiotic_aneuploidy", "triploidy", "parental_triploidy"],
            phenotype="maternal_meiotic_aneuploidy",
            parent="mother",
            dataset_type="discovery",           
        ),


# -------- Functions to determine run GWAS on maternal and paternal for each phenotype -------- #
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


rule vcf2bed:
    """Take vcf for each chromosome, convert to bed as a temp file for use in the GWAS"""
    input:
        vcf_input=vcf_fp + "opticall_concat_{chrom}.norm.b38.vcf.gz",
    output:
        bedfile=temp(gwas_results + "opticall_concat_{chrom}.norm.b38.bed"),
        bimfile=temp(gwas_results + "opticall_concat_{chrom}.norm.b38.bim"),
        famfile=temp(gwas_results + "opticall_concat_{chrom}.norm.b38.fam"),
        logfile=temp(gwas_results + "opticall_concat_{chrom}.norm.b38.log"),
    params:
        outfix=lambda wildcards: gwas_results
        + f"opticall_concat_{wildcards.chrom}.norm.b38",
    threads: 24
    shell:
        "plink2 --vcf {input.vcf_input} --keep-allele-order --double-id --make-bed --threads {threads} --out {params.outfix}"


rule generate_phenotypes:
    """Make file for each phenotype"""
    input:
        rscript="scripts/phenotypes/{phenotype}.R",
        ploidy_calls="/data/rmccoy22/natera_spectrum/karyohmm_outputs/compiled_output/natera_embryos_v2.karyohmm_v14.bph_sph_trisomy.071023.tsv.gz",
        metadata=metadata,
    output:
        phenotype_file="results/phenotypes/{phenotype}_by_{parent}.csv",
    wildcard_constraints:
        phenotype="maternal_meiotic_aneuploidy|haploidy|triploidy|embryo_count|parental_triploidy",
        parent="mother|father",
    params:
        bayes_factor_cutoff=2,
        nullisomy_min=5,
        ploidy_max=11,
        ploidy_min=12,
    run:
        command = "Rscript --vanilla {input.rscript} {output.phenotype_file} {wildcards.parent}"

        if wildcards.phenotype == "maternal_meiotic_aneuploidy":
            command += " {input.ploidy_calls} {params.bayes_factor_cutoff} {params.nullisomy_min} {params.ploidy_max}"
        elif wildcards.phenotype in ["haploidy", "triploidy", "parental_triploidy"]:
            command += " {input.ploidy_calls} {params.bayes_factor_cutoff} {params.ploidy_min}"
        elif wildcards.phenotype == "embryo_count": 
            command += " {input.metadata}"

        shell(command)


rule run_gwas:
    """Run GWAS for each set of parameters"""
    input:
        gwas_rscript="scripts/gwas/gwas_all.R",
        metadata=metadata,
        bed=rules.vcf2bed.output.bedfile,
        discovery_test=general_outputs_fp
        + "discover_validate_split_{parent}.txt",
        parental_pcs=rules.run_plink_pca.output.eigenvec,
        pheno=rules.generate_phenotypes.output.phenotype_file,
        bim=rules.vcf2bed.output.bimfile,
    output:
        gwas_output=gwas_results + "gwas_{phenotype}_by_{parent}_{dataset_type}_{chrom}.txt",
    threads: 32
    wildcard_constraints:
        dataset_type="discovery|test",
        phenotype="maternal_meiotic_aneuploidy|triploidy|haploidy|embryo_count|parental_triploidy",
        parent="mother|father",
    shell:
        "Rscript --vanilla {input.gwas_rscript} {input.metadata} {input.bed} {input.discovery_test} {input.parental_pcs} {input.pheno} {input.bim} {wildcards.dataset_type} {wildcards.phenotype} {wildcards.parent} {threads} {output.gwas_output}"

rule merge_chroms:
    """Create single file for each phenotype, merging all chromosomes"""
    input:
        expand(
            gwas_results + "gwas_{phenotype}_by_{parent}_{dataset_type}_{chrom}.txt",
            phenotype=phenotypes,
            parent=parents,
            dataset_type=dataset_type,
            chrom=range(1, 23),  
        ),
    output:
        merged_file = gwas_results + "gwas_{phenotype}_{parent}_total.txt.gz",
    shell:
        "cat {input} | gzip > {output.merged_file}"
        
