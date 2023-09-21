#!python3

# Usage: nohup snakemake -p --cores 48 --snakefile gwas.smk > nohup_date.out 2>&1 &
# Optional: add -j 12 to submit as 12 jobs, etc.
# Optional: add -n to do a dry run
# Run from the filepath /scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/

# -------- Setting variables and paths for pre-GWAS processing steps ------- #
king_exec = "~/code/king"
king_outputs_fp = "results/"
vcf_fp = "/data/rmccoy22/natera_spectrum/genotypes/opticall_parents_031423/genotypes/"
alleleorder_fp = "results/opticall_concat_total.norm.b38.alleleorder"
discovery_validate_R = "scripts/discovery_test_split.R"
metadata = (
    "/data/rmccoy22/natera_spectrum/data/summary_metadata/spectrum_metadata_merged.csv"
)
fam_file = "/data/rmccoy22/natera_spectrum/genotypes/opticall_parents_031423/genotypes/opticall_concat_total.norm.b38.fam"
discovery_validate_out_fp = "results/"
pcs_out = "results/parental_genotypes_pcs/"


# -------- Setting key variables/paths for running GWAS across phenotypes in the Natera dataset ------- #
genotype_files = (
    "/data/rmccoy22/natera_spectrum/genotypes/opticall_parents_031423/genotypes/"
)
ploidy_calls = "/data/rmccoy22/natera_spectrum/karyohmm_outputs/compiled_output/natera_embryos.karyohmm_v11.052723.tsv.gz"
phenotype_scripts = "scripts/phenotypes/"
phenotype_results = "results/phenotypes/"
gwas_scripts = "scripts/gwas/"
gwas_results = "results/gwas/"


# Define the chromosomes that you will be running the pipeline on ...
chroms = range(1, 24)


# -------- Rule all to run whole pipeline -------- #
rule all:
    input:
        expand(
            gwas_results + "gwas_{phenotype}_by_{parent}_{dataset_type}_{chrom}.txt",
            phenotype="maternal_meiotic_aneuploidy",
            parent="mother",
            dataset_type="discovery",
            chrom=22,
        ),


# -------- Functions to determine run GWAS on maternal and paternal for each phenotype, and plot -------- #
rule run_king:
    """Reformat parental genotypes vcf and run king to identify related individuals"""
    input:
        vcf_input=vcf_fp + "opticall_concat_total.norm.b38.vcf.gz",
    output:
        king_bed=temp(
            king_outputs_fp + "opticall_concat_total.norm.b38.alleleorder.bed"
        ),
        king_bim=temp(
            king_outputs_fp + "opticall_concat_total.norm.b38.alleleorder.bim"
        ),
        king_fam=temp(
            king_outputs_fp + "opticall_concat_total.norm.b38.alleleorder.fam"
        ),
        king_log=temp(
            king_outputs_fp + "opticall_concat_total.norm.b38.alleleorder.log"
        ),
    params:
        plink_outfix=king_outputs_fp + "opticall_concat_total.norm.b38.alleleorder",
        king_outfix=king_outputs_fp,
    shell:
        """
        plink --vcf {input.vcf_input} --double-id --allow-extra-chr --make-bed --out {params.plink_outfix} 
        {king_exec} -b {output.king_bed} --unrelated --degree 2 --prefix {params.king_outfix}
        """


rule discovery_validate_split:
    """Split families into discovery/validation sets for use in GWAS"""
    input:
        discovery_validate_R=discovery_validate_R,
        metadata=metadata,
        fam_file=fam_file,
        king_to_remove=king_outputs_fp + "unrelated_toberemoved.txt",
    output:
        discovery_validate_maternal=discovery_validate_out_fp
        + "discover_validate_split_mother.txt",
        discovery_validate_paternal=discovery_validate_out_fp
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
        rscript=phenotype_scripts + "{phenotype}.R",
        ploidy_calls=ploidy_calls,
        metadata=metadata,
    output:
        phenotype_file=phenotype_results + "{phenotype}_by_{parent}.csv",
    wildcard_constraints:
        phenotype="maternal_meiotic_aneuploidy|haploidy|triploidy|embryo_count",
        parent="mother|father",
    params:
        nullisomy_min=5,
        ploidy_max=3,
        ploidy_min=20,
    run:
        command = "Rscript --vanilla {input.rscript} {output.phenotype_file} {wildcards.parent}"

        if wildcards.phenotype == "maternal_meiotic_aneuploidy":
            command += " {input.ploidy_calls} {params.nullisomy_min} {params.ploidy_max}"
        elif wildcards.phenotype in ["haploidy", "triploidy"]:
            command += " {input.ploidy_calls} {params.ploidy_min}"
        elif wildcards.phenotype == "embryo_count": 
            command += " {input.metadata}"

        shell(command)


rule run_gwas:
    """Run each GWAS"""
    input:
        gwas_rscript=gwas_scripts + "gwas_all.R",
        metadata=metadata,
        bed=rules.vcf2bed.output.bedfile,
        discovery_test=discovery_validate_out_fp
        + "discover_validate_split_{parent}.txt",
        parental_pcs=rules.run_plink_pca.output.eigenvec,
        pheno=rules.generate_phenotypes.output.phenotype_file,
        bim=rules.vcf2bed.output.bimfile,
        fam=rules.vcf2bed.output.famfile,
    output:
        gwas_output=gwas_results
        + "gwas_{phenotype}_by_{parent}_{dataset_type}_{chrom}.txt",
    wildcard_constraints:
        dataset_type="discovery|test",
        phenotype="maternal_meiotic_aneuploidy|embryo_count|triploidy|haploidy",
        parent="mother|father",
    shell:
        "Rscript --vanilla {input.gwas_rscript} {input.metadata} {input.bed} {input.discovery_test} {input.parental_pcs} {input.pheno} {input.bim} {wildcards.dataset_type} {wildcards.phenotype} {output.gwas_output}"
