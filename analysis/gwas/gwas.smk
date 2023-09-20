#!python3

# Usage: nohup snakemake -p --cores 48 --snakefile gwas.smk > nohup_date.out 2>&1 &
# Optional: add -j 12 to submit as 12 jobs, etc.

# -------- Setting variables and paths for pre-GWAS processing steps ------- #
king_exec = "~/code/king"
king_outputs_fp = (
    "/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/results/"
)
vcf_fp = "/data/rmccoy22/natera_spectrum/genotypes/opticall_parents_031423/genotypes/"
alleleorder_fp = "/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/results/opticall_concat_total.norm.b38.alleleorder"
discovery_validate_R = "/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/scripts/discovery_validate_split.R"
metadata = (
    "/data/rmccoy22/natera_spectrum/data/summary_metadata/spectrum_metadata_merged.csv"
)
fam_file = "/data/rmccoy22/natera_spectrum/genotypes/opticall_parents_031423/genotypes/opticall_concat_total.norm.b38.fam"
discovery_validate_out_fp = (
    "/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/results/"
)
plot_discovery_validate = "/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/results/discovery_validation_split_covariates.pdf"
pcs_out = "/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/results/parental_genotypes_pcs/"


# -------- Setting key variables/paths for running GWAS across phenotypes in the Natera dataset ------- #
linker_file = ""  # two-column file with mom's arrayID with associate dad arrayID for use in paternal DNA gwas
parent_pca = "genotypes_pca.eigenvec"
genotype_files = (
    "/data/rmccoy22/natera_spectrum/genotypes/opticall_parents_031423/genotypes/"
)
gwas_Rscript_maternal = "/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/scripts/maternal_gwas.R"
gwas_Rscript_paternal = "/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/scripts/paternal_gwas.R"
phenotype_input_path = "/natera_spectrum/gwas/phenotyping/phenotype_files/"  # update this for Rockfish and match to output from other file
gwas_output_path = "/natera_spectrum/gwas/output_files/"


# Define the chromosomes that you will be running the pipeline on ...
chroms = range(1, 24)


# -------- Rule all to run whole pipeline -------- #
rule all:
    input:
        expand(
            gwas_output_path + "maternal_meiotic_mat_{dataset_type}_{chrom}.txt",
            dataset_type="discovery",
            chrom=chroms,
        ),


# -------- Functions to determine run GWAS on maternal and paternal for each phenotype, and plot -------- #
rule run_king:
    """Reformat parental genotypes vcf and run king to identify related individuals"""
    input:
        vcf_input=vcf_fp + "opticall_concat_total.norm.b38.vcf.gz",
        alleleorder_fp=alleleorder_fp,
    output:
        alleleorder_bed=alleleorder_fp + ".bed",
        king_outputs=king_outputs_fp,
    shell:
        """
        plink --vcf {input.vcf_input} --double-id --allow-extra-chr --make-bed --out {input.alleleorder_fp} 
        {king_exec} -b {output.alleleorder_bed} --unrelated --degree 2 {output.king_outputs}
        """


rule discovery_validate_split:
    """Split families into discovery/validation sets for use in GWAS"""
    input:
        discovery_validate_R=discovery_validate_R,
        metadata=metadata,
        fam_file=fam_file,
        king_to_remove=king_outputs_fp + "kingunrelated_toberemoved.txt",
    output:
        discovery_validate_maternal=discovery_validate_out_fp
        + "discover_validate_split_mat.txt",
        discovery_validate_paternal=discovery_validate_out_fp
        + "discover_validate_split_pat.txt",
        plot_discovery_validate=discovery_validate_out_fp
        + "discover_validate_plots.pdf",
    shell:
        "Rscript --vanilla {input.discovery_validate_R} {input.metadata} {input.fam_file} {input.king_to_remove} {output.discovery_validate_maternal} {output.discovery_validate_paternal} {output.plot_discovery_validate}"


rule compute_PCs:
    """Run plink to get principal components for parental genotypes"""
    input:
        vcf_input=vcf_fp + "opticall_concat_total.norm.b38.vcf.gz",
    output:
        pcs_out=pcs_output,
    shell:
        "plink --vcf {input.vcf_input} --double-id --allow-extra-chr --pca --out {output.pcs_out}"


rule vcf2bed:
    """Take vcf for each chromosome, convert to bed as a temp file for use in the gwas"""
    input:
        vcf_input=vcf_fp + "opticall_concat_{chrom}.norm.b38.vcf.gz",
    output:
        bedfile=temp(vcf_input + "opticall_concat_{chrom}.norm.b38.bed"),
        bimfile=temp(vcf_input + "opticall_concat_{chrom}.norm.b38.bim"),
        famfile=temp(vcf_input + "opticall_concat_{chrom}.norm.b38.fam"),
        logfile=temp(vcf_input + "opticall_concat_{chrom}.norm.b38.log"),
    params:
        out_prefix=lambda wildcards: vcf_fp + f"opticall_concat_{wildcards.chrom}.norm.b38",
    shell:
        "plink --vcf {input.vcf_input} --keep-allele-order --double-id --make-bed --out {params.out_prefix}"


rule gwas_maternal_m_meiotic:
    """Maternal GWAS with maternal meiotic error"""
    input:
        metadata,
        gwas_Rscript_maternal=gwas_Rscript_maternal,
        parental_pcs=rules.compute_PCs.output.pcs_out,
        bed=rules.vcf2bed.output.bedfile,
        bim=rules.vcf2bed.output.bimfile,
        pheno=phenotype_input_path + "ploidy_by_fam_{dataset_type}.txt",
        discovery_test=rules.discovery_validate_split.output.discovery_validate_maternal,
    output:
        gwas_maternal_m_meiotic_out=gwas_output_path
        + "maternal_meiotic_mat_{dataset_type}_{chrom}.txt",
    wildcard_constraints:
        dataset_type="discovery|test"
    shell:
        "Rscript --vanilla {input.gwas_Rscript_maternal} {input.parental_pcs} {input.bed} {input.bim} {input.pheno} {input.metadata} {input.discovery_test} {wildcards.dataset_type} {output.gwas_maternal_m_meiotic_out}"
