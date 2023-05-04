#!python3

import 


# Usage: nohup snakemake -p --cores 48 --snakefile gwas.smk > nohup_date.out 2>&1 &
# Optional: add -j 12 to submit as 12 jobs, etc.

# -------- Setting variables and paths for pre-GWAS processing steps ------- #
king_exec <- "~/code/king"
vcf_input <- "/data/rmccoy22/natera_spectrum/genotypes/opticall_parents_031423/genotypes/opticall_concat_total.norm.b38.vcf.gz"
alleleorder_fp <- "/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/results/opticall_concat_total.norm.b38.alleleorder"
pcs_output <- "/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/results/parental_genotypes_pcs"
plot_discovery_validate <- "/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/results/discovery_validation_split_covariates.pdf"
discovery_validate_maternal <- "/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/results/discover_validate_split_f.txt"
discovery_validate_paternal <- "/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/results/discover_validate_split_m.txt"


# -------- Setting key variables/paths for running GWAS across phenotypes in the Natera dataset ------- #
metadata = "spectrum_metadata_merged.csv"
linker_file = "" # two-column file with mom's arrayID with associate dad arrayID for use in paternal DNA gwas 
parent_pca = ".eigenvec"
gwas_Rscript_maternal = "/maternal_gwas.R"
gwas_Rscript_maternal = "/paternal_gwas.R"
phenotype_input_path = "/natera_spectrum/gwas/phenotyping/phenotype_files/" #update this for Rockfish and match to output from other file
gwas_output_path = "/natera_spectrum/gwas/output_files"



# -------- Functions to determine run GWAS on maternal and paternal for each phenotype, and plot -------- #
rule run_king: 
  """Reformat parental genotypes vcf and run king to identify related individuals"""
  input: 
    vcf_input = vcf_input,
    alleleorder_fp = alleleorder_fp,
  output: 
    alleleorder_bed = alleleorder_fp + ".bed",
  shell: 
    """
    plink --vcf {input.vcf_input} --double-id --allow-extra-chr --make-bed --out {input.alleleorder_fp}
    {king_exec} -b {output.alleleorder_bed} --unrelated --degree 2
    """

rule compute_PCs: 
  """Run plink to get principal components for parental genotypes""" 
  input: 
    vcf_input = vcf_input,
  output: 
    pcs = pcs_output
  shell:
    "plink --vcf {input.vcf_input} --double-id --allow-extra-chr --pca --out {output.pcs}"


rule discovery_test_split: 
  """ """
  shell: 
    "Rscript " 

rule run_maternal: 
  """ """ 
  input:
    gwas_Rscript_maternal,
  output: 
    gwas_maternal_all_sites=gwas_output_path + 
    gwas_maternal_MAF=gwas_output_path + 
  params: 
    ## 
  shell: 
    "Rscript {input.gwas_Rscript_maternal} {output.gwas_maternal_all_sites} {output.gwas_maternal_MAF}
    
rule run_paternal: 
  """ """ 
  input:
    ## 
  output: 
    ## 
  params: 
    ## 
  shell: 
    ## 

rule plot_manhattan_qq: 
  """ """ 
  input:
    ## 
  output: 
    ## 
  params: 
    ## 
  shell: 
    ## 
