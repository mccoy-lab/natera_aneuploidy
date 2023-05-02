#!python3

import 


# Usage: nohup snakemake -p --cores 48 --snakefile gwas.smk > nohup_date.out 2>&1 &
# Optional: add -j 12 to submit as 12 jobs, etc.

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
  """ """
  shell: 
    """
    plink --vcf /scratch16/rmccoy22/abiddan1/natera_spectrum/genotyping/opticall_calls/opticall_concat_total.norm.b38.vcf.gz --double-id --allow-extra-chr --make-bed --out ./opticall_concat_total.norm.b38.alleleorder
    ~/code/king -b ./opticall_concat_total.norm.b38.alleleorder.bed --unrelated --degree 2
    """

rule compute_PCs: 
  """ """ 
  shell:
    "plink --vcf {input.genotypes.b38.vcf.gz} data/rmccoy22/natera_spectrum/genotypes/opticall_genotyping/opticall_concat_total.norm.b38.vcf.gz --double-id --allow-extra-chr --pca --out genotypes_pca"


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
