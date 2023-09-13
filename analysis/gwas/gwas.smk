#!python3

import 


# Usage: nohup snakemake -p --cores 48 --snakefile gwas.smk > nohup_date.out 2>&1 &
# Optional: add -j 12 to submit as 12 jobs, etc.

# -------- Setting variables and paths for pre-GWAS processing steps ------- #
king_exec <- "~/code/king"
king_outputs_fp <- "/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/results/"
vcf_input <- "/data/rmccoy22/natera_spectrum/genotypes/opticall_parents_031423/genotypes/opticall_concat_total.norm.b38.vcf.gz"
alleleorder_fp <- "/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/results/opticall_concat_total.norm.b38.alleleorder"
discovery_validate_R <- "/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/scripts/discovery_validate_split.R"
metadata = "/scratch16/rmccoy22/scarios1/natera_aneuploidy/data/spectrum_metadata_merged.csv"
fam_file <- "/data/rmccoy22/natera_spectrum/genotypes/opticall_parents_031423/genotypes/opticall_concat_total.norm.b38.fam"
discovery_validate_out_fp <- "/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/results/"
discovery_validate_maternal <- "/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/results/discover_validate_split_f.txt"
discovery_validate_paternal <- "/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/results/discover_validate_split_m.txt"
plot_discovery_validate <- "/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/results/discovery_validation_split_covariates.pdf"
pcs_out <- "/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/results/parental_genotypes_pcs"


# -------- Setting key variables/paths for running GWAS across phenotypes in the Natera dataset ------- #
linker_file = "" # two-column file with mom's arrayID with associate dad arrayID for use in paternal DNA gwas 
parent_pca = ".eigenvec"
genotype_files = "/data/rmccoy22/natera_spectrum/genotypes/opticall_parents_031423/genotypes/"
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
    king_outputs = king_outputs_fp,
  shell: 
    """
    plink --vcf {input.vcf_input} --double-id --allow-extra-chr --make-bed --out {input.alleleorder_fp}
    {king_exec} -b {output.alleleorder_bed} --unrelated --degree 2 {output.king_outputs}
    """

rule discovery_test_split: 
  """Split families into discovery/validation sets for use in GWAS"""
  input:
    discovery_validate_R = discovery_validate_R,
    metadata = metadata,
    fam_file = fam_file,
    king_to_remove = king_outputs_fp + "kingunrelated_toberemoved.txt",
  output: 
    discovery_validate_maternal = discovery_validate_out_fp + "discover_validate_split_mat.txt",
    discovery_validate_paternal = discovery_validate_out_fp + "discover_validate_split_pat.txt",
    plot_discovery_validate = discovery_validate_out_fp + "discover_validate_plots.pdf",
  shell: 
    "Rscript {input.discovery_validate_R} {input.metadata} {input.fam_file} {input.king_to_remove} {output.discovery_validate_maternal} {output.discovery_validate_paternal} {output.plot_discovery_validate}" 

rule compute_PCs: 
  """Run plink to get principal components for parental genotypes""" 
  input: 
    vcf_input = vcf_input,
  output: 
    pcs_out = pcs_output,
  shell:
    "plink --vcf {input.vcf_input} --double-id --allow-extra-chr --pca --out {output.pcs}" 
    
rule gwas_maternal_m_meiotic: 
  """Run GWAS on association with maternal meiotic error and maternal genotypes""" 
  input:
    gwas_Rscript_maternal,
    parental_pcs = pcs_out + parent_pca,
    bed = genotype_files + "opticall_concat_total.norm.b38.bed",
    bim = genotype_files + "opticall_concat_total.norm.b38.bim",
    pheno = phenotype_input_path + "ploidy_by_fam_discovery.txt",
    metadata,
    discovery_validate_maternal,
    out_fname = gwas_output_path + maternal_meiotic_mat_discovery_{chrom}.txt
  output: 
    out_fname,
  params: 
    ## 
  shell: 
    "Rscript {input.gwas_Rscript_maternal} {input.parental_pcs} {input.bed} {input.bim} {input.pheno} {input.metadata} {input.discovery_validate_maternal} {input.out_fname} {output.out_fname}"
    
rule gwas_paternal_m_meiotic: 
  """Run GWAS on association with maternal meiotic error and paternal genotypes""" 
  input:
    gwas_Rscript_paternal,
    parental_pcs = pcs_out,
  output: 
    gwas_paternal_all_sites = gwas_output_path + "paternal_m_meiotic.txt",
    gwas_paternal_MAF = gwas_output_path + "paternal_m_meiotic_MAF.txt",
  params: 
    ## 
  shell: 
    "Rscript {input.gwas_Rscript_paternal} {input.parental_pcs} {output.gwas_paternal_all_sites} {output.gwas_paternal_MAF}"

rule plot_manhattan_qq: 
  """Create qq plot for maternal meiotic""" 
  input:
    ## 
  output: 
    ## 
  params: 
    ## 
  shell: 
    ## 
