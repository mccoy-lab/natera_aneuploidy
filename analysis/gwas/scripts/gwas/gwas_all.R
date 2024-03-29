## GWAS for phenotypes from Natera 

# load libraries
library(data.table)
library(BEDMatrix)
library(dplyr)
library(pbmcapply)
library(purrr)

# Usage: /scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/scripts/gwas/gwas_all.R \
# "/data/rmccoy22/natera_spectrum/data/summary_metadata/spectrum_metadata_merged.csv" \
# "/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/results/gwas/opticall_concat_21.norm.b38.bed" \
# "/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/results/discover_validate_split_mother.txt" \
# "/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/results/parental_genotypes_pcs/parental_genotypes.eigenvec" \
# "/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/results/phenotypes/haploidy_by_mother.csv" \
# "/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/results/gwas/opticall_concat_21.norm.b38.bim" \
# "discovery" \ # dataset type 
# "haploidy" \ # phenotype 
# "mother" \ # parent 
# 16 \ # number of threads 
# "/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/results/gwas/gwas_haploidy_by_mother_discovery_9.tsv"

# get command line arguments 
args <- commandArgs(trailingOnly = TRUE)
# metadata for all years 
metadata <- args[1]
# parental genotypes as described https://github.com/mccoy-lab/natera_spectrum/tree/main/genotyping
bed <- args[2]
# discovery vs. test split for applicable parent
discovery_test <- args[3] 
# principal components for parents 
pcs <- args[4]
# phenotype file of interest
phenotype <- args[5]
# file containing snp name and map positions, corresponds with bed file 
bim <- args[6]
# string variable noting whether to GWAS using the discovery or test set for parents 
dataset_type <- args[7] 
# string variable naming the phenotype being considered (e.g., maternal_meiotic)
phenotype_name <- args[8] 
# string variable noting which parent to GWAS 
parent <- args[9]
# number of threads to use in mc.cores for GWAS 
threads <- as.numeric(args[10]) 
# output file name 
out_fname <- args[11]

# source Rscript with functions ``
source("/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/scripts/gwas/helper_functions/gwas_helper_functions.R")

# read in files 
metadata <- fread(metadata)
bed <- BEDMatrix(bed)
discovery_test <- fread(discovery_test)
pcs <- fread(pcs)
colnames(pcs)[1] <- "array"
phenotype <- fread(phenotype)
bim <- fread(bim) %>% 
  setnames(., c("chr", "snp_id", "drop", "pos", "ref", "alt"))


# Conduct GWAS across all sites 
gwas_results_dt <- run_gwas(dataset_type, discovery_test, metadata, bed, bim, pcs, phenotype, parent, phenotype_name, threads) 

# write to file
write.table(gwas_results_dt, out_fname, append = FALSE, sep = "\t", dec = ".",
            row.names = FALSE, col.names = TRUE, quote = FALSE)

