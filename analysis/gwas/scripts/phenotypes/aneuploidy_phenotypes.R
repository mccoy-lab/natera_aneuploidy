## Make phenotype file for embryos affected by aneuploidy phenotypes 

# load libraries
library(data.table)
library(tidyr)
library(dplyr)

# Usage: 
# /scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/scripts/phenotypes/aneuploidy_phenotypes.R \ 
# /data/rmccoy22/natera_spectrum/karyohmm_outputs/compiled_output/natera_embryos.karyohmm_v18.bph_sph_trisomy.full_annotation.112023.filter_bad_trios.tsv.gz \
# mother \ # parent to group phenotype by  
# /data/rmccoy22/natera_spectrum/data/summary_metadata/spectrum_metadata_merged.csv \
# maternal_meiotic \ # phenotype name 
# TRUE 
# 2 \ remove chr with bayes factor > bayes_factor_cutoff
# 5 \ remove embryos that had more chr with cn = 0 for than nullisomy_threshold 
# 0.9 \ minimum posterior probability for each cn call 
# 3 \ max number of affected chr to count for maternal meiotic phenotype
# 15 \ min number of affected chr to count for whole genome gain/loss
# /scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/results/phenotypes/maternal_meiotic_aneuploidy_by_mother.csv \

# get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# assign arguments to variables 
# ploidy calls for each embryo's chromosomes, from karyoHMM
ploidy_calls <- args[1]
# parent to measure phenotype
parent <- args[2]
# metadata to filter for day5 embryos
metadata <- args[3]
# phenotype 
phenotype <- args[4]
# whether to keep only day 5 embryos 
filter_day_5 <- args[5]
# minimum bayes factor for filtering
bayes_factor_cutoff <- as.numeric(args[6])
# how many cn = 0 before it's failed amplification 
nullisomy_threshold <- as.numeric(args[7])
# min posterior probability of cn state 
min_prob <- as.numeric(args[8])
# max number of chr for maternal meiotic pheno
max_meiotic <- as.numeric(args[9])
# min number of chr for whole genome gain/loss 
min_ploidy <- as.numeric(args[10])
# output file name
out_fname <- args[11]

# source Rscript with functions `filter_data`, `day5_only`, 
# `count_ploidy_by_parent`, and `run_phenotype`
source("/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/scripts/phenotypes/helper_functions/get_ploidy.R")

# read in embryos and metadata 
ploidy_calls <- fread(ploidy_calls)
metadata <- fread(metadata)

# keep only high-quality embryos (remove noisy, high-bayes factor, and
# potential failed amplification embryos) and day 5 embryos
# count number of aneuploid and non-aneuploid embryos per parent 
ploidy_counts_by_parent <- run_phenotype(ploidy_calls, parent, metadata, 
                                         phenotype, filter_day_5, 
                                         bayes_factor_cutoff, 
                                         nullisomy_threshold, min_prob, 
                                         max_meiotic, min_ploidy)

# write to file 
write.csv(ploidy_counts_by_parent, out_fname, row.names = FALSE)
