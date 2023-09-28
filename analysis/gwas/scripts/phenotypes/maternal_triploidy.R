## Make phenotype file for embryos affected by maternal triploidy (cn = 3m for a number above the threshold)

## load libraries
library(data.table)
library(tidyr)
library(dplyr)

# Usage: 
# /scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/scripts/phenotypes/maternal_triploidy.R \ 
# /scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/results/phenotypes/maternal_triploidy_by_mother.csv \
# mother \ 
# /data/rmccoy22/natera_spectrum/karyohmm_outputs/compiled_output/natera_embryos_v2.karyohmm_v14.bph_sph_trisomy.071023.tsv.gz \ 
# 2 \ 
# 10 \ 
# maternal_triploidy

# get command line arguments
args <- commandArgs(trailingOnly = TRUE)
# output file name
out_fname <- args[1]
# parent to measure phenotype
parent <- args[2]
# ploidy calls from karyohmm
embryos <- args[3]
# minimum bayes factor for filtering
bayes_factor_cutoff <- as.numeric(args[4])
# minimum number of chromosomes affected to declare phenotype 
ploidy_threshold <- as.numeric(args[5])
# name of phenotype 
phenotype <- args[6]

# source Rscript with functions `filter_data` and `count_ploidy_by_parent`
source("helper_functions/get_ploidy.R")

# read in and filter data 
embryos <- fread(embryos)

# generate phenotypes 
ploidy_counts_by_parent <- run_phenotype(embryos, bayes_factor_cutoff, parent, phenotype, ploidy_threshold)

# write to file 
write.csv(ploidy_counts_by_parent, out_fname, row.names = FALSE)
