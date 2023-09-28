## load libraries
library(data.table)
library(tidyr)
library(dplyr)

# Usage: 
# /scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/scripts/phenotypes/maternal_haploidy.R \ 
# /scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/results/phenotypes/maternal_haploidy_by_mother.csv \
# mother \ 
# /data/rmccoy22/natera_spectrum/karyohmm_outputs/compiled_output/natera_embryos_v2.karyohmm_v14.bph_sph_trisomy.071023.tsv.gz \ 
# 2 \ 
# 10 \ 
# maternal_haploidy

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

# check input 
if (!(parent %in% c("mother", "father"))) {
  stop("Invalid 'parent' argument. Use 'mother' or 'father'.")
}

# source Rscript with functions `filter_data` and `count_ploidy_by_parent`
source("get_ploidy.R")

# read in and filter data 
embryos <- fread(embryos)
embryos_filtered <- filter_data(embryos, bayes_factor_cutoff)

# group ploidy by respective parent 
ploidy_counts_by_parent <- count_ploidy_by_parent(embryos, !!as.name(parent), phenotype, ploidy_threshold)
colnames(ploidy_counts_by_parent)[1] <- "array"

# write to file 
write.csv(ploidy_counts_by_parent, out_fname, row.names = FALSE)
