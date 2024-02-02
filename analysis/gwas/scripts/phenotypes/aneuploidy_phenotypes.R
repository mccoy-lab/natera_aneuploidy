## Make phenotype file for embryos affected by aneuploidy phenotypes 

# load libraries
library(data.table)
library(tidyr)
library(dplyr)

# Usage: 
# /scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/scripts/phenotypes/aneuploidy_phenotypes.R \ 
# /scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/results/phenotypes/maternal_meiotic_aneuploidy_by_mother.csv \
# mother \ #parent to sort results by 
# /data/rmccoy22/natera_spectrum/karyohmm_outputs/compiled_output/natera_embryos.karyohmm_v18.bph_sph_trisomy.full_annotation.112023.filter_bad_trios.tsv.gz \
# 2 \ #remove chromosomes with bayes factor > 2
# 5 \ #5 or more chromosomes at cn=0 is considered failed amplification 
# /data/rmccoy22/natera_spectrum/data/summary_metadata/spectrum_metadata_merged.csv \
# maternal_meiotic \ #phenotype name 
# TRUE \ 
# 0.9 \ #minimum posterior probability for cn state 
# 3 \ #3 or more aneuploid chromosomes is not considered "maternal aneuploidy" but rather another ploidy \

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
# maximum number of chromosomes that are allowed nullisomies
# anything more is considered failed amplification
nullisomy_threshold <- as.numeric(args[5])
# metadata to filter for day5 embryos
metadata <- args[6]
# phenotype 
phenotype <- args[7]
# whether to keep only day 5 embryos 
filter_day_5 <- args[8]
# min posterior probability of cn state 
min_prob <- args[9]
# maximum number of chromosomes that are allowed aneuploidy
# anything more is a whole-genome gain/loss
ploidy_threshold <- as.numeric(args[10])

# source Rscript with functions `filter_data`, `day5_only`, 
# `count_ploidy_by_parent`, and `run_phenotype`
source("/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/scripts/phenotypes/helper_functions/get_ploidy.R")

# read in embryos and metadata 
embryos <- fread(embryos)
metadata <- fread(metadata)

# keep only high-quality embryos (remove noisy, high-bayes factor, and
# potential failed amplification embryos) and day 5 embryos
# count number of aneuploid and non-aneuploid embryos per parent 
ploidy_counts_by_parent <- run_phenotype(embryos, !!as.name(parent), metadata, 
                                         phenotype, filter_day_5 = TRUE, 
                                         bayes_factor_cutoff = 2, 
                                         nullisomy_threshold = 5, 
                                         ploidy_threshold = 3, min_prob = 0.9)

# write to file 
write.csv(ploidy_counts_by_parent, out_fname, row.names = FALSE)
