## Make phenotype file for embryos affected by maternal meiotic aneuploidy 
# (cn = 1p or 3m for a number of chromosomes between the thresholds)

# load libraries
library(data.table)
library(tidyr)
library(dplyr)

# Usage: 
# /scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/scripts/phenotypes/maternal_meiotic_aneuploidy.R \ 
# /scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/results/phenotypes/maternal_meiotic_aneuploidy_by_mother.csv \
# mother \
# /data/rmccoy22/natera_spectrum/karyohmm_outputs/compiled_output/natera_embryos.karyohmm_v18.bph_sph_trisomy.full_annotation.112023.filter_bad_trios.tsv.gz \
# 2 \
# 5 \ # 5 or more chromosomes at cn=0 is considered failed amplification 
# 3 \ # 3 or more aneuploid chromosomes is not considered "maternal aneuploidy" but rather another ploidy (number of chromosomes greater than which the embryo is not just "aneuploid" but rather has an entire ploidy) 
# /data/rmccoy22/natera_spectrum/data/summary_metadata/spectrum_metadata_merged.csv
# maternal_meiotic

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
# maximum number of chromosomes that are allowed aneuploidy
# anything more is a whole-genome gain/loss
ploidy_threshold <- as.numeric(args[6])
# metadata to filter for day5 embryos
metadata <- args[7]
# phenotype 
phenotype <- args[8]

# source Rscript with functions `filter_data`, `day5_only`, 
# `remove_null_wg`, and `count_ploidy_by_parent`
source("/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/scripts/phenotypes/helper_functions/get_ploidy.R")

# read in embryos and metadata 
embryos <- fread(embryos)
metadata <- fread(metadata)

# keep only day 5 embryos
embryos <- day5_only(embryos, metadata)
# keep only embryos that meet bf cutoff 
embryos <- filter_data(embryos, bayes_factor_cutoff)

# remove embryos with failed amplification 
embryos_filtered <- remove_failed_amp(embryos, !!as.name(parent),
                                       nullisomy_threshold)

# remove embryos with whole genome gain/loss 
embryos_filtered <- remove_wholegenome_gainloss(embryos_filtered,
                                                !!as.name(parent), 
                                                ploidy_threshold)

# count maternal meiotic aneuploidies per embryo, based on parent
# group ploidy by respective parent 
ploidy_counts_by_parent <- count_ploidy_by_parent(embryos_filtered, 
                                                  !!as.name(parent), 
                                                  phenotype, 
                                                  ploidy_threshold)
colnames(ploidy_counts_by_parent)[1] <- "array"

# write to file 
write.csv(ploidy_counts_by_parent, out_fname, row.names = FALSE)
