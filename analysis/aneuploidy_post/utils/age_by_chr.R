# Age models for aneuploidy of each chromosome

# load libraries
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2) 

# Usage: /aneuploidy_post/utils/age_by_chr.R \ 
# "/data/rmccoy22/natera_spectrum/karyohmm_outputs/compiled_output/natera_embryos.karyohmm_v30a.bph_sph_trisomy.full_annotation.031624.tsv.gz" \
# "/data/rmccoy22/natera_spectrum/data/summary_metadata/spectrum_metadata_merged.csv"
# 100 


# get command line arguments 
args <- commandArgs(trailingOnly = TRUE)
ploidy_calls <- args[1]
metadata <- args[2]
# minimum number of patients in an age group for that age to be included 
min_age_count <- args[3]


# get functions `filter_data` and `day5_only`
setwd(".")
source("../../gwas/scripts/phenotypes/helper_functions/phenotyping_helper_functions.R")

# filter data 
ploidy_calls_filtered <- filter_data(ploidy_calls, parent, 
                                     bayes_factor_cutoff = 2, 
                                     nullisomy_threshold = 5, min_prob = 0.9) 
ploidy_calls_filtered <- day5_only(ploidy_calls_filtered, metadata)


# age effect of specific chr (15, 16, 21, 22) 
# table where each mother is in it 22 times, once for each chromosome 
# columns are maternal array, weighted age across all embryos, number with trisomy for that chr, number with monosomy for that chr, total 

aneu_by_chr <- ploidy_calls_filtered %>% 
  group_by(mother, chrom) %>% 
  summarise(num_monosomies = sum(bf_max_cat == "1p"), 
            num_trisomies = sum(bf_max_cat == "1m"))

aneu_by_chr_age <- merge(aneu_by_chr, )









