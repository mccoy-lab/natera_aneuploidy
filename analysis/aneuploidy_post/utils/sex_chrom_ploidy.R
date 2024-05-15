# Incidence of each sex chrom karyotype 
# Sex ratio of euploid embryos 

# load libraries
library(data.table)
library(dplyr)
library(ggplot2) 

# Usage: /aneuploidy_post/utils/sex_chrom_ploidy.R \ 
# "/data/rmccoy22/natera_spectrum/karyohmm_outputs/compiled_output/natera_embryos.karyohmm_v30a.sex_embryos.031624.tsv.gz" \
# "/data/rmccoy22/natera_spectrum/data/summary_metadata/spectrum_metadata_merged.csv"  

# get command line arguments 
args <- commandArgs(trailingOnly = TRUE)
# ploidy calls 
sex_chrom_calls <- args[1]
# metadata to filter to day 5s 
metadata <- args[2]

# read in data 
sex_chrom_calls <- fread(sex_chrom_calls)
metadata <- fread(metadata)

# get function `day5_only`
setwd(".")
source("../../gwas/scripts/phenotypes/helper_functions/phenotyping_helper_functions.R")

# filter data to keep only day 5 embryos 
sex_chrom_calls <- day5_only(sex_chrom_calls, metadata)

# add function to phenotying helper functions that filters the sex chrom data 


# plot karyotype status by year of case (2014-2020) 



# combine calls for x and y chromosomes 
sex_chrom_calls$karyotype <- paste0(sex_chrom_calls$x_maxBFcat, 
                                    sex_chrom_calls$y_maxBFcat)
# remove embryos for which there are no calls 
sex_chrom_calls <- sex_chrom_calls[sex_chrom_calls$karyotype != "nannan"]

# plot all karyotypes 
ggplot(sex_chrom_calls, aes(x = karyotype)) + 
  geom_bar()

# plot aneuploid karyotypes 
ggplot(sex_chrom_calls[!(sex_chrom_calls$karyotype %in% c("x1my1", "x2y0")),], 
       aes(x = karyotype)) + 
  geom_bar()
