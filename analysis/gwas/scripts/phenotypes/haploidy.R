# load libraries
library(data.table)
library(tidyr)
library(dplyr)

# Usage: 
# /scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/scripts/phenotypes/haploidy.R \ 
# /scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/results/phenotypes/haploidy_by_mother.csv \
# mother \
# /data/rmccoy22/natera_spectrum/karyohmm_outputs/compiled_output/natera_embryos_v2.karyohmm_v14.bph_sph_trisomy.071023.tsv.gz \ 
# 2 \
# 20

# get command line arguments
args <- commandArgs(trailingOnly = TRUE)
out_fname <- args[1]
parent <- args[2]
input_data <- args[3]
bayes_factor_cutoff <- as.numeric(args[4])
haploidy_threshold <- as.numeric(args[5])

# read in data
input_data <- fread(input_data)
# keep only rows that have probabilities for all 6 cn states
embryos <- input_data[complete.cases(input_data[,c("0", "1m", "1p", "2", "3m", "3p")]),]
# filter bayes factors 
embryos <- embryos[embryos$bf_max > bayes_factor_cutoff,]

# find max posterior probability 
selected_columns <- c("0", "1m", "1p", "2", "3m", "3p")
highest_values <- apply(embryos[, selected_columns, with = FALSE], 1, function(x) max(x, na.rm = TRUE))
embryos[, highest := highest_values]
# add column for most likely copy number 
embryos[, putative_cn := colnames(embryos[, 7:12])[apply(embryos[, 7:12], 1, which.max)]]
# create column for each chromosome number
embryos$chromosome <- gsub("chr", "", embryos$chrom) %>% as.integer()

# define function to count haploid embryos per parent column
count_haploid_by_parent <- function(data, parent_column, haploidy_threshold) {
  data %>%
    group_by({{ parent_column }}, child) %>%
    summarise(num_monosomies = sum(putative_cn == "1m" | putative_cn == "1p")) %>%
    mutate(is_haploid = if_else(num_monosomies >= haploidy_threshold, "aneu_true", "aneu_false")) %>%
    count(is_haploid) %>%
    pivot_wider(names_from = is_haploid, values_from = n, values_fill = 0) %>%
    replace(is.na(.), 0)
}

# call function with the specified parent column
haploid_counts_by_parent <- count_haploid_by_parent(embryos, !!as.name(parent), haploidy_threshold)
colnames(haploid_counts_by_parent)[1] <- "array"

# write to file 
write.csv(haploid_counts_by_parent, out_fname, row.names = FALSE)
