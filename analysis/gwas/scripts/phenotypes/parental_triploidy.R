# load libraries
library(data.table)
library(tidyr)
library(dplyr)

# Usage: 
# /scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/scripts/phenotypes/parental_triploidy.R \ 
# /scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/results/phenotypes/parental_triploidy_by_mother.csv \
# mother \
# /data/rmccoy22/natera_spectrum/karyohmm_outputs/compiled_output/natera_embryos_v2.karyohmm_v14.bph_sph_trisomy.071023.tsv.gz \ 
# 2 \
# 20 \

# get command line arguments
args <- commandArgs(trailingOnly = TRUE)
out_fname <- args[1]
parent <- args[2]
input_data <- args[3]
bayes_factor_cutoff <- as.numeric(args[4])
triploidy_threshold <- as.numeric(args[5])

if (!(parent %in% c("mother", "father"))) {
  stop("Invalid 'parent' argument. Use 'mother' or 'father'.")
}

# read in data
input_data <- fread(input_data)
# keep only rows that have probabilities for all 6 cn states
embryos <- input_data[complete.cases(input_data[,c("0", "1m", "1p", "2", "3m", "3p")]),]
# filter bayes factors 
embryos <- embryos[embryos$bf_max > bayes_factor_cutoff,]


# find max posterior probability 
#selected_columns <- c("0", "1m", "1p", "2", "3m", "3p")
#highest_values <- apply(embryos[, selected_columns, with = FALSE], 1, function(x) max(x, na.rm = TRUE))
#embryos[, highest := highest_values]
# add column for most likely copy number 
embryos[, putative_cn := colnames(embryos[, 7:12])[apply(embryos[, 7:12], 1, which.max)]]
# create column for each chromosome number
embryos$chromosome <- gsub("chr", "", embryos$chrom) %>% as.integer()

# define function to count triploid embryos per parent column
count_triploid_by_parent <- function(data, triploidy_threshold, parent) {

  if(parent == "mother") {
  	data %>%
    group_by(mother, child) %>%
    summarise(num_trisomies = sum(putative_cn == "3m")) %>%
    mutate(is_triploid = if_else(num_trisomies >= triploidy_threshold, "aneu_true", "aneu_false")) %>%
    count(is_triploid) %>%
    pivot_wider(names_from = is_triploid, values_from = n, values_fill = 0) %>%
    replace(is.na(.), 0)
  } else if (parent == "father") {
  	data %>%
    group_by(father, child) %>%
    summarise(num_trisomies = sum(putative_cn == "3p")) %>%
    mutate(is_triploid = if_else(num_trisomies >= triploidy_threshold, "aneu_true", "aneu_false")) %>%
    count(is_triploid) %>%
    pivot_wider(names_from = is_triploid, values_from = n, values_fill = 0) %>%
    replace(is.na(.), 0)
  }
  # data %>%
  #   group_by({{ parent_column }}, child) %>%
  #   summarise(num_trisomies = if(parent == "mother") sum(putative_cn == "3m")
  #                           else if(parent == "father") sum(putative_cn == "3p")) %>%
  #   mutate(is_triploid = if_else(num_trisomies >= triploidy_threshold, "aneu_true", "aneu_false")) %>%
  #   count(is_triploid) %>%
  #   pivot_wider(names_from = is_triploid, values_from = n, values_fill = 0) %>%
  #   replace(is.na(.), 0)
}

# call function with the specified parent column
#triploid_counts_by_parent <- count_triploid_by_parent(embryos, !!as.name(parent), triploidy_threshold, parent)
triploid_counts_by_parent <- count_triploid_by_parent(embryos, triploidy_threshold, parent)
colnames(triploid_counts_by_parent)[1] <- "array"

# write to file 
write.csv(triploid_counts_by_parent, out_fname, row.names = FALSE)
