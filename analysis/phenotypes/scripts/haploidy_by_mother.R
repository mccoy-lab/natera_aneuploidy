# load libraries
library(data.table)
library(tidyr)

# Usage: ./haploidy_by_mother.R \ 
# /data/rmccoy22/natera_spectrum/karyohmm_outputs/compiled_output/natera_embryos.karyohmm_v11.052723.tsv.gz \ 
# /scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/phenotypes/haploid_count_mother.csv \
# 20

# get command line arguments
args <- commandArgs(trailingOnly = TRUE)
input_data <- args[1]
out_fname <- args[2]
haploidy_threshold <- args[3]

# read in data
input_data <- fread(input_data)
# keep only rows that have probabilities for all 6 cn states
embryos <- input_data[complete.cases(input_data[,c("0", "1m", "1p", "2", "3m", "3p")]),]


# find max posterior probability 
selected_columns <- c("0", "1m", "1p", "2", "3m", "3p")
highest_values <- apply(embryos[, selected_columns, with = FALSE], 1, function(x) max(x, na.rm = TRUE))
embryos[, highest := highest_values]
# add column for most likely copy number 
embryos[, putative_cn := colnames(embryos[, 10:15])[apply(embryos[, 10:15], 1, which.max)]]
# create column for each chromosome number
embryos$chromosome <- gsub("chr", "", embryos$chrom) %>% as.integer()

# count triploid embryos per mother
haploid_counts_by_mother <- embryos %>% 
  group_by(mother, child) %>% 
  summarise(num_monosomies = sum(putative_cn == "1m" | putative_cn == "1p")) %>% 
  mutate(is_haploid = if_else(num_monosomies >= haploidy_threshold, "true", "false")) %>% 
  count(is_haploid) %>%
  pivot_wider(names_from = is_haploid, values_from = n, values_fill = 0) %>% 
  replace(is.na(.), 0)

# write to file 
write.csv(haploid_counts_by_mother, out_fname, row.names = FALSE)
