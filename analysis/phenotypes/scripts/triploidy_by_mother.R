# load libraries
library(data.table)
library(tidyr)

# Usage: ./triploidy_by_mother.R \ 
# /data/rmccoy22/natera_spectrum/karyohmm_outputs/compiled_output/natera_embryos.karyohmm_v11.052723.tsv.gz \ 
# /scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/phenotypes/triploid_count_mother.csv

# get command line arguments
args <- commandArgs(trailingOnly = TRUE)
input_data <- args[1]
out_fname <- args[2]

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
triploid_counts_by_mother <- embryos %>% 
  group_by(mother, child) %>% 
  summarise(num_trisomies = sum(putative_cn == "3m" | putative_cn == "3p")) %>% 
  mutate(is_triploid = if_else(num_trisomies >= 20, "true", "false")) %>% 
  count(is_triploid) %>%
  pivot_wider(names_from = is_triploid, values_from = n, values_fill = 0) %>% 
  replace(is.na(.), 0)

# write to file 
write.csv(triploid_counts_by_mother, out_fname, row.names = FALSE)
