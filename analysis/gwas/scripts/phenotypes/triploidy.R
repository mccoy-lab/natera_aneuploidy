# load libraries
library(data.table)
library(tidyr)
library(dplyr)

# Usage: ./triploidy_by_mother.R \ 
# /scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/phenotypes/triploidy_by_mother.csv \
# mother \
# /data/rmccoy22/natera_spectrum/karyohmm_outputs/compiled_output/natera_embryos.karyohmm_v11.052723.tsv.gz \ 
# 20

# get command line arguments
args <- commandArgs(trailingOnly = TRUE)
out_fname <- args[1]
parent <- args[2]
input_data <- args[3]
triploidy_threshold <- args[4]

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

# Define a function to count triploid embryos per parent column
count_triploid_by_parent <- function(data, parent_column, triploidy_threshold) {
  data %>%
    group_by({{ parent_column }}, child) %>%
    summarise(num_trisomies = sum(putative_cn == "3m" | putative_cn == "3p")) %>%
    mutate(is_triploid = if_else(num_trisomies >= triploidy_threshold, "aneu_true", "aneu_false")) %>%
    count(is_triploid) %>%
    pivot_wider(names_from = is_triploid, values_from = n, values_fill = 0) %>%
    replace(is.na(.), 0)
}

# Call the function with the specified parent column
triploid_counts_by_parent <- count_triploid_by_parent(embryos, !!as.name(parent), triploidy_threshold)
colnames(triploid_counts_by_parent)[1] <- "array"

# write to file 
write.csv(triploid_counts_by_parent, out_fname, row.names = FALSE)
