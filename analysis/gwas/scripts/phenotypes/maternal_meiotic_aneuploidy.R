# load libraries
library(data.table)
library(tidyr)
library(dplyr)

# Usage: ./maternal_meiotic_aneuploidy.R \ 
# /scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/phenotypes/maternal_meiotic_aneuploidy_mother.csv \
# mother \
# /data/rmccoy22/natera_spectrum/karyohmm_outputs/compiled_output/natera_embryos.karyohmm_v11.052723.tsv.gz \ 
# 5 \ # 5 or more chromosomes at cn=0 is considered failed amplification 
# 3 # 3 or more aneuploid chromosomes is not considered "maternal aneuploidy" but rather another ploidy 

# get command line arguments
args <- commandArgs(trailingOnly = TRUE)
out_fname <- args[1]
parent <- args[2]
input_data <- args[3]
nullisomy_threshold <- args[4]
ploidy_threshold <- args[5] 
# number of chromosomes greater than which the embryo is not just "aneuploid" but rather has an entire ploidy


# read in data
input_data <- fread(input_data)
# keep only rows that have probabilities for all 6 cn states
embryos <- input_data[complete.cases(input_data[,c("0", "1m", "1p", "2", "3m", "3p")]),]
# filter bayes factors 
embryos <- embryos[embryos$bf_max > 2,]

# find max posterior probability 
selected_columns <- c("0", "1m", "1p", "2", "3m", "3p")
highest_values <- apply(embryos[, selected_columns, with = FALSE], 1, function(x) max(x, na.rm = TRUE))
embryos[, highest := highest_values]
# add column that says what the copy number is 
embryos[, putative_cn := colnames(embryos[, 10:15])[apply(embryos[, 10:15], 1, which.max)]]
# create new column that is just the chromosome number 
embryos$chromosome <- gsub("chr", "", embryos$chrom) %>% as.integer()


# remove embryos with failed amplification, triploidies, or haploidies 
# grab embryos with 5 or more nullisomies (suggests failed amplification)
count_nullisomies <- embryos %>% 
  group_by({{parent}}, child) %>% 
  summarise(num_nullisomies = sum(putative_cn == "0"))
successful_amp <- count_nullisomies[count_nullisomies$num_nullisomies < nullisomy_threshold,]
# grab triploid embryos 
count_triploidies <- embryos %>% 
  group_by({{parent}}, child) %>% 
  summarise(num_trisomies = sum(putative_cn == "3m" | putative_cn == "3p")) 
non_trip <- count_triploidies[count_triploidies$num_trisomies < ploidy_threshold,]              
# grab haploid embryos (this would also catch isoUPD embryos)
count_haploidies <- embryos %>% 
  group_by({{parent}}, child) %>% 
  summarise(num_monosomies = sum(putative_cn == "1m" | putative_cn == "1p"))
non_hap <- count_haploidies[count_haploidies$num_monosomies < ploidy_threshold,]
# remove failed amp, triploid, haploid embryos 
embryos_filtered <- embryos[embryos$child %in% successful_amp$child & embryos$child %in% non_trip$child & embryos$child %in% non_hap$child]


# count maternal meiotic aneuploidies per embryo
embryo_counts_by_parent <- embryos_filtered %>% 
  group_by({{parent}}, child) %>% 
  summarise(mat_aneu = sum(putative_cn == "3m" | putative_cn == "1p")) %>% 
  mutate(is_aneu = if_else(mat_aneu > 0, "true", "false")) %>% 
  count(is_aneu) %>%
  pivot_wider(names_from = is_aneu,  values_from = n, values_fill = 0, names_prefix = "aneu_") %>% 
  replace(is.na(.), 0)


# write to file 
write.csv(embryo_counts_by_parent, out_fname, row.names = FALSE)
