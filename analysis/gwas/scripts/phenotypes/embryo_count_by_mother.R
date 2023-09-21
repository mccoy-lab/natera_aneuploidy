## Make embryo count phenotype for use in GWAS (weighted maternal age for use as covariate)

# load libraries
library(data.table)
library(dplyr)
library(tidyr)

# Usage: ./embryo_count_by_mother.R \ 
# /scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/phenotypes/embryo_count_by_mother.csv \
# /data/rmccoy22/natera_spectrum/data/summary_metadata/spectrum_metadata_merged.csv \

# accept args from snakemake
args = commandArgs(trailingOnly = TRUE)
out_fname <- args[1]
metadata_fp <- args[2]

# read in metadata 
metadata <- fread(metadata_fp)

# create new column that tags every individual affiliated with each mother even if in different caseIDs
metadata_merged_array <- metadata %>%
  mutate(array_id_merged = ifelse(family_position == "mother", paste0(array, "_merged"), NA_character_)) %>%
  group_by(casefile_id) %>%
  fill(array_id_merged) %>%
  ungroup() 

# create dataframe that calculates the weighted age (based on age at each embryo), number of embryos, and number of visits 
weighted_ages <- metadata_merged_array %>%
  filter(family_position == "child") %>%
  group_by(array_id_merged) %>%
  summarise(weighted_age = sum(patient_age) / n(),
            child_count = n(), 
            num_visits = length(unique(patient_age))) %>% 
  as.data.frame()

# remove "_merged" from array column to allow easier downstream intersection 
weighted_ages$array <- gsub('_merged', '', weighted_ages$array_id_merged)
  
# write to file 
write.csv(weighted_ages, out_fname)