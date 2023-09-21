## Make embryo count phenotype for use in GWAS (weighted maternal age for use as covariate)

# load libraries
library(data.table)
library(dplyr)
library(tidyr)

# Usage: ./embryo_count_by_mother.R \ 
# /scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/phenotypes/embryo_count_by_mother.csv \
# mother \
# /data/rmccoy22/natera_spectrum/data/summary_metadata/spectrum_metadata_merged.csv \

# accept args from snakemake
args = commandArgs(trailingOnly = TRUE)
out_fname <- args[1]
parent <- args[2]
metadata <- args[3]

# read in metadata 
metadata <- fread(metadata)

create_merged_array_and_weighted_ages <- function(metadata, parent) {

  # Create new column that tags every individual affiliated with each parent even if in different caseIDs
  metadata_merged_array <- metadata %>%
    mutate(
      array_id_merged = ifelse(family_position == parent, paste0(array, "_merged"), NA_character_)
    ) %>%
    group_by(casefile_id) %>%
    fill(array_id_merged) %>%
    ungroup()

  # Create dataframe that calculates the weighted age, number of embryos, and number of visits
  weighted_ages <- metadata_merged_array %>%
    filter(family_position == "child") %>%
    group_by(array_id_merged) %>%
    summarise(
      weighted_age = sum(patient_age) / n(),
      num_embryos = n(),
      num_visits = length(unique(patient_age))
    ) %>%
    as.data.frame()

  # Remove "_merged" from array column to allow easier downstream intersection
  weighted_ages$array <- gsub('_merged', '', weighted_ages$array_id_merged)

  return(weighted_ages)
}

result <- create_merged_array_and_weighted_ages(metadata, parent)

  
# write to file 
write.csv(result, out_fname)