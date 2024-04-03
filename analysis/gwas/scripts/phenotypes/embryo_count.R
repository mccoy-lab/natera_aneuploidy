## Make embryo count phenotype for use in GWAS (weighted maternal age for use as covariate)

# load libraries
library(data.table)
library(tidyr)
library(dplyr)

# Usage: 
# /scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/scripts/phenotypes/embryo_count.R \ 
# /data/rmccoy22/natera_spectrum/data/summary_metadata/spectrum_metadata_merged.csv \
# /data/rmccoy22/natera_spectrum/karyohmm_outputs/compiled_output/natera_embryos.karyohmm_v30a.bph_sph_trisomy.full_annotation.031624.tsv.gz \
# mother \
# /scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/results/phenotypes/embryo_count_by_mother.csv \

# accept args from snakemake
args = commandArgs(trailingOnly = TRUE)
metadata <- args[1]
ploidy_calls <- args[2]
parent <- args[3]
out_fname <- args[4]

# read in metadata 
metadata <- fread(metadata)
# read in ploidy calls
ploidy_calls <- fread(ploidy_calls)

# source Rscript with functions `filter_data`, `day5_only`
source("/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/scripts/phenotypes/helper_functions/get_ploidy.R")

# get weighted maternal age across embryos, number of cycles, 
# number of embryos, and number of euploid (all chr cn=2) and 
# number of aneuploid (1 or more chr with other than cn=2)
count_embryos_by_parent <- function(ploidy_calls, metadata, parent) {
  
  # Count number of chromosomes for which cn is not 2
  cn <- c("0", "1m", "1p", "3m", "3p")
  
  # Count number of embryos for which each mother is aneuploid 
  aneuploidy_counts <- ploidy_calls %>%
    group_by(get(parent), child) %>%
    summarise(num_affected = sum(bf_max_cat %in% cn)) %>% 
    mutate(is_ploidy = ifelse(num_affected > 0, "aneu_true", "aneu_false")) %>%
    count(is_ploidy) %>%
    pivot_wider(names_from = is_ploidy, values_from = n, values_fill = 0) %>%
    rename("aneuploid" := aneu_true, "euploid" := aneu_false)
  
  # match column names with external data 
  colnames(aneuploidy_counts)[1] <- "array"
  
  # create new column that tags every individual affiliated with each parent even if in different caseIDs
  metadata_merged_array <- metadata %>%
    mutate(
      array_id_merged = ifelse(family_position == parent, paste0(array, "_merged"), NA_character_)
    ) %>%
    group_by(casefile_id) %>%
    fill(array_id_merged) %>%
    ungroup()
  
  # create dataframe that calculates the weighted age, number of embryos, and number of visits
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
  
  # Merge aneuploidy counts and weighted ages dataframes
  merged_table <- merge(weighted_ages, aneuploidy_counts, by = "array", all.x = TRUE)
  
  return(merged_table)
}

# Usage
combined_result <- combined_function(ploidy_calls, metadata, parent)

# write to file 
write.csv(combined_result, out_fname)