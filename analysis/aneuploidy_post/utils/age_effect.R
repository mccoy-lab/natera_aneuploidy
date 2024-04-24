# Plot errors affecting maternal and paternal errors by maternal or paternal age

# load libraries
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2) 

# Usage: /aneuploidy_post/utils/age_effect.R \ 
# "/data/rmccoy22/natera_spectrum/karyohmm_outputs/compiled_output/natera_embryos.karyohmm_v30a.bph_sph_trisomy.full_annotation.031624.tsv.gz" \
# "/data/rmccoy22/natera_spectrum/data/summary_metadata/spectrum_metadata_merged.csv"
# 100 


# get command line arguments 
args <- commandArgs(trailingOnly = TRUE)
ploidy_calls <- args[1]
metadata <- args[2]
# minimum number of patients in an age group for that age to be included 
min_age_count <- args[3]


# get functions `filter_data` and `day5_only`
setwd(".")
source("../../gwas/scripts/phenotypes/helper_functions/phenotyping_helper_functions.R")

# filter data 
ploidy_calls_filtered <- filter_data(ploidy_calls, parent, 
                                     bayes_factor_cutoff = 2, 
                                     nullisomy_threshold = 5, min_prob = 0.9) 
ploidy_calls_filtered <- day5_only(ploidy_calls_filtered, metadata)


# get embryo assignments in maternal age tranches 
aneuploidy_rate_by_age <- function(ploidy_calls, metadata, parent, 
                                   min_age_count) {
  
  # count chr of each embryo with maternal aneu, paternal aneu, or disomy
  embryo_counts <- ploidy_calls %>%
    group_by(child) %>%
    summarise(
      num_maternal = sum(bf_max_cat %in% c("1p", "3m")),
      num_paternal = sum(bf_max_cat %in% c("1m", "3p")),
      num_disomy = sum(bf_max_cat == "2")
    )
  # rename array to merge 
  colnames(embryo_counts)[1] <- "array"
  
  # count chr classifications for each embryo
  embryo_counts <- embryo_counts %>% 
    mutate(
      type = case_when(
        num_maternal > 0 & num_paternal == 0 ~ "maternal_aneu",
        num_maternal == 0 & num_paternal > 0 ~ "paternal_aneu",
        num_maternal == 0 & num_paternal == 0 & num_disomy > 15 ~ "euploid",
        num_maternal > 0 & num_paternal > 0 ~ "complex aneu",
        TRUE ~ "other"
      )
    )
  
  # merge aneuploidy counts and weighted ages dataframes
  merged_table <- merge(embryo_counts, 
                        metadata[,c("array", "patient_age", "partner_age")], 
                        by = "array", all.x = TRUE)
  
  # remove donors because their ages are not what's recorded 
  if (parent == "mother") {
    merged_table$rounded_age <- round(merged_table$patient_age)
  } else if (parent == "father") {
    merged_table$rounded_age <- round(merged_table$partner_age)
  }

  # remove patients that have NA for age 
  merged_table <- merged_table[complete.cases(merged_table$rounded_age),]
  
  # get aneuploidy rates for each age 
  aneu_rates <- merged_table %>%
    group_by(rounded_age) %>%
    summarise(
      mat_aneu_rate = sum(num_maternal) / 
        sum(num_maternal + num_paternal + num_disomy),
      pat_aneu_rate = sum(num_paternal) / 
        sum(num_maternal + num_paternal + num_disomy)
    )
  
  # keep only ages with at least x individuals in the dataset 
  # count number of occurrences of each rounded age
  frequency_counts <- merged_table %>%
    count(rounded_age)
  aneu_rates <- merge(aneu_rates, frequency_counts, by = "rounded_age")
  # keep only ages with at least x individuals in the dataset 
  aneu_rates <- aneu_rates %>% 
    filter(rounded_age %in% 
             frequency_counts$rounded_age[frequency_counts$n >= min_age_count])

  return(aneu_rates)
}

# plot error rates by maternal age 
error_rates_by_patient_age <- aneuploidy_rate_by_age(
  ploidy_calls_filtered, metadata, "mother", min_age_count)

ggplot(error_rates_by_patient_age, aes(x = rounded_age)) +
  geom_line(aes(y = mat_aneu_rate, color = "Maternal Error")) +
  geom_point(aes(y = mat_aneu_rate, color = "Maternal Error")) +
  geom_line(aes(y = pat_aneu_rate, color = "Paternal Error")) +
  geom_point(aes(y = pat_aneu_rate, color = "Paternal Error")) +
  labs(x = "Maternal Age", y = "Aneuploidy Ratio", color = "Error Type") +
  ggtitle("Maternal and Paternal Errors by Maternal Age") +
  scale_color_manual(values = c("Maternal Error" = "red", 
                                "Paternal Error" = "blue")) +
  theme_minimal() 


# plot error rates by paternal age 
error_rates_by_partner_age <- aneuploidy_rate_by_age(
  ploidy_calls_filtered, metadata, "father", min_age_count)
ggplot(error_rates_by_partner_age, aes(x = rounded_age)) +
  geom_line(aes(y = mat_aneu_rate, color = "Maternal Error")) +
  geom_point(aes(y = mat_aneu_rate, color = "Maternal Error")) +
  geom_line(aes(y = pat_aneu_rate, color = "Paternal Error")) +
  geom_point(aes(y = pat_aneu_rate, color = "Paternal Error")) +
  labs(x = "Maternal Age", y = "Aneuploidy Ratio", color = "Error Type") +
  ggtitle("Maternal and Paternal Errors by Paternal Age") +
  scale_color_manual(values = c("Maternal Error" = "red", 
                                "Paternal Error" = "blue")) +
  theme_minimal() 




# Plot ploidy ratios within each mother, by maternal age

# Get weighted maternal age across embryos, number of cycles,
# number of embryos, and number of euploid (all chr cn=2) and
# number of aneuploid (1 or more chr with other than cn=2)
count_embryo_error_types <- function(ploidy_calls, metadata, min_age_count) {
  
  # count chr of each embryo with maternal aneu, paternal aneu, or disomy
  embryo_counts <- ploidy_calls_filtered %>%
    group_by(mother, child) %>%
    summarise(
      num_maternal = sum(bf_max_cat %in% c("1p", "3m")),
      num_paternal = sum(bf_max_cat %in% c("1m", "3p")),
      num_disomy = sum(bf_max_cat == "2")
    ) 
  
  # count embryo ploidy per mother 
  mother_summary <- embryo_counts %>%
    group_by(mother) %>%
    summarise(
      maternal_error = sum(num_maternal > 0),
      paternal_error = sum(num_paternal > 0),
      euploid = sum(num_maternal == 0 & num_paternal == 0 & num_disomy >= 15)
    )
  
  # Create new column that tags every individual affiliated with each
  # parent, even if in different caseIDs
  metadata_merged_array <- metadata %>%
    mutate(
      array_id_merged = ifelse(family_position == parent,
                               paste0(array, "_merged"), NA_character_)
    ) %>%
    group_by(casefile_id) %>%
    fill(array_id_merged) %>%
    ungroup()
  
  # rename column to allow merge 
  colnames(mother_summary)[1] <- "array"
  
  # Create dataframe that calculates the weighted age, number of embryos,
  # and number of visits
  weighted_ages <- metadata_merged_array %>%
    filter(family_position == "child" & array %in% embryo_counts$child) %>%
    group_by(array_id_merged) %>%
    summarise(
      weighted_age = sum(patient_age) / n(),
      num_embryos = n(),
      num_visits = length(unique(patient_age))
    ) %>%
    as.data.frame()
  
  # Remove "_merged" from array column to allow easier downstream intersection
  weighted_ages$array <- gsub("_merged", "", weighted_ages$array_id_merged)
  # Remove array_id_merged column
  weighted_ages <- subset(weighted_ages, select = -c(array_id_merged))
  
  # Merge aneuploidy counts and weighted ages dataframes
  merged_table <- merge(weighted_ages, mother_summary, by = "array",
                        all.x = TRUE)
  
  # Remove any mothers who had zero embryos through karyoHMM 
  merged_table <- merged_table[(merged_table$maternal_error +
                                  merged_table$paternal_error + 
                                  merged_table$euploid) != 0, ]
  
  # Round maternal ages for plotting
  merged_table$rounded_age <- round(merged_table$weighted_age)
  
  # Calculate proportion aneuploid 
  merged_table$mat_error_ratio <- 
    merged_table$maternal_error / 
    (merged_table$maternal_error + merged_table$paternal_error + 
       merged_table$euploid)
  merged_table$pat_error_ratio <- 
    merged_table$paternal_error / 
    (merged_table$maternal_error + merged_table$paternal_error + 
       merged_table$euploid)
  
  # Remove any rows that are entirely NA 
  merged_table <- merged_table[complete.cases(merged_table$mat_error_ratio),]
  
  # Remove any mother for which no ratio could be calculated 
  # (removed entirely during embryo filtering for phenotypes) 
  # calculate average proportion aneu for each age 
  average_proportions <- merged_table %>%
    group_by(rounded_age) %>%
    summarise(avg_mat_error_ratio = mean(mat_error_ratio), 
              avg_pat_error_ratio = mean(pat_error_ratio))
  
  # keep only ages with at least x individuals in the dataset 
  # count number of occurrences of each rounded age
  frequency_counts <- merged_table %>%
    count(rounded_age)
  # keep only ages with at least x individuals in the dataset 
  average_proportions_filtered <- average_proportions %>% 
    filter(rounded_age %in% 
             frequency_counts$rounded_age[frequency_counts$n >= min_age_count])
  # remove individuals with NA for rounded age 
  # (due to NAs for patient or child in metadata)
  average_proportions_filtered <- 
    average_proportions_filtered[
      !is.na(average_proportions_filtered$rounded_age),]
  
  return(average_proportions_filtered)
}

errors_by_age <- count_embryo_error_types(ploidy_calls_filtered, 
                                          metadata, min_age_count)

ggplot(errors_by_age, aes(x = rounded_age)) +
  geom_line(aes(y = avg_mat_error_ratio, color = "Maternal Error")) +
  geom_point(aes(y = avg_mat_error_ratio, color = "Maternal Error")) +
  geom_line(aes(y = avg_pat_error_ratio, color = "Paternal Error")) +
  geom_point(aes(y = avg_pat_error_ratio, color = "Paternal Error")) +
  labs(x = "Weighted Age", y = "Aneuploidy Ratio by Mother", color = "Error Type") +
  ggtitle("Maternal and Paternal Errors by Maternal Age") +
  scale_color_manual(values = c("Maternal Error" = "red", "Paternal Error" = "blue")) +
  theme_minimal()

