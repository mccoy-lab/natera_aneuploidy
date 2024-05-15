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

# if running locally 
ploidy_calls <- fread("/data/rmccoy22/natera_spectrum/karyohmm_outputs/compiled_output/natera_embryos.karyohmm_v30a.bph_sph_trisomy.full_annotation.031624.tsv.gz")
metadata <- fread("/data/rmccoy22/natera_spectrum/data/summary_metadata/spectrum_metadata_merged.csv")
min_age_count <- 20

# filter data 
ploidy_calls_filtered <- filter_data(ploidy_calls, parent, 
                                     bayes_factor_cutoff = 2, 
                                     nullisomy_threshold = 5, min_prob = 0.9) 
ploidy_calls_filtered <- day5_only(ploidy_calls_filtered, metadata)



# Plot ploidy ratios within each mother, by maternal age

# Get weighted maternal age across embryos, number of cycles,
# number of embryos, and number of euploid (all chr cn=2) and
# number of aneuploid (1 or more chr with other than cn=2)

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
    partner_age_weighted = sum(partner_age) / n(),
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

# Remove any mothers who had zero embryos with karyoHMM calls 
merged_table <- merged_table[(merged_table$maternal_error +
                                merged_table$paternal_error + 
                                merged_table$euploid) != 0, ]

# Remove egg donors (age doesn't show genetic age of donor)
merged_table <- merged_table[!(merged_table$array %in%
                                 metadata[metadata$egg_donor == "yes",]$array),]

# Round maternal ages for plotting
merged_table$rounded_age <- round(merged_table$weighted_age)
merged_table$rounded_partner_age <- round(merged_table$partner_age_weighted)

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

# keep only ages with at least x individuals in the dataset 
# count number of occurrences of each rounded age
frequency_counts <- merged_table %>%
  count(rounded_age)
# keep only ages with at least x individuals in the dataset 
merged_table <- merged_table %>% 
  filter(rounded_age %in% 
           frequency_counts$rounded_age[frequency_counts$n >= min_age_count])


# get se for each parent  
avg_mat_error <- merged_table %>%
  group_by(rounded_age) %>%
  summarise(avg_mat_error = mean(mat_error_ratio),
            se_mat = sd(mat_error_ratio) / sqrt(n()),
            avg_pat_error = mean(pat_error_ratio),
            se_pat = sd(pat_error_ratio) / sqrt(n()))
frequency_counts_paternal <- merged_table %>%
  count(rounded_partner_age)
avg_pat_error <- merged_table %>%
  filter(rounded_partner_age %in% frequency_counts_paternal$rounded_partner_age[frequency_counts_paternal$n >= min_age_count]) %>%
  group_by(rounded_partner_age) %>%
  summarise(avg_pat_error = mean(pat_error_ratio),
            se_pat = sd(pat_error_ratio) / sqrt(n()))



# Make glm data for maternal 
merged_table$not_maternal_error <- merged_table$paternal_error + merged_table$euploid

m0 <- glm(data = merged_table, formula = cbind(maternal_error, not_maternal_error) ~ weighted_age, family = "quasibinomial")
m1 <- glm(data = merged_table, formula = cbind(maternal_error, not_maternal_error) ~ poly(weighted_age, 2), family = "quasibinomial")
anova(m1, m0, test = "Chisq")

predicted_data <- data.table(weighted_age = 26:45)
predicted_data[, avg_mat_error := predict(m1, newdata = predicted_data, type = "response")]
predicted_data[, rounded_age := weighted_age]


# Plot maternal points 
ggplot(avg_mat_error, aes(x = rounded_age)) +
  geom_point(aes(y = avg_mat_error), color = "#CC79A7") +
  theme_minimal() + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank()) + 
  ylim(0, 1)
  
# Plot maternal points with glm 
ggplot(avg_mat_error, aes(x = rounded_age)) +
  geom_point(aes(y = avg_mat_error), color = "#CC79A7") +
  geom_line(data = predicted_data, aes(y = avg_mat_error), color = "#CC79A7") +
  theme_minimal() + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank()) + 
  ylim(0, 1)

# Plot maternal points with glm and paternal points 
ggplot(avg_mat_error, aes(x = rounded_age)) +
  geom_point(aes(y = avg_mat_error), color = "#CC79A7") +
  geom_line(data = predicted_data, aes(y = avg_mat_error), color = "#CC79A7") +
  geom_point(aes(y = avg_pat_error), color = "#0072B2") +
  theme_minimal() + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank()) + 
  ylim(0, 1)



# Make glm data for paternal 
merged_table_pat <- merged_table[!is.na(merged_table$partner_age_weighted),]
merged_table_pat$not_paternal_error <- merged_table_pat$maternal_error + merged_table_pat$euploid

m2 <- glm(data = merged_table_pat, formula = cbind(paternal_error, not_paternal_error) ~ partner_age_weighted, family = "quasibinomial")
m3 <- glm(data = merged_table_pat, formula = cbind(paternal_error, not_paternal_error) ~ poly(partner_age_weighted, 2), family = "quasibinomial")
anova(m2, m3, test = "Chisq")

predicted_data_pat <- data.table(partner_age_weighted = 28:53)
predicted_data_pat[, avg_pat_error := predict(m2, newdata = predicted_data_pat, type = "response")]
predicted_data_pat[, rounded_partner_age := partner_age_weighted]

# Plot paternal points 
ggplot(avg_pat_error, aes(x = rounded_partner_age)) +
  geom_point(aes(y = avg_pat_error), color = "#0072B2") +
  theme_minimal() + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank()) + 
  ylim(0, 1)

# Plot paternal points with glm 
ggplot(avg_pat_error, aes(x = rounded_partner_age)) +
  geom_point(aes(y = avg_pat_error), color = "#0072B2") +
  geom_line(data = predicted_data_pat, aes(y = avg_pat_error), color = "#0072B2") +
  theme_minimal() + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank()) + 
  ylim(0, 1)
