## Make phenotype file for embryos affected by aneuploidy phenotypes

# load libraries
library(data.table)
library(tidyr)
library(dplyr)

# Usage:
# /scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/scripts/phenotypes/aneuploidy_phenotypes.R \
# /data/rmccoy22/natera_spectrum/karyohmm_outputs/compiled_output/natera_embryos.karyohmm_v30a.bph_sph_trisomy.full_annotation.031624.tsv.gz \
# mother \ # parent to group phenotype by
# /scratch16/rmccoy22/abiddan1/natera_segmental/analysis/segmental_qc/results/tables/segmental_calls_postqc.tsv.gz \ # segmental aneuploidy calls to remove chrom 
# /data/rmccoy22/natera_spectrum/data/summary_metadata/spectrum_metadata_merged.csv \
# maternal_meiotic_aneuploidy \ # phenotype name
# TRUE
# 2 \ remove chr with bayes factor > bayes_factor_cutoff
# 5 \ remove embryos that had more chr with cn = 0 for than nullisomy_threshold
# 0.9 \ minimum posterior probability for each cn call
# 3 \ max number of affected chr to count for maternal meiotic phenotype
# 15 \ min number of affected chr to count for whole genome gain/loss
# /scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/results/phenotypes/maternal_meiotic_aneuploidy_by_mother.csv \

# get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# assign arguments to variables
# ploidy calls for each embryo's chromosomes, from karyoHMM
ploidy_calls <- args[1]
# segmental calls to distinguish from whole chromosomal aneuploidies 
segmental_calls <- args[2]
# parent to measure phenotype
parent <- args[3]
# metadata to filter for day5 embryos
metadata <- args[4]
# phenotype
phenotype <- args[5]
# whether to keep only day 5 embryos
filter_day_5 <- as.logical(args[6])
# minimum bayes factor for filtering
bayes_factor_cutoff <- as.numeric(args[7])
# how many cn = 0 before it's failed amplification
nullisomy_threshold <- as.numeric(args[8])
# min posterior probability of cn state
min_prob <- as.numeric(args[9])
# max number of chr for maternal meiotic pheno
max_meiotic <- as.numeric(args[10])
# min number of chr for whole genome gain/loss
min_ploidy <- as.numeric(args[11])
# output file name
out_fname <- args[12]


# Function to filter embryos by quality
filter_data <- function(ploidy_calls, parent, segmental_calls, 
                        bayes_factor_cutoff = 2, nullisomy_threshold = 5, 
                        min_prob = 0.9) {
  
  # Confirm that bayes_factor_cutoff is numeric and positive
  if (!is.numeric(bayes_factor_cutoff) || bayes_factor_cutoff <= 0) {
    stop("Invalid 'bayes_factor_cutoff' argument.
    	     It should be a positive numeric value.")
  }
  
  # Confirm that nullisomy_threshold is numeric and positive
  if (!is.numeric(nullisomy_threshold) || nullisomy_threshold <= 0) {
    stop("Invalid 'nullisomy_threshold' argument.
    	     It should be a positive numeric value.")
  }
  
  # Confirm that min_prob is numeric and positive
  if (!is.numeric(min_prob) || min_prob <= 0 || min_prob > 1) {
    stop("Invalid 'min_prob' argument.
    	     It should be a positive numeric value.")
  }
  
  # Keep each embryo only once (one call for each chromosome)
  ploidy_calls <- ploidy_calls %>%
    distinct(child, chrom, .keep_all = TRUE)

  
  # Remove embryos that have noise more than 3sd from mean
  ploidy_calls <- ploidy_calls[ploidy_calls$embryo_noise_3sd == FALSE, ]
  
  # Remove embryos with failed amplification
  # Count number of chromosomes called as nullisomies for each embryo
  count_nullisomies <- ploidy_calls %>%
    group_by(get(parent), child) %>%
    summarise(num_nullisomies = sum(bf_max_cat == "0"))
  # Identify embryos with fewer nullisomies than the threshold
  successful_amp <- count_nullisomies[count_nullisomies$num_nullisomies
                                      < nullisomy_threshold, ]
  # Keep only embryos without failed amplification
  ploidy_calls <- ploidy_calls[ploidy_calls$child %in% successful_amp$child, ]
  
  # Keep only chrom that have probabilities for all 6 cn states
  ploidy_calls <- ploidy_calls[complete.cases(
    ploidy_calls[,c("0", "1m", "1p", "2", "3m", "3p")]), ]
  
  # Keep only rows that met the threshold for bayes factor qc
  ploidy_calls <- ploidy_calls[ploidy_calls$bf_max > bayes_factor_cutoff, ]
  
  # Add column that checks whether the max posterior is greater than threshold
  ploidy_calls <- ploidy_calls %>%
    mutate(high_prob = pmax(`0`, `1m`, `1p`, `2`, `3m`, `3p`) > min_prob)
  # Keep only chromosomes with sufficiently high probability cn call
  ploidy_calls <- ploidy_calls[ploidy_calls$high_prob == TRUE, ]
  
  # Keep only chromosomes that are not affected by post-QC segmental aneu
  # Add column that makes ID of mother-father-child-chrom in whole chr 
  ploidy_calls$uid <- paste0(ploidy_calls$mother, "+", ploidy_calls$father, "+",
                             ploidy_calls$child, "+", ploidy_calls$chrom)
  # Add column that makes ID of mother-father-child-chrom in segmental
  segmental_calls$uid <- paste0(segmental_calls$mother, "+", 
                                segmental_calls$father, "+", 
                                segmental_calls$child, "+", 
                                segmental_calls$chrom)
  # Remove from ploidy calls any chromosomes in segmentals 
  ploidy_calls <- ploidy_calls[!ploidy_calls$uid %in% segmental_calls$uid,]
  
  return(ploidy_calls)
}

# Function to keep only day 5 embryos
day5_only <- function(ploidy_calls, metadata) {
  
  # Intersect embryo with metadata (few_cells = day 5)
  ploidy_calls <- ploidy_calls[ploidy_calls$child %in% 
                                 metadata[metadata$sample_scale ==
                                            "few_cells", ]$array, ]
  return(ploidy_calls)
}

# Create table with array, number of embryos, number of visits, weighted age,
# and aneuploidy counts (if ploidy phenotype) 
make_phenotype <- function(metadata, parent, phenotype, ploidy_calls, 
                           max_meiotic = 3,
                           min_ploidy = 15) {
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
  
  # Create dataframe that calculates the weighted age, number of embryos,
  # and number of visits
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
  weighted_ages$array <- gsub("_merged", "", weighted_ages$array_id_merged)
  # Remove array_id_merged column
  weighted_ages <- subset(weighted_ages, select = -c(array_id_merged))
  # Rename for return if not ploidy merging 
  return <- weighted_ages 
  
  # If phenotype is a ploidy, count ploidies by parent and merge with age 
  if (grepl("ploidy", phenotype)) {
    # Select ploidy status based on phenotype and parent 
    if (phenotype == "triploidy" & parent == "mother") {
      cn <- "3m"
    } else if (phenotype == "triploidy" & parent == "father") {
      cn <- "3p"
    } else if (phenotype == "haploidy" & parent == "mother") {
      cn <- "1p"
    } else if (phenotype == "haploidy" & parent == "father") {
      cn <- "1m"
    } else if (phenotype == "maternal_meiotic_aneuploidy") {
      cn <- c("3m", "1p")
    } else if (phenotype == "complex_aneuploidy") {
      cn <- c("0", "1m", "1p", "3m", "3p")
    }
    
    # Count ploidies based on phenotype definition
    result <- ploidy_calls %>%
      group_by(get(parent), child) %>%
      summarise(num_affected = sum(bf_max_cat %in% cn),
                unique_bf_max_cat =
                  n_distinct(bf_max_cat[bf_max_cat %in% cn])) %>%
      mutate(
        is_ploidy = case_when(
          phenotype == "maternal_meiotic_aneuploidy" ~ 
            ifelse(num_affected > 0 & num_affected < max_meiotic,
                   "aneu_true", "aneu_false"),
          phenotype == "complex_aneuploidy" ~ ifelse(num_affected > 0
                                                     & unique_bf_max_cat >= 2,
                                                     "aneu_true", "aneu_false"),
          TRUE ~ ifelse(num_affected > min_ploidy, "aneu_true", "aneu_false")
        )
      ) %>%
      count(is_ploidy) %>%
      pivot_wider(names_from = is_ploidy, values_from = n, values_fill = 0)
    
    # Update column name to match with external data
    colnames(result)[1] <- "array"
    
    # Merge aneuploidy counts and weighted ages dataframes
    return <- merge(weighted_ages, result, by = "array",
                          all.x = TRUE)
    
  }
  
  return(return)
}

# Generate file for phenotype of interest
run_phenotype <- function(ploidy_calls, parent, segmental_calls, metadata, 
                          phenotype, filter_day_5 = TRUE, 
                          bayes_factor_cutoff = 2, nullisomy_threshold = 5, 
                          min_prob = 0.9, max_meiotic = 3, min_ploidy = 15) {
  
  # Filter embryo data
  ploidy_calls <- filter_data(ploidy_calls, parent, segmental_calls, 
                              bayes_factor_cutoff, nullisomy_threshold, 
                              min_prob)
  
  # Keep only day 5 embryos
  if (filter_day_5 == TRUE) {
    ploidy_calls <- day5_only(ploidy_calls, metadata)
  }
  
  # Compute phenotype 
  pheno_output <- make_phenotype(metadata, parent, phenotype, ploidy_calls, 
                                 max_meiotic, min_ploidy)
  
  return(pheno_output)
}


# Read in embryos and metadata
ploidy_calls <- fread(ploidy_calls)
segmental_calls <- fread(segmental_calls)
metadata <- fread(metadata)

# Keep only high-quality embryos (remove noisy, high-bayes factor, and
# Potential failed amplification embryos) and day 5 embryos
# Count number of aneuploid and non-aneuploid embryos per parent
pheno_by_parent <- run_phenotype(ploidy_calls, parent, segmental_calls, 
                                 metadata, phenotype, filter_day_5, 
                                 bayes_factor_cutoff, nullisomy_threshold, 
                                 min_prob, max_meiotic, min_ploidy)

# Write phenotype info to file
write.csv(pheno_by_parent, out_fname, row.names = FALSE)
