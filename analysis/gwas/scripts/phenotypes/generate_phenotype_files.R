## Make phenotype file for embryos affected by aneuploidy phenotypes

# load libraries
library(data.table)
library(tidyr)
library(dplyr)

# Usage:
# /scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/scripts/phenotypes/generate_phenotype_files.R \
# /data/rmccoy22/natera_spectrum/karyohmm_outputs/compiled_output/natera_embryos.karyohmm_v30a.bph_sph_trisomy.full_annotation.031624.tsv.gz \
# /scratch16/rmccoy22/abiddan1/natera_segmental/analysis/segmental_qc/results/tables/segmental_calls_postqc.tsv.gz \ # segmental aneuploidy calls to remove chrom 
# mother \ # parent to group phenotype by
# /data/rmccoy22/natera_spectrum/data/summary_metadata/spectrum_metadata_merged.csv \
# maternal_meiotic_aneuploidy \ # phenotype name
# TRUE
# 2 \ remove chr with bayes factor > bayes_factor_cutoff
# 5 \ remove embryos that had more chr with cn = 0 for than nullisomy_threshold
# 0.9 \ minimum posterior probability for each cn call
# 5 \ max number of affected chr to count for maternal meiotic phenotype
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
                        bayes_factor_cutoff = 2, filter_day_5 = TRUE, 
                        nullisomy_threshold = 5, 
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
  
  # Remove day 3 embryos 
  if (filter_day_5 == TRUE) {
    ploidy_calls <- ploidy_calls[ploidy_calls$day3_embryo == FALSE,]
  }
  
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

# create phenotypes for trais, both aneuploidy and metadata-based 
# (maternal meiotic aneuploidy, triploidy, haploidy, embryo count, maternal age)  
make_phenotype <- function(metadata, parent, phenotype, ploidy_calls, 
                           segmental_calls, bayes_factor_cutoff = 2, 
                           filter_day_5 = TRUE, nullisomy_threshold = 5, 
                           max_meiotic = 5, min_ploidy = 15) {
  
  # assign all members of each family to the same mother, across casefile IDs 
  metadata_mothers <- metadata %>%
    mutate(mother_id = ifelse(family_position == "mother", array, NA)) %>%
    fill(mother_id, .direction = "downup")
  
  # assign visit and keep only day 5 embryos
  child_data <- metadata_mothers %>%
    filter(family_position == 'child', sample_scale == 'few_cells') %>%
    arrange(mother_id, patient_age) %>%
    group_by(mother_id) %>%
    mutate(visit_id = dense_rank(patient_age)) %>%
    ungroup()
  
  # count number of visits per mother 
  num_visits <- child_data %>%
    group_by(mother_id) %>%
    summarise(num_visits = n_distinct(visit_id)) %>%
    ungroup()
  
  # for aneuploidy phenotypes, count aneuploid/euploid embryos per visit
  # and merge with parental age 
  if (grepl("ploidy", phenotype)) {
    # select ploidy status based on phenotype and parent 
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
    
    # filter karyohmm data (quality control, day 5, posterior probabilities)
    ploidy_calls <- filter_data(ploidy_calls, parent, segmental_calls, 
                                bayes_factor_cutoff, filter_day_5,
                                nullisomy_threshold, min_prob)
    
    # create a lookup table for array and visit_id
    visit_lookup <- child_data %>%
      dplyr::select(array, visit_id) %>%
      distinct()
    
    # add visit number to each embryo
    ploidy_calls <- ploidy_calls %>%
      left_join(visit_lookup, by = c("mother" = "array")) %>%
      rename(visit_id_mother = visit_id) %>%
      left_join(visit_lookup, by = c("father" = "array")) %>%
      rename(visit_id_father = visit_id) %>%
      left_join(visit_lookup, by = c("child" = "array")) %>%
      rename(visit_id_child = visit_id)
    
    # combine the visit_id columns across family members
    ploidy_calls <- ploidy_calls %>%
      mutate(visit_id = 
               coalesce(visit_id_mother, visit_id_father, visit_id_child)) %>%
      dplyr::select(-visit_id_mother, -visit_id_father, -visit_id_child)
    
    # count number of aneuploid and euploid embryos in each visit 
    result <- ploidy_calls %>%
      group_by(mother, father, child, visit_id) %>%
      summarise(num_affected = sum(bf_max_cat %in% cn), unique_bf_max_cat = 
                  n_distinct(bf_max_cat[bf_max_cat %in% cn])) %>%
      mutate(
        is_ploidy = case_when(
          phenotype == "maternal_meiotic_aneuploidy" ~ 
            ifelse(num_affected > 0 & num_affected <= max_meiotic, 
                   "aneu_true", "aneu_false"),
          phenotype == "complex_aneuploidy" ~ 
            ifelse(num_affected > 0 & unique_bf_max_cat >= 2, 
                   "aneu_true", "aneu_false"),
          TRUE ~ 
            ifelse(num_affected > min_ploidy, "aneu_true", "aneu_false")
        )
      ) %>%
      count(visit_id, is_ploidy) %>%
      pivot_wider(names_from = is_ploidy, values_from = n, values_fill = list(n = 0))
    
    # track parental age at each visit 
    age_produced <- child_data %>%
      group_by(mother_id, visit_id) %>%
      summarise(patient_age = mean(patient_age, na.rm = TRUE),
                partner_age = mean(partner_age, na.rm = TRUE)) %>%
      ungroup()
    
    # group by family 
    result2 <- result %>%
      group_by(mother, father, visit_id) %>%
      summarise(
        aneu_true = sum(aneu_true),     
        aneu_false = sum(aneu_false),   
        total_embryos = sum(aneu_true + aneu_false)  
      ) %>%
      left_join(age_produced, by = c("mother" = "mother_id", "visit_id")) %>%
      left_join(num_visits, by = c("mother" = "mother_id")) %>%
      ungroup()
  } else {
    # calculate phenotypes based on metadata (embryo count, maternal age) 
    mother_summary <- child_data %>%
      group_by(mother_id, visit_id) %>%
      summarise(
        num_embryos = n(),                     # Count number of children per visit
        patient_age = first(patient_age),
        partner_age = first(partner_age) # Capture the patient age for the visit
      ) %>%
      ungroup()
    
    # count number of visits per mother 
    mother_summary <- mother_summary %>%
      left_join(num_visits, by = "mother_id")
  }
}

# Generate file for phenotype of interest
run_phenotype <- function(ploidy_calls, parent, segmental_calls, metadata, 
                          phenotype, filter_day_5 = TRUE, 
                          bayes_factor_cutoff = 2, nullisomy_threshold = 5, 
                          min_prob = 0.9, max_meiotic = 5, min_ploidy = 15) {
  
  # Filter embryo data
  ploidy_calls <- filter_data(ploidy_calls, parent, segmental_calls, 
                              bayes_factor_cutoff, filter_day_5, 
                              nullisomy_threshold, min_prob)
  
  # Compute phenotype 
  pheno_output <- make_phenotype(metadata, parent, phenotype, ploidy_calls, 
                                 max_meiotic, min_ploidy)
  
  return(pheno_output)
}


# Read in embryos and metadata
ploidy_calls <- fread(ploidy_calls)
segmental_calls <- fread(segmental_calls)
metadata <- fread(metadata)

# Run phenotype 
pheno_by_parent <- make_phenotype(metadata, parent, phenotype, ploidy_calls, 
                                  segmental_calls, bayes_factor_cutoff = 2, 
                                  filter_day_5 = TRUE, nullisomy_threshold = 5, 
                                  max_meiotic = 5, min_ploidy = 15)

# Write phenotype info to file
write.csv(pheno_by_parent, out_fname, row.names = FALSE)
