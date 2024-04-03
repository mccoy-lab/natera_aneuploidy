## Functions for use in generating aneuploidy phenotypes 
# /scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/scripts/phenotypes/aneuploidy_phenotypes.R

# function to filter embryos by quality
filter_data <- function(ploidy_calls, parent, bayes_factor_cutoff = 2,
                        nullisomy_threshold = 5, min_prob = 0.9) {
  
  # confirm that bayes_factor_cutoff is numeric and positive
  if (!is.numeric(bayes_factor_cutoff) || bayes_factor_cutoff <= 0) {
    stop("Invalid 'bayes_factor_cutoff' argument. 
    	     It should be a positive numeric value.")
  }
  
  # confirm that nullisomy_threshold is numeric and positive
  if (!is.numeric(nullisomy_threshold) || nullisomy_threshold <= 0) {
    stop("Invalid 'nullisomy_threshold' argument. 
    	     It should be a positive numeric value.")
  }
  
  # confirm that min_prob is numeric and positive
  if (!is.numeric(min_prob) || min_prob <= 0) {
    stop("Invalid 'min_prob' argument. 
    	     It should be a positive numeric value.")
  }
  
  # keep each embryo only once (one call for each chromosome)
  ploidy_calls <- ploidy_calls %>%
    distinct(child, chrom, .keep_all = TRUE)

  
  # remove embryos that have noise more than 3sd from mean
  ploidy_calls <- ploidy_calls[ploidy_calls$embryo_noise_3sd == FALSE,]
  
  # remove embryos with failed amplification
  # count number of chromosomes called as nullisomies for each embryo 
  count_nullisomies <- ploidy_calls %>% 
    group_by(get(parent), child) %>% 
    summarise(num_nullisomies = sum(bf_max_cat == "0"))
  # identify embryos with fewer nullisomies than the threshold
  successful_amp <- count_nullisomies[count_nullisomies$num_nullisomies 
                                      < nullisomy_threshold,]
  # keep only embryos without failed amplification 
  ploidy_calls <- ploidy_calls[ploidy_calls$child %in% successful_amp$child,]
  
  # keep only chrom that have probabilities for all 6 cn states
  ploidy_calls <- ploidy_calls[complete.cases(
    ploidy_calls[,c("0", "1m", "1p", "2", "3m", "3p")]),]
  
  # keep only rows that met the threshold for bayes factor qc
  ploidy_calls <- ploidy_calls[ploidy_calls$bf_max > bayes_factor_cutoff,]
  
  # add column that checks whether the max posterior is greater than threshold
  ploidy_calls <- ploidy_calls %>%
    mutate(high_prob = pmax(`0`, `1m`, `1p`, `2`, `3m`, `3p`) > min_prob)
  # keep only chromosomes with sufficiently high probability cn call
  ploidy_calls <- ploidy_calls[ploidy_calls$high_prob == TRUE,]
  
  return(ploidy_calls)
}

# function to keep only day 5 embryos 
day5_only <- function(ploidy_calls, metadata) {
  
  # intersect embryo with metadata (few_cells = day 5)
  ploidy_calls <- ploidy_calls[ploidy_calls$child %in% 
                                 metadata[metadata$sample_scale == 
                                            "few_cells", ]$array,]
  return(ploidy_calls)
}


# count number of embryos per parent affected by each phenotype 
count_ploidy_by_parent <- function(ploidy_calls, parent, phenotype, 
                                   max_meiotic = 3,
                                   min_ploidy = 15) {
  
  # confirm that max_meiotic is numeric and positive
  if (!is.numeric(max_meiotic) || max_meiotic <= 0) {
    stop("Invalid 'max_meiotic' argument. 
    	     It should be a positive numeric value.")
  }
  
  # confirm that min_ploidy is numeric and positive
  if (!is.numeric(min_ploidy) || min_ploidy <= 0) {
    stop("Invalid 'min_ploidy' argument. 
    	     It should be a positive numeric value.")
  }
  
  # choose relevant aneuploidies based on phenotype 
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
  } else {
    stop("Invalid 'phenotype' argument.")
  }
  
  # Count ploidies based on phenotype definition 
  result <- ploidy_calls %>%
    group_by(get(parent), child) %>%
    summarise(num_affected = sum(bf_max_cat %in% cn), 
              unique_bf_max_cat = 
                n_distinct(bf_max_cat[bf_max_cat %in% cn])) %>% 
    mutate(
      is_ploidy = case_when(
        phenotype == "maternal_meiotic_aneuploidy" ~ ifelse(
          num_affected > 0 & num_affected < max_meiotic,
          "aneu_true", "aneu_false"),
        phenotype == "complex_aneuploidy" ~ ifelse(
          num_affected > 0  & unique_bf_max_cat >= 2, 
          "aneu_true", "aneu_false"),
        TRUE ~ ifelse(num_affected > min_ploidy, 
          "aneu_true", "aneu_false")
      )
    ) %>%
    count(is_ploidy) %>%
    pivot_wider(names_from = is_ploidy, values_from = n, values_fill = 0)
  
  # match column names with external data 
  colnames(result)[1] <- "array"

  return(result)
}

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


# generate the phenotype of interest 
run_phenotype <- function(ploidy_calls, parent, metadata, phenotype,
                          filter_day_5 = TRUE, bayes_factor_cutoff = 2, 
                          nullisomy_threshold = 5, min_prob = 0.9,
                          max_meiotic = 3, min_ploidy = 15) {
  
  # filter embryo data 
  ploidy_calls <- filter_data(ploidy_calls, parent,
                                  bayes_factor_cutoff, nullisomy_threshold,
                                  min_prob)
  
  # keep only day 5 embryos
  if (filter_day_5 == TRUE) {
    ploidy_calls <- day5_only(ploidy_calls, metadata)
  }
  
  # if aneuploidy phenotype, group ploidy by respective parent 
  # else if embryo count phenotype, count embryos (total, euploid, and aneu) 
  # by respective parent 
  if (grepl("ploidy", phenotype)) {
    pheno_output <- count_ploidy_by_parent(ploidy_calls, parent, 
                                           phenotype, max_meiotic)
  } else {
    pheno_output <- count_embryos_by_parent(ploidy_calls, metadata, parent)
  }
  
  return(pheno_output)
}
