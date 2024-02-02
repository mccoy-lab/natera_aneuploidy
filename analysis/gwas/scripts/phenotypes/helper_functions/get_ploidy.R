## Functions for use in generating aneuploidy phenotypes 
# /scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/scripts/phenotypes/aneuploidy_phenotypes.R

# function to filter embryos by quality
filter_data <- function(embryos, parent, bayes_factor_cutoff = 2,
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
  
  # keep each embryo only once 
  # consider grouping by mother, child 
  embryos <- embryos %>%
    group_by(child) %>%
    slice(1:22) %>%
    ungroup()
  
  # remove embryos that have noise more than 3sd from mean
  embryos <- embryos[embryos$embryo_noise_3sd == FALSE,]
  
  # remove embryos with failed amplification
  # count number of chromosomes called as nullisomies for each embryo 
  count_nullisomies <- embryos %>% 
    group_by({{parent}}, child) %>% 
    summarise(num_nullisomies = sum(bf_max_cat == "0"))
  # identify embryos with fewer nullisomies than the threshold
  successful_amp <- count_nullisomies[count_nullisomies$num_nullisomies 
                                      < nullisomy_threshold,]
  # keep only embryos without failed amplification 
  embryos <- embryos[embryos$child %in% successful_amp$child,]
  
  # keep only chrom that have probabilities for all 6 cn states
  embryos <- embryos[complete.cases(
    embryos[,c("0", "1m", "1p", "2", "3m", "3p")]),]
  
  # keep only rows that met the threshold for bayes factor qc
  embryos <- embryos[embryos$bf_max > bayes_factor_cutoff,]
  
  # add column that checks whether the max posterior is greater than threshold
  embryos <- embryos %>%
    mutate(high_prob = pmax(`0`, `1m`, `1p`, `2`, `3m`, `3p`) > min_prob)
  # keep only chromosomes with sufficiently high probability cn call
  embryos <- embryos[embryos$high_prob == TRUE,]
  
  return(embryos)
}

# function to keep only day 5 embryos 
day5_only <- function(embryos, metadata) {
  
  # intersect embryo with metadata (few_cells = day 5)
  day5_embryos <- embryos[embryos$child %in% metadata[metadata$sample_scale == 
                                                        "few_cells", ]$array,]
  return(day5_embryos)
}

# function to identify and remove embryos with whole genome gain/loss 
remove_wholegenome_gainloss <- function(embryos, parent, 
                                        max_meiotic = 3) {
  
  # count number of single-chr aneuploidies in each embryo
  count_aneuploidies <- embryos %>% 
    group_by({{parent}}, child) %>% 
    summarise(num_aneuploidies = sum(bf_max_cat == "3m" | 
                                       bf_max_cat == "3p" |
                                       bf_max_cat == "1m" |
                                       bf_max_cat == "1p"))
  # keep only embryos with number of aneuploidies below the threshold 
  non_ploid <- count_aneuploidies[count_aneuploidies$num_aneuploidies 
                                   < max_meiotic,]              
  
  # remove embryos affected by multiple aneuploidies 
  embryos_filtered <- embryos[embryos$child %in% non_ploid$child,]
  
  return(embryos_filtered)
}

# count number of embryos per parent affected by each phenotype 
count_ploidy_by_parent <- function(embryos, parent, phenotype, 
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
  if (phenotype == "maternal_triploidy") {
    cn <- "3m"
  } else if (phenotype == "paternal_triploidy") {
    cn <- "3p"
  } else if (phenotype == "maternal_haploidy") {
    cn <- "1p"
  } else if (phenotype == "paternal_haploidy") {
    cn <- "1m"
  } else if (phenotype == "maternal_meiotic") {
    cn <- c("3m", "1p")
  } else if (phenotype == "complex_aneuploidy") {
    cn <- c("0", "1m", "1p", "3m", "3p")
  } else {
    stop("Invalid 'phenotype' argument.")
  }
  
  # add that complex aneu includes at least two unique types of aneu 
  # Count ploidies based on phenotype definition 
  result <- embryos_filtered %>%
    group_by({{parent}}, child) %>%
    summarise(num_affected = sum(bf_max_cat %in% cn),
              unique_bf_max_cat = n_distinct(bf_max_cat[bf_max_cat %in% cn])) %>%
    mutate(
      is_ploidy = case_when(
        phenotype == "maternal_meiotic" ~ ifelse(
          num_affected > 0 & num_affected < max_meiotic,
          "aneu_true", "aneu_false"),
        phenotype == "complex_aneuploidy" ~ ifelse(
          num_affected > 0  & unique_bf_max_cat >= 2, 
          "aneu_true", "aneu_false"),
        TRUE ~ ifelse(num_affected > min_ploidy, "aneu_true", "aneu_false")
      )
    ) %>%
    count(is_ploidy) %>%
    pivot_wider(names_from = is_ploidy, values_from = n, values_fill = 0)
  
  # match column names with external data 
  colnames(result)[1] <- "array"

  return(result)
}


# generate the phenotype of choice 
run_phenotype <- function(embryos, parent, metadata, phenotype,
                          filter_day_5 = TRUE, bayes_factor_cutoff = 2, 
                          nullisomy_threshold = 5, max_meiotic = 3,
                          min_ploidy = 15, min_prob = 0.9) {
  
  # filter embryo data 
  embryos_filtered <- filter_data(embryos, !!as.name(parent),
                                  bayes_factor_cutoff, nullisomy_threshold,
                                  min_prob)
  
  # keep only day 5 embryos
  if (filter_day_5 == TRUE) {
    embryos_filtered <- day5_only(embryos_filtered, metadata)
  }
  
  # if not considering a whole genome gain/loss phenotype (haploid, triploid),
  # remove whole genome gain/loss embryos
  if (phenotype == "maternal_meiotic") {
    embryos_filtered <- remove_wholegenome_gainloss(embryos_filtered,
                                                    !!as.name(parent), 
                                                    max_meiotic)
  } 
  
  # group ploidy by respective parent 
  ploidy_counts_by_parent <- count_ploidy_by_parent(embryos_filtered, 
                                                    !!as.name(parent),
                                                    phenotype, max_meiotic)
  
  return(ploidy_counts_by_parent)
}
