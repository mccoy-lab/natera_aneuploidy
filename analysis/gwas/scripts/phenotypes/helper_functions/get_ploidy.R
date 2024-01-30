## Script called in each large-scale ploidy phenotype R file

# function to filter embryos by quality
filter_data <- function(embryos, bayes_factor_cutoff) {
  
  # confirm that bayes_factor_cutoff is numeric and positive
  if (!is.numeric(bayes_factor_cutoff) || bayes_factor_cutoff <= 0) {
    stop("Invalid 'bayes_factor_cutoff' argument. 
    	     It should be a positive numeric value.")
  }
  
  # remove embryos that have noise more than 3sd from mean
  embryos <- embryos[embryos$embryo_noise_3sd == FALSE,]
  
  # keep only rows that have probabilities for all 6 cn states
  embryos <- embryos[complete.cases(
    embryos[,c("0", "1m", "1p", "2", "3m", "3p")]),]
  
  # keep only rows that met the threshold for bayes factor qc
  embryos <- embryos[embryos$bf_max > bayes_factor_cutoff,]
  
  # keep each embryo only once 
  embryos <- embryos %>%
    group_by(child) %>%
    slice(1:22) %>%
    ungroup()
  
  return(embryos)
}

# function to keep only day 5 embryos 
day5_only <- function(embryos, metadata) {
  
  # intersect embryo with metadata (few_cells = day 5)
  day5 <- metadata[metadata$sample_scale == "few_cells",]
  day5_embryos <- embryos[embryos$child %in% day5$array,]
  
  return(day5_embryos)
}

# function to remove embryos that failed amplification
remove_failed_amp <- function(embryos, parent_col, nullisomy_threshold) {
  
  # count number of chromosomes called as nullisomies for each embryo 
  count_nullisomies <- embryos %>% 
    group_by({{parent_col}}, child) %>% 
    summarise(num_nullisomies = sum(bf_max_cat == "0"))
  # keep only embryos with fewer nullisomies than the threshold
  successful_amp <- count_nullisomies[count_nullisomies$num_nullisomies 
                                      < nullisomy_threshold,]
  
  # remove failed amplification embryos 
  embryos_filtered <- embryos[embryos$child %in% successful_amp$child,]
  
  # return remaining embryos
  return(embryos_filtered)
}

# function to identify and remove embryos with whole genome gain/loss 
remove_wholegenome_gainloss <- function(embryos, parent_col, ploidy_threshold) {
  
  # count number of single-chr aneuploidies in each embryo
  count_aneuploidies <- embryos %>% 
    group_by({{parent_col}}, child) %>% 
    summarise(num_aneuploidies = sum(bf_max_cat == "3m" | 
                                       bf_max_cat == "3p" |
                                       bf_max_cat == "1m" |
                                       bf_max_cat == "1p"))
  # keep only embryos with number of aneuploidies below the threshold 
  non_ploidy <- count_aneuploidies[count_aneuploidies$num_aneuploidies 
                                   < ploidy_threshold,]              
  
  # remove triploid, haploid embryos 
  embryos_filtered <- embryos[embryos$child %in% non_ploidy$child,]
  
  # return embryos with this info removed 
  return(embryos_filtered)
}

# count number of embryos per parent affected by each phenotype 
count_ploidy_by_parent <- function(embryos, parent_col, phenotype, 
                                   ploidy_threshold) {
  
  # check input 
  if (!(parent %in% c("mother", "father"))) {
    stop("Invalid 'parent' argument. Use 'mother' or 'father'.")
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
  } else {
    stop("Invalid 'phenotype' argument.")
  }
  
  # group ploidy by respective parent
  result <- embryos %>%
    group_by({{parent_col}}, child) %>%
    summarise(num_affected = sum(bf_max_cat %in% cn)) %>%
    mutate(is_ploidy = if_else(num_affected > 0, 
                               "aneu_true", "aneu_false")) %>%
    count(is_ploidy) %>%
    pivot_wider(names_from = is_ploidy, 
                values_from = n, values_fill = 0)
  
  return(result)
}

count_complex_ploidy_by_parent <- function(embryos, parent_col, parent) {
  # check input 
  if (!(parent %in% c("mother", "father"))) {
    stop("Invalid 'parent' argument. Use 'mother' or 'father'.")
  }
  
  # group ploidy by respective parent 
  result <- embryos %>%
    group_by({{parent}}, child) %>%
    summarise(parents_affected = 
                sum(any(bf_max_cat %in% c("3m", "1p")) 
                    & any(bf_max_cat %in% c("1m", "3p")))) %>% 
    mutate(is_ploidy = if_else(parents_affected > 1, 
                               "aneu_true", "aneu_false")) %>%
    count(is_ploidy) %>%
    pivot_wider(names_from = is_ploidy, 
                values_from = n, values_fill = 0) %>%
    replace(is.na(.), 0)
  
  return(result)
}

run_phenotype <- function(embryos, bayes_factor_cutoff, parent, 
                          phenotype, ploidy_threshold, metadata, nullisomy_threshold) {
  # filter embryo data 
  embryos_filtered <- filter_data(embryos, bayes_factor_cutoff)
  # keep only day 5 embryos
  embryos_filtered <- day5_only(embryos_filtered, metadata)
  # remove nullisomies (amplification errors)
  embryos_filtered <- remove_nullisomies(embryos_filtered, !!as.name(parent),
                                         nullisomy_threshold)
  
  # if not considering a whole genome gain/loss phenotype (haploid, triploid),
  # remove whole genome gain/loss embryos
  if (phenotype == "maternal_meiotic") {
    embryos_filtered <- remove_wholegenome_gainloss(embryos_filtered,
                                                    !!as.name(parent), 
                                                    ploidy_threshold)
  } 
  
  # group ploidy by respective parent 
  ploidy_counts_by_parent <- 
    count_ploidy_by_parent(embryos_filtered, !!as.name(parent), 
                           phenotype, ploidy_threshold, parent)
  colnames(ploidy_counts_by_parent)[1] <- "array"
  
  return(ploidy_counts_by_parent)
}
