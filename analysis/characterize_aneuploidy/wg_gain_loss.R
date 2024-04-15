# Whole genome gain or loss by parent 

# load libraries
library(data.table)
library(dplyr)
library(ggplot2) 

# load data 
ploidy_calls <- fread("/data/rmccoy22/natera_spectrum/karyohmm_outputs/compiled_output/natera_embryos.karyohmm_v30a.bph_sph_trisomy.full_annotation.031624.tsv.gz")

# filter data as per phenotyping 
#source("/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/phenotypes/generate_phenotypes.R")
# functions also below 
## Functions from source 
filter_data <- function(ploidy_calls, parent, bayes_factor_cutoff = 2,
                        nullisomy_threshold = 5, min_prob = 0.9) {
  
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


# filter data 
parent <- "mother"
metadata <- fread("/data/rmccoy22/natera_spectrum/data/summary_metadata/spectrum_metadata_merged.csv")
ploidy_calls_filtered <- filter_data(ploidy_calls, parent, bayes_factor_cutoff = 2,
                                     nullisomy_threshold = 5, min_prob = 0.9) 
ploidy_calls_filtered <- day5_only(ploidy_calls_filtered, metadata)

# calculate number of whole genome gain/loss by parent 
whole_genome <- function(ploidy_calls_filtered, parent, ploidy_threshold = 15) {
  
  # identify embryos with whole genome gain/loss 
  if (parent == "mother") {
    triploidy <- "3m"
    haploidy <- "1p"
  } else if (parent == "father") {
    triploidy <- "3p"
    haploidy <- "1m"
  }
  
  # find which individuals have whole genome gain/loss
  counts <- ploidy_calls_filtered %>%
    group_by(child) %>%
    summarise(num_affected_chromosomes_3 = sum(bf_max_cat == triploidy), 
              num_affected_chromosomes_1 = sum(bf_max_cat == haploidy)) %>% 
    ungroup()
  triploidies <- counts[counts$num_affected_chromosomes_3 >= ploidy_threshold,]
  haploidies <- counts[counts$num_affected_chromosomes_1 >= ploidy_threshold,]
  
  # count number of embryos with whole genome gain/loss
  num_triploidies <- nrow(triploidies)
  num_haploidies <- nrow(haploidies)
  
  output <- list(num_triploidies, num_haploidies)
  return(output)
}

maternal <- whole_genome(ploidy_calls_filtered, "mother")
paternal <- whole_genome(ploidy_calls_filtered, "father")


