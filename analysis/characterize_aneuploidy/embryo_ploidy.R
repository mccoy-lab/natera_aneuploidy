# Number of affected chr per embryo 
# Aneuploidy incidence (by type) per chr

# load libraries
library(data.table)
library(dplyr)
library(ggplot2) 

# load data 
ploidy_calls <- fread("/data/rmccoy22/natera_spectrum/karyohmm_outputs/compiled_output/natera_embryos.karyohmm_v30a.bph_sph_trisomy.full_annotation.031624.tsv.gz")

# filter data as per phenotyping 
source("/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/phenotypes/generate_phenotypes.R")
# functions also below 

# filter data (using functions from phenotyping and pasted below) 
parent <- "mother"
metadata <- fread("/data/rmccoy22/natera_spectrum/data/summary_metadata/spectrum_metadata_merged.csv")
ploidy_calls_filtered <- filter_data(ploidy_calls, parent, bayes_factor_cutoff = 2,
                                  nullisomy_threshold = 5, min_prob = 0.9) 
ploidy_calls_filtered <- day5_only(ploidy_calls_filtered, metadata)


# function to plot number of affected chr per embryo and return benchmarks 
affected_chr_per_embryo <- function(ploidy_calls_filtered, 
                                    max_meiotic = 3, 
                                    min_ploidy = 15) {
  
  # count affected chr per embryo 
  affected_chr_per_embryo <- ploidy_calls_filtered %>%
    group_by(child) %>%
    summarise(num_affected_chromosomes = sum(bf_max_cat != 2)) %>%
    ungroup()
  
  # plot count of affected chromosomes per embryo
  p1 <- ggplot(affected_chr_per_embryo, aes(x = num_affected_chromosomes)) +
    geom_bar() +
    labs(title = "Number of Affected Chromosomes per Embryo",
         x = "Number of Affected Chromosomes",
         y = "Count of Embryos") + 
    theme_minimal()
  
  # count number of single-chr affected
  num_meiotic <- length(which(affected_chr_per_embryo$num_affected_chromosomes <= max_meiotic))
  # count number of putative whole-genome affects 
  num_ploidy <- length(which(affected_chr_per_embryo$num_affected_chromosomes >= min_ploidy))
  
  # return output as list 
  output <- list(p1, num_meiotic, num_ploidy)
  return(output)
}

output <- affected_chr_per_embryo(ploidy_calls_filtered)
affected_chr_per_embryo <- output[1]
num_meiotic <- output[2]
num_ploidy <- output[3]



# function to plot error type by chromosome 
plot_ploidy_by_chr <- function(ploidy_calls_filtered, 
                               as_proportion = TRUE, 
                               aneuploid_only = TRUE, 
                               remove_whole_genome = TRUE, 
                               min_ploidy = 15) {
  # create new column with just numeric value of chromosome to improve plotting
  # remove the prefix "chr" from chromosome names
  ploidy_calls_filtered$chr <- gsub("chr", "", ploidy_calls_filtered$chrom)
  # convert chromosome names to factors and reorder them numerically
  ploidy_calls_filtered$chr <- factor(ploidy_calls_filtered$chr, levels = as.character(1:22))
  
  # remove embryos experiencing whole genome gain/loss 
  if (remove_whole_genome == TRUE) {
    affected_chr_per_embryo <- ploidy_calls_filtered %>%
      group_by(child) %>%
      summarise(num_affected_chromosomes = sum(bf_max_cat != 2)) %>%
      ungroup()
    whole_genome <- affected_chr_per_embryo[affected_chr_per_embryo$num_affected_chromosomes < min_ploidy,]
    ploidy_calls_filtered <- ploidy_calls_filtered[ploidy_calls_filtered$child %in% whole_genome$child,]
  }
  
  # remove euploid chromosomes 
  if (aneuploid_only == TRUE) {
    ploidy_calls_filtered <- ploidy_calls_filtered[ploidy_calls_filtered$bf_max_cat != 2,]
  }
  
  if (as_proportion == TRUE)
    # plot ploidy status by chr as proportion 
    p1 <- ggplot(ploidy_calls_filtered, aes(x = chr, fill = factor(bf_max_cat))) + 
      geom_bar(position = "fill") +
      labs(title = "Ploidy Status by Chromosome",
           x = "Chromosome",
           y = "Count") +
    theme_minimal()
  else {
    # plot ploidy status by chr as count 
    p1 <- ggplot(ploidy_calls_filtered, aes(x = chr, fill = factor(bf_max_cat))) +
      geom_bar(stat = "count") +
      labs(title = "Ploidy Status by Chromosome",
           x = "Chromosome",
           y = "Count") + 
      theme_minimal()
  }
  
  return(p1)
}

plot <- plot_ploidy_by_chr(ploidy_calls_filtered)



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
