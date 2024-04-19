# Number of affected chr per embryo 
# Aneuploidy incidence (by type) per chr

# load libraries
library(data.table)
library(dplyr)
library(ggplot2) 

# Usage: /aneuploidy_post/utils/embryo_ploidy.R \ 
# "/data/rmccoy22/natera_spectrum/karyohmm_outputs/compiled_output/natera_embryos.karyohmm_v30a.bph_sph_trisomy.full_annotation.031624.tsv.gz" \
# "mother"
# "/data/rmccoy22/natera_spectrum/data/summary_metadata/spectrum_metadata_merged.csv" 
# 15 # threshold number of chr to count as whole-genome 
# "/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/aneuploidy_post/utils/test.pdf"

# get command line arguments 
args <- commandArgs(trailingOnly = TRUE)
# ploidy calls 
ploidy_calls <- args[1]
# parent to get gain/loss by 
parent <- args[2]
# metadata to remove day3 embryos 
metadata <- args[3]
# threshold number of affected chr for whole genome ploidy 
ploidy_threshold <- args[4]
# out plot for whole genome gain/loss counts 
out_ploidy_by_chr <- args[5]


# read in data 
ploidy_calls <- fread(ploidy_calls)
metadata <- fread(metadata)

# get functions `filter_data` and `day5_only`
setwd(".")
source("../../gwas/scripts/phenotypes/helper_functions/phenotyping_helper_functions.R")

# filter data (using functions from phenotyping and pasted below) 
ploidy_calls_filtered <- filter_data(ploidy_calls, parent, 
                                     bayes_factor_cutoff = 2, 
                                     nullisomy_threshold = 5, min_prob = 0.9) 
ploidy_calls_filtered <- day5_only(ploidy_calls_filtered, metadata)


# function to plot number of affected chr per embryo and return benchmarks 
affected_chr_per_embryo <- function(ploidy_calls_filtered, 
                                    max_meiotic = 3, 
                                    ploidy_threshold = 15) {
  
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
  num_meiotic <- length(which(affected_chr_per_embryo$num_affected_chromosomes 
                              <= max_meiotic))
  # count number of putative whole-genome affects 
  num_ploidy <- length(which(affected_chr_per_embryo$num_affected_chromosomes 
                             >= ploidy_threshold))
  
  # return output as list 
  output <- list(p1, num_meiotic, num_ploidy)
  return(output)
}

output <- affected_chr_per_embryo(ploidy_calls_filtered)
affected_chr_per_embryo_plot <- output[1]
num_meiotic <- output[2]
num_ploidy <- output[3]


# function to plot error type by chromosome 
plot_ploidy_by_chr <- function(ploidy_calls_filtered, 
                               as_proportion = TRUE, 
                               aneuploid_only = TRUE, 
                               remove_whole_genome = TRUE, 
                               ploidy_threshold = 15) {
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
    whole_genome <- affected_chr_per_embryo[affected_chr_per_embryo$num_affected_chromosomes < ploidy_threshold,]
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

ploidy_by_chr <- plot_ploidy_by_chr(ploidy_calls_filtered)

# write plot to pdf 
pdf(out_ploidy_by_chr)
ggsave(out_ploidy_by_chr, ploidy_by_chr)
