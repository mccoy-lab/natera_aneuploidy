# Count and identify genome gains or losss by parent 

# load libraries
library(data.table)
library(dplyr)
library(ggplot2) 

# Usage: /aneuploidy_post/utils/wg_gain_loss.R \ 
# "/data/rmccoy22/natera_spectrum/karyohmm_outputs/compiled_output/natera_embryos.karyohmm_v30a.bph_sph_trisomy.full_annotation.031624.tsv.gz" \
# "mother"
# "/data/rmccoy22/natera_spectrum/data/summary_metadata/spectrum_metadata_merged.csv" 
# 15 # threshold number of chr to count as whole-genome 
# "/utils/results/whole_genome_mother.txt"
# "/utils/results/affected_trios_mother.txt"

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
# outfile for whole genome gain/loss counts 
out_whole_genome <- args[5]
# outfile for affected trios 
out_affected_trios <- args[6]

# read in data 
ploidy_calls <- fread(ploidy_calls)
metadata <- fread(metadata)

# get functions `filter_data` and `day5_only`
source("/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/scripts/phenotypes/helper_functions/phenotyping_helper_functions.R")

# filter ploidy calls to keep only high quality data and day5 embryos 
ploidy_calls_filtered <- filter_data(ploidy_calls, parent, 
                                     bayes_factor_cutoff = 2,
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
  
  output <- list(triploidies, haploidies)
  return(output)
}

whole_genome_gain_loss <- whole_genome(ploidy_calls_filtered, parent, 
                                       ploidy_threshold)

write.csv(whole_genome_gain_loss, out_whole_genome, 
          row.names = FALSE, quote = FALSE)



# find affected trios 
affected_trios <- function(ploidy_calls_filtered, whole_genome_gain_loss, 
                           parent) {
  triploidies <- ploidy_calls_filtered[ploidy_calls_filtered$child %in% 
                                         whole_genome_gain_loss[1][[1]]$child,] %>%
    select(mother, father, child) %>% 
    distinct() %>% 
    mutate(ploidy = "triploidy") %>% 
    mutate(parent = parent)
  
  haploidies <- ploidy_calls_filtered[ploidy_calls_filtered$child %in% 
                                         whole_genome_gain_loss[2][[1]]$child,] %>%
    select(mother, father, child) %>% 
    distinct() %>% 
    mutate(ploidy = "haploidy") %>% 
    mutate(parent = parent)
  
  ploidies <- rbind(triploidies, haploidies)
  
}

# compute affected trios and write to file 
affected_trios_output <- affected_trios(ploidy_calls_filtered, whole_genome_gain_loss, parent)
write.csv(affected_trios_output, out_affected_trios, row.names = FALSE, quote = FALSE)

