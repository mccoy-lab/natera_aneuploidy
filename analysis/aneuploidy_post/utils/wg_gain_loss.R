# Whole genome gain or loss by parent 

# load libraries
library(data.table)
library(dplyr)
library(ggplot2) 

# load data 
ploidy_calls <- fread("/data/rmccoy22/natera_spectrum/karyohmm_outputs/compiled_output/natera_embryos.karyohmm_v30a.bph_sph_trisomy.full_annotation.031624.tsv.gz")

bad_trios <- fread("/data/rmccoy22/natera_spectrum/karyohmm_outputs/compiled_output/invalid_trios.010324.txt", header = FALSE)

ploidy_calls_no_bad_trios <- ploidy_calls[!(ploidy_calls$child %in% bad_trios$V3),]

ploidy_calls_v18 <- fread("/data/rmccoy22/natera_spectrum/karyohmm_outputs/compiled_output/natera_embryos.karyohmm_v18.bph_sph_trisomy.full_annotation.112023.filter_bad_trios.tsv.gz")
ploidy_calls_v14 <- fread("/data/rmccoy22/natera_spectrum/karyohmm_outputs/compiled_output/natera_embryos_v2.karyohmm_v14.070923.tsv.gz")

# filter data as per phenotyping 
source("/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/scripts/phenotypes/helper_functions/phenotyping_helper_functions.R")


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
  
  output <- list(triploidies, haploidies)
  return(output)
}

maternal <- whole_genome(ploidy_calls_filtered, "mother")
paternal <- whole_genome(ploidy_calls_filtered, "father")

# get affected trios - paternal 
affected_trios_pat_trip <- ploidy_calls[ploidy_calls$child %in% paternal[1][[1]]$child,] %>%
  select(mother, father, child) %>% 
  distinct()
affected_trios_pat_hap <- ploidy_calls[ploidy_calls$child %in% paternal[2][[1]]$child,] %>%
  select(mother, father, child) %>% 
  distinct()

# affected paternal trios 
affected_paternal_trios <- rbind(affected_trios_pat_trip, affected_trios_pat_hap)
write.csv(affected_paternal_trios, "/scratch16/rmccoy22/scarios1/sandbox/paternal_wholegenome.csv", row.names = FALSE, quote = FALSE)


# get affected trios - maternal
affected_trios_mat_trip <- ploidy_calls[ploidy_calls$child %in% maternal[1][[1]]$child,] %>%
  select(mother, father, child) %>% 
  distinct()
affected_trios_mat_hap <- ploidy_calls[ploidy_calls$child %in% maternal[2][[1]]$child,] %>%
  select(mother, father, child) %>% 
  distinct()

# affected maternal trios 
affected_maternal_trios <- rbind(affected_trios_mat_trip, affected_trios_mat_hap)
write.csv(affected_maternal_trios, "/scratch16/rmccoy22/scarios1/sandbox/maternal_wholegenome.csv", row.names = FALSE, quote = FALSE)



