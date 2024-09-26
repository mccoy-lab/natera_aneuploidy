## Make phenotype file for embryos affected by aneuploidy phenotypes

# load libraries
library(data.table)
library(tidyr)
library(dplyr)

# Usage:
# /scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/scripts/phenotypes/generate_phenotype_files.R \
# /data/rmccoy22/natera_spectrum/karyohmm_outputs/compiled_output/natera_embryos.karyohmm_v30a.bph_sph_trisomy.full_annotation.031624.tsv.gz \
# /scratch16/rmccoy22/abiddan1/natera_segmental/analysis/segmental_qc/results/tables/segmental_calls_postqc_refined.tsv.gz \ # segmental aneuploidy calls to remove chrom 
# mother \ # parent to group phenotype by
# /data/rmccoy22/natera_spectrum/data/summary_metadata/spectrum_metadata_merged.csv \
# maternal_meiotic_aneuploidy \ # phenotype name
# TRUE \ whether to keep only day 5 embryos 
# 2 \ keep only chroms with bayes factor greater than the threshold for bayes factor qc
# 5 \ remove embryos that had more chr with cn = 0 for than nullisomy_threshold
# 0.9 \ minimum posterior probability for each cn call
# 5 \ max number of affected chr to count for maternal meiotic phenotype
# 15 \ min number of affected chr to count for whole genome gain/loss
# TRUE \ whether to filter out putative mosaics 
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
# phenotype name
phenotype_name <- args[5]
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
# whether to filter out mosaic embryos
filter_mosaics <- args[12]
# output file name
out_fname <- args[13]

# assign chromosome variable if phenotype is single-chromosome aneuploidy 
if (grepl("^chr[0-9]+_aneuploidy$", phenotype_name)) {
  # extract the chr number 
  chromosome <- as.integer(gsub("[^0-9]", "", phenotype_name))
} else {
  chromosome <- NULL
}


# Function to filter ploidy calls by quality for aneuploidy phenotypes 
filter_data <- function(ploidy_calls, parent, segmental_calls, 
                        bayes_factor_cutoff = 2, filter_day_5 = TRUE, 
                        nullisomy_threshold = 5, 
                        min_prob = 0.9, filter_mosaics = TRUE, 
                        mosaic_threshold = 0.340) {
  
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
  
  # Keep only rows with bayes factor greater than the threshold for 
  # bayes factor qc
  ploidy_calls <- ploidy_calls[ploidy_calls$bf_max > bayes_factor_cutoff, ]
  
  # Keep only chromosomes with sufficiently high probability cn call
  ploidy_calls <- ploidy_calls[ploidy_calls$post_max > min_prob, ]
  
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
  
  # Remove trisomies that are likely mosaic 
  if (filter_mosaics == TRUE) {
    ploidy_calls <- ploidy_calls[!(ploidy_calls$bf_max_cat %in% c("3m", "3p") & 
         (ploidy_calls$post_bph_centro < mosaic_threshold & 
            ploidy_calls$post_bph_noncentro < mosaic_threshold)), ]
  }
  
  return(ploidy_calls)
}

# create phenotypes for traits, both aneuploidy and metadata-based 
# (maternal meiotic aneuploidy, triploidy, haploidy, embryo count, maternal age)  
make_phenotype <- function(metadata, parent, phenotype_name, ploidy_calls, 
                           segmental_calls, bayes_factor_cutoff = 2, 
                           filter_day_5 = TRUE, nullisomy_threshold = 5, 
                           max_meiotic = 5, min_ploidy = 15, 
                           filter_mosaics = TRUE, mosaic_threshold = 0.340, 
                           chromosome = NULL) {
  
  # assign all members of each family to the same mother, across casefile IDs 
  metadata_mothers <- metadata %>%
    mutate(mother_id = ifelse(family_position == "mother", array, NA)) %>%
    fill(mother_id, .direction = "downup")
  
  # if values for embryo or father ages were NA, fill age from mother 
  # associated with that casefileID
  metadata_mothers <- metadata_mothers %>%
    group_by(casefile_id) %>%
    mutate(
      patient_age = ifelse(family_position %in% c("child", "father") & 
                             is.na(patient_age) & egg_donor != "yes", 
                           patient_age[family_position == "mother"], 
                           patient_age),
      partner_age = ifelse(family_position %in% c("child", "father") & 
                             is.na(partner_age) & sperm_donor != "yes", 
                           partner_age[family_position == "mother"], 
                           partner_age)
    ) %>%
    ungroup()
  
  # assign visit and keep only day 5 embryos
  child_data <- metadata_mothers %>%
    filter(family_position == "child", sample_scale == "few_cells") %>%
    arrange(mother_id, patient_age) %>%
    group_by(mother_id) %>%
    mutate(visit_id = dense_rank(patient_age)) %>%
    ungroup() %>%
    distinct(array, .keep_all = TRUE)
  
  # count number of visits per mother 
  num_visits <- child_data %>%
    group_by(mother_id) %>%
    summarise(num_visits = n_distinct(visit_id)) %>%
    ungroup()
  
  # track parental age at each visit 
  age_produced <- child_data %>%
    group_by(mother_id, visit_id) %>%
    summarise(patient_age_cycle = mean(patient_age, na.rm = TRUE),
              partner_age_cycle = mean(partner_age, na.rm = TRUE), 
              egg_donor = first(egg_donor),
              sperm_donor = first(sperm_donor),
              year = first(year)) %>%
    ungroup()
  
  # for aneuploidy phenotypes, count aneuploid/euploid embryos per visit
  if (grepl("ploidy", phenotype_name)) {
    # select ploidy status based on phenotype and parent 
    if (phenotype_name == "triploidy" & parent == "mother") {
      cn <- "3m"
    } else if (phenotype_name == "triploidy" & parent == "father") {
      cn <- "3p"
    } else if (phenotype_name == "haploidy" & parent == "mother") {
      cn <- "1p"
    } else if (phenotype_name == "haploidy" & parent == "father") {
      cn <- "1m"
    } else if (grepl("maternal_meiotic_aneuploidy", phenotype_name) | 
               grepl("^chr[0-9]+_aneuploidy$", phenotype_name)) {
      cn <- c("3m", "1p")
    } else if (phenotype_name == "complex_aneuploidy") {
      cn <- c("0", "1m", "1p", "3m", "3p")
    } 
    
    # filter karyohmm data (quality control, day 5, posterior probabilities)
    ploidy_calls <- filter_data(ploidy_calls, parent, segmental_calls, 
                                bayes_factor_cutoff, filter_day_5,
                                nullisomy_threshold, min_prob, filter_mosaics,
                                mosaic_threshold)
    
    # create a lookup table for array and visit_id to match each child's 
    # array with its ploidy calls 
    visit_lookup <- child_data %>%
      dplyr::select(array, visit_id)
    
    # add visit number to each embryo
    ploidy_calls <- ploidy_calls %>%
      left_join(visit_lookup, by = c("child" = "array"))
    
    # if phenotype is single-chromosome, keep only ploidy calls from that chr
    if (grepl("^chr[0-9]+_aneuploidy$", phenotype_name)) {
      ploidy_calls <- ploidy_calls[as.integer(sub("chr", "", 
                                          ploidy_calls$chrom)) == chromosome, ]
    }
    
    # count number of aneuploid and euploid embryos in each visit 
    result <- ploidy_calls %>%
      group_by(mother, father, child, visit_id) %>%
      summarise(num_affected = sum(bf_max_cat %in% cn), 
                unique_bf_max_cat = 
                  n_distinct(bf_max_cat[bf_max_cat %in% cn])) %>%
      mutate(
        is_ploidy = case_when(
          grepl("maternal_meiotic_aneuploidy", phenotype_name) ~ 
            ifelse(num_affected > 0 & num_affected <= max_meiotic, 
                   "aneu_true", "aneu_false"),
          grepl("^chr[0-9]+_aneuploidy$", phenotype_name) ~ 
            ifelse(num_affected == 1, 
                   "aneu_true", "aneu_false"),
          phenotype_name == "complex_aneuploidy" ~ 
            ifelse(num_affected > 0 & unique_bf_max_cat >= 2, 
                   "aneu_true", "aneu_false"),
          TRUE ~ 
            ifelse(num_affected > min_ploidy, "aneu_true", "aneu_false")
        )
      ) %>%
      count(visit_id, is_ploidy) %>%
      pivot_wider(names_from = is_ploidy, values_from = n, 
                  values_fill = list(n = 0))
    
    # group by family 
    mother_summary <- result %>%
      group_by(mother, father, visit_id) %>%
      summarise(
        aneu_true = sum(aneu_true),     
        aneu_false = sum(aneu_false),   
        total_embryos = sum(aneu_true + aneu_false)  
      ) %>%
      left_join(num_visits, by = c("mother" = "mother_id")) %>%
      left_join(age_produced, 
                by = c("mother" = "mother_id", "visit_id" = "visit_id")) %>%
      ungroup()
    
  } else {
    # calculate phenotypes based on metadata (embryo count, maternal age) 
    mother_summary <- child_data %>%
      group_by(mother_id, visit_id) %>%
      summarise(
        num_embryos = n() # Count number of children per visit
      ) %>%
      left_join(num_visits, by = c("mother_id" = "mother_id")) %>%
      left_join(age_produced, 
                by = c("mother_id" = "mother_id", "visit_id" = "visit_id")) %>%
      ungroup()
    
  }
  
  # change name of first column to mother
  colnames(mother_summary)[1] <- "mother"
  
  return(mother_summary)
}

# Read in embryos and metadata
ploidy_calls <- fread(ploidy_calls)
segmental_calls <- fread(segmental_calls)
metadata <- fread(metadata)

# Run phenotype 
pheno_by_parent <- make_phenotype(metadata, parent, phenotype_name, ploidy_calls, 
                                  segmental_calls, bayes_factor_cutoff, 
                                  filter_day_5, nullisomy_threshold, 
                                  max_meiotic, min_ploidy, filter_mosaics, 
                                  mosaic_threshold = 0.340, chromosome)

# Write phenotype info to file
write.csv(pheno_by_parent, out_fname, row.names = FALSE)
