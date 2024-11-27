## For top hits from regular model, run GWAS again using a linear mixed model
## GWAS for phenotypes from Natera

# =================
# author: Sara A. Carioscia, Biology Dept., Johns Hopkins University
# email: scarios1@jhu.edu
# 2024
# =================

# load libraries
library(data.table)
library(BEDMatrix)
library(dplyr)
library(pbmcapply)
library(purrr)
# library for glmer
library(lme4)
# library for ztpoisson 
library(countreg)

# Usage: /scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/scripts/gwas/glmmer.R \
# "/data/rmccoy22/natera_spectrum/data/summary_metadata/spectrum_metadata_merged.csv"
# "/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/results/gwas/subsets/spectrum_imputed_chr21_rehead_filter_cpra_18.bed" \
# "/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/results/gwas/intermediate_files/discover_validate_split_mother.txt" \
# "/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/results/parental_genotypes_pcs/parental_genotypes.eigenvec" \
# "/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/results/phenotypes/chr22_aneuploidy_by_mother.csv" \
# "/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/results/gwas/subsets/spectrum_imputed_chr21_rehead_filter_cpra_18.bim" \
# "discovery" \ # dataset_type
# "haploidy" \ # phenotype_name
# "mother" \ # parent
# 16 \ # number of threads
# "/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/results/gwas/gwas_lmm_haploidy_by_mother_discovery_9.tsv"

# get command line arguments
args <- commandArgs(trailingOnly = TRUE)
# metadata for all years
metadata <- args[1]
# parental genotypes as described
# https://github.com/mccoy-lab/natera_spectrum/tree/main/genotyping
bed <- args[2]
# discovery vs. test split for applicable parent
discovery_test <- args[3]
# principal components for parents
pcs <- args[4]
# phenotype file of interest
phenotype <- args[5]
# file containing snp name and map positions, corresponds with bed file
bim <- args[6]
# string variable noting whether to use the discovery or test set for parents
dataset_type <- args[7]
# string variable naming the phenotype being considered (e.g., maternal_meiotic)
phenotype_name <- args[8]
# string variable noting which parent to GWAS
parent <- args[9]
# number of threads to use in mc.cores for GWAS
threads <- as.numeric(args[10])
# output file name
out_fname <- args[11]


# Subset data to include only discovery or test
discovery_test_split <- function(dataset_type, discovery_test, metadata, bed) {
  
  # Check dataset type and subset accordingly
  if (dataset_type == "discovery") {
    dataset <- discovery_test[discovery_test$is_discovery == TRUE,]
  } else if (dataset_type == "test") {
    dataset <- discovery_test[discovery_test$is_discovery == FALSE,]
  } else {
    stop("Invalid 'dataset_type' argument.")
  }
  
  # Get arrays belonging to mothers in dataset type
  metadata_set <- metadata[metadata$array %in% dataset$array,]
  # add column for array_array to match rownames of bed
  metadata_set$array_array <- paste0(metadata_set$array, "_",
                                     metadata_set$array)
  
  # Subset bed file to just the mothers in the relevant set 
  bed_indices <- which(rownames(bed) %in% metadata_set$array_array)
  bed_dataset <- bed[bed_indices,]
  
  # Return bed file subset by dataset type
  return(bed_dataset)
}


# Function to pre-process genotype info for each site
get_gt <- function(bed, bed_dataset_indices, snp_index, metadata, 
                   phenotype, pcs, parent, dataset_type) {
  
  # Make genotype file
  gt <- data.table(names(bed[bed_dataset_indices, snp_index]),
                   unname(bed[bed_dataset_indices, snp_index])) %>%
    setnames(., c("array", "alt_count")) %>%
    # Apply different rules based on the presence of "CYTO12AB" in the array
    .[, array := ifelse(grepl("CYTO12AB", array),
                        sub("(.*?_.*?_.*?)_.*", "\\1", array),  # Regex for CYTO12AB cases
                        sub("(.*?_.*?)_.*", "\\1", array))]    # Regex for other cases
  
  # Set array column to be the parent of interest 
  if (parent == "mother") {
    phenotype[, array := mother]
  } else if (parent == "father") {
    phenotype[, array := father]
  } else {
    stop("Invalid parent value. It should be either 'mother' or 'father'.")
  }
  
  # Keep only cycles from mothers that are in the discovery set 
  if (dataset_type == "discovery") {
    phenotype <- phenotype[phenotype$array %in% 
                                       discovery_test[is_discovery == TRUE]$array,]
  } else {
    phenotype <- phenotype[phenotype$array %in% 
                             discovery_test[is_discovery == FALSE]$array,]
  }
  
  # Merge pcs for each phenotype 
  merged_with_pcs <- merge(phenotype, pcs, by = "array", all.x = TRUE)
  
  # Merge in the gt info 
  gt <- merge(merged_with_pcs, gt, by = "array", all.x = TRUE)
  
  # Assign egg and sperm donor ages 
  # get average egg donor age from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7530253/
  # get average sperm donor age from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9118971/#:~:text=Donors%20were%20aged%2027%20years,aged%2030%20years%20and%20younger.
  gt[egg_donor == "yes", patient_age_cycle := as.numeric(25)]
  gt[sperm_donor == "yes", partner_age_cycle := as.numeric(27)]
  
  return(gt)
}


# Function to make GWAS model based on phenotype
make_model <- function(phenotype_name) {
  
  # Choose variable columns and model family based on phenotype name
  if (grepl("ploidy", phenotype_name)) {
    response_variable <- "cbind(aneu_true, aneu_false)"
    family <- "binomial"
  } else if (phenotype_name == "embryo_count") {
    response_variable <- "num_embryos"
    family <- "ztpoisson"
  } else if (phenotype_name == "maternal_age") {
    response_variable <- "patient_age_cycle"
    family <- "gaussian"
  } else if (phenotype_name == "sex_ratio") {
    response_variable <- "cbind(XY, XX)"
    family <- "binomial"
  } else {
    stop("Invalid 'phenotype_name' argument.")
  }
  
  # Set formula string to include both parental ages, both donor statuses, PCs, 
  # and response variable (omit patient_age as a covariate when maternal 
  # age is the phenotype)
  formula_string <- paste0(response_variable, " ~ (1 | array) + PC1 + PC2 + PC3 + ",
                           "PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + ", 
                           "PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + ",
                           "PC19 + PC20 + (egg_donor == 'yes') + ",
                           "(sperm_donor == 'yes') + ", 
                           "scale(partner_age_cycle) + alt_count", 
                           collapse = "")
  
  # For all phenotypes other than maternal age, include patient_age as covariate
  if (phenotype_name != "maternal_age") {
    formula_string <- paste0(formula_string, " + scale(patient_age_cycle)", 
                             collapse = "")
  }
  
  # If the phenotype includes "age_interaction", create another formula in 
  # which there's an interaction effect between age and alt_count 
  if (grepl("age_interaction", phenotype_name)) {
    formula_string2 <- 
      gsub("alt_count \\+ scale\\(patient_age_cycle\\)", 
           "alt_count * scale(patient_age_cycle)", formula_string)
  } else {
    # if not age interaction just keep formula_string2 as is
    formula_string2 <- formula_string
  }
  
  # Return model for use in GWAS, including age interaction formula if necessary
  return(list(formula_string = formula_string, family = family,
              formula_string2 = formula_string2))
}


# Calculate GWAS at each genomic site
gwas_per_site <- function(snp_index, bed, bim, pcs, phenotype,
                          bed_dataset_indices, metadata, phenotype_name, 
                          parent, model, dataset_type) {
  
  # Get characteristics for each site
  snp_name <- colnames(bed)[snp_index]
  snp_pos <- bim[snp_index]$pos
  
  # Get genotype info for each site
  gt <- get_gt(bed, bed_dataset_indices, snp_index, metadata, phenotype, pcs,
               parent, dataset_type)
  
  # Make GWAS model
  formula_string <- model$formula_string
  family <- model$family
  model1 <- glmer(formula_string, family = family, nAGQ = 0, 
                  control = glmerControl(optimizer = "nloptwrap"),
                  data = gt)
  m1 <- model1 %>% summary()
  
  # Get info for GWAS output
  coef <- data.table(term = rownames(m1$coefficients), m1$coefficients)
  
  # Calculate alternate allele frequency
  gt_filtered  <- gt %>% filter(!is.na(alt_count))
  alt_af <- sum(gt_filtered$alt_count) / (2 * nrow(gt_filtered))
  
  # Produce GWAS output for each site
  output <- data.table(snp = snp_name,
                       pos = snp_pos,
                       beta = unlist(coef[term == "alt_count", 2]),
                       se = unlist(coef[term == "alt_count", 3]),
                       t = unlist(coef[term == "alt_count", 4]),
                       p.value = unlist(coef[term == "alt_count", 5]),
                       af = alt_af)
  
  # For age interaction, do the GWAS for the age interaction model as well, 
  # then do the ANOVA to compare the models and bind the output 
  if (grepl("age_interaction", phenotype_name)) {
    # Make GWAS model for age interaction 
    formula_string2 <- model$formula_string2
    model2 <- glmer(formula_string2, family = family, nAGQ = 0, 
                    control = glmerControl(optimizer = "nloptwrap"),
                    data = gt) 
    m2 <- model2 %>% 
      summary()
    # Get p-value for age-altcount interaction 
    coef2 <- data.table(term = rownames(m2$coefficients), m2$coefficients)
    
    # Run ANOVA to compare the two models 
    model_comparison <- anova(model1, model2)
    
    # Output the data 
    output <- data.table(snp = snp_name,
                         pos = snp_pos,
                         p.value_age_interaction = unlist(coef2[term == "alt_count:scale(patient_age_cycle)", 5]),
                         p.value_comparison = model_comparison[2,8])
  
  # Return GWAS for a given site
  return(output)
}


# Define a function to process a chunk of SNPs
gwas_per_chunk <- function(snp_indices, bed, bim, pcs, phenotype,
                           bed_dataset_indices, metadata, phenotype_name, 
                           parent, model, dataset_type) {
  chunk_results <- lapply(snp_indices, function(snp_index) {
    gwas_per_site(snp_index, bed, bim, pcs, phenotype, 
                  bed_dataset_indices, metadata, 
                  phenotype_name, parent, model, dataset_type)
  })
  
  # Combine the results for all SNPs in this chunk
  return(rbindlist(chunk_results))
}


# Run GWAS across dataset
run_gwas <- function(dataset_type, discovery_test, metadata, bed, bim, pcs, 
                     phenotype, phenotype_name, parent, 
                     threads = 32) {
  
  # Subset bed file corresponding to correct dataset type
  bed_dataset <- discovery_test_split(dataset_type, discovery_test, 
                                      metadata, bed)
  
  # Get indices to execute function
  bed_dataset_indices <- 1:nrow(bed_dataset)
  
  # Make GWAS model
  model <- make_model(phenotype_name)
  
  # Calculate GWAS across each chunk of sites
  snp_indices <- 1:ncol(bed_dataset)
  snp_chunks <- split(snp_indices, ceiling(seq_along(snp_indices) / 100))
  
  # Run GWAS on each chunk 
  gwas_results <- pbmclapply(snp_chunks,
                             function(chunk) gwas_per_chunk(chunk, bed_dataset, 
                                                            bim, pcs, 
                                                            phenotype,
                                                            bed_dataset_indices,
                                                            metadata, 
                                                            phenotype_name, parent, 
                                                            model, dataset_type),
                             mc.cores = threads)
  # Bind output across all sites
  gwas_results_dt <-
    rbindlist(gwas_results[unlist(map(gwas_results, is.data.table))]) %>%
    .[!is.na(p.value)] %>%
    .[, c("snp_id", "effect_allele") := tstrsplit(snp, "_", fixed = TRUE)] %>%
    merge(., bim[, -"snp_id"], by = "pos") %>%
    setorder(., p.value)
  
  # Return GWAS output across all sites
  return(gwas_results_dt)
}


# read in files
metadata <- fread(metadata)
bed <- BEDMatrix(bed)
discovery_test <- fread(discovery_test)
pcs <- fread(pcs)
colnames(pcs)[1] <- "array"
phenotype <- fread(phenotype)
bim <- fread(bim) %>%
  setnames(., c("chr", "snp_id", "drop", "pos", "ref", "alt"))


# conduct GWAS across all sites
gwas_results_dt <- run_gwas(dataset_type, discovery_test, metadata, bed, bim,
                            pcs, phenotype, phenotype_name, parent,
                            threads)

# write to file
write.table(gwas_results_dt, out_fname, append = FALSE, sep = "\t", dec = ".",
            row.names = FALSE, col.names = FALSE, quote = FALSE)