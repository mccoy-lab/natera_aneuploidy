## GWAS for phenotypes from Natera

# load libraries
library(data.table)
library(BEDMatrix)
library(dplyr)
library(pbmcapply)
library(purrr)

# Usage: /scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/scripts/gwas/gwas_all.R \
# "/data/rmccoy22/natera_spectrum/data/summary_metadata/spectrum_metadata_merged.csv"
# "/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/results/gwas/subsets/spectrum_imputed_chr21_rehead_filter_cpra_18.bed" \
# "/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/results/discover_validate_split_mother.txt" \
# "/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/results/parental_genotypes_pcs/parental_genotypes.eigenvec" \
# "/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/results/phenotypes/haploidy_by_mother.csv" \
# "/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/results/gwas/subsets/spectrum_imputed_chr21_rehead_filter_cpra_18.bim" \
# "discovery" \ # dataset_type
# "haploidy" \ # phenotype_name
# "mother" \ # parent
# 16 \ # number of threads
# "/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/results/gwas/gwas_haploidy_by_mother_discovery_9.tsv"

# get command line arguments
args <- commandArgs(trailingOnly = TRUE)
# metadata for all years
metadata <- args[1]
# parental genotypes as described
#https://github.com/mccoy-lab/natera_spectrum/tree/main/genotyping
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

# source Rscript with helper functions
#source("helper_functions/gwas_helper_functions.R")

## Functions for use in executing Natera GWAS
# /scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/scripts/gwas/gwas_all.R

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
  
  # Get caseIDs belonging to mothers in dataset type
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
                   phenotype, pcs) {
  
  # Make genotype file
  gt <- data.table(names(bed[bed_dataset_indices, snp_index]),
                   unname(bed[bed_dataset_indices, snp_index])) %>%
    setnames(., c("array", "alt_count")) %>%
    .[, array := sub("(.*?_.*?)_.*", "\\1", array)]
  
  gt <- merge(gt, metadata, by = "array") %>%
    merge(phenotype, by = "array") %>%
    merge(pcs, by = "array") %>%
    .[!duplicated(array)]
  
  # Get whether each mother is an egg donor
  gt$egg_donor_factor <- factor(gt$egg_donor, levels = c("", "yes"))
  
  return(gt)
}


# Function to make GWAS model based on parent and phenotype
make_model <- function(parent, phenotype_name) {
  
  # Choose variable columns and model family based on phenotype name
  if (grepl("ploidy", phenotype_name)) {
    response_variable <- "cbind(aneu_true, aneu_false)"
    family <- "quasibinomial"
  } else if (phenotype_name == "embryo_count") {
    response_variable <- "num_embryos"
    family <- "quasipoisson"
  } else if (phenotype_name == "maternal_age") {
    response_variable <- "weighted_age"
    family <- "gaussian"
  } else if (phenotype_name == "sex_ratio") {
    response_variable <- "cbind(XY, XX)"
    family <- "quasibinomial"
  } else {
    stop("Invalid 'phenotype_name' argument.")
  }
  
  # Set formula string to include age column (unless age phenotype) and 
  # response variable
  formula_string <- paste0(response_variable, " ~ PC1 + PC2 + PC3 + PC4 + ", 
                             "PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + ", 
                             "PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 +", 
                             " PC19 + PC20 + alt_count + egg_donor_factor", 
                           collapse = "")
  
  if (phenotype_name != "maternal_age") {
   formula_string <- paste0(formula_string, " + weighted_age", collapse = "")
  }
   
  # Return model for use in GWAS
  return(list(formula_string = formula_string, family = family))
}


# Calculate GWAS at each genomic site
gwas_per_site <- function(snp_index, bed, bim, pcs, phenotype,
                          bed_dataset_indices, metadata, parent,
                          phenotype_name) {
  
  # Get characteristics for each site
  snp_name <- colnames(bed)[snp_index]
  snp_chr <- bim[snp_index]$chr
  snp_pos <- bim[snp_index]$pos
  
  # Get genotype info for each site
  gt <- get_gt(bed, bed_dataset_indices, snp_index, metadata, phenotype, pcs)
  
  # Make GWAS model
  model <- make_model(parent, phenotype_name)
  formula_string <- model$formula_string
  family <- model$family
  m1 <- glm(formula_string, family = family, data = gt) %>%
    summary()
  
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

  # Return GWAS for a given site
  return(output)
}


# Run GWAS across dataset
run_gwas <- function(dataset_type, discovery_test, metadata, bed, bim, pcs, 
                     phenotype, parent, phenotype_name, threads = 32) {
  
  # Subset bed file corresponding to correct dataset type
  bed_dataset <- discovery_test_split(dataset_type, discovery_test, 
                                      metadata, bed)
  
  # Get indices to execute function
  bed_dataset_indices <- 1:nrow(bed_dataset)
  
  # Calculate GWAS across each site
  gwas_results <- pbmclapply(1:ncol(bed_dataset),
                             function(x) gwas_per_site(x, bed_dataset, bim,
                                                       pcs, phenotype,
                                                       bed_dataset_indices,
                                                       metadata, parent,
                                                       phenotype_name),
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
                            pcs, phenotype, parent, phenotype_name, threads)

# write to file
write.table(gwas_results_dt, out_fname, append = FALSE, sep = "\t", dec = ".",
            row.names = FALSE, col.names = TRUE, quote = FALSE)
