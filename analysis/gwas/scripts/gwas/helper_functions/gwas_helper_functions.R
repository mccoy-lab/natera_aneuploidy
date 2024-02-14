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
get_gt <- function(pcs, phenotype, bim, bed, bed_dataset_indices, snp_index, metadata) {
  
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


# Function to make GWAS model based on phenotype and parent 
make_model <- function(parent, phenotype_name, gt) {
  
  # Determine the age covariate column based on the "parent" argument
  if (parent == "mother") {
    age_column <- "patient_age"
  } else if (parent == "father") {
    age_column <- "partner_age"
  } else {
    stop("Invalid 'parent' argument.")
  }
  
  # Choose phenotype columns and model family based on phenotype 
  if (grepl("ploidy", phenotype_name)) {
    response_variable <- "cbind(aneu_true, aneu_false)"
    family <- "quasibinomial"
  } else if (phenotype_name == "embryo_count") {
    response_variable <- "num_embryos"
    family <- "quasipoisson"
  } else if (phenotype_name == "maternal_age") {
    response_variable <- "weighted_age"
    family <- "gaussian"
  } else {
    stop("Invalid 'phenotype_name' argument.")
  }
  
  
  
  # Set formula string to include age column and response variable
  formula_string <- paste0(response_variable, " ~ PC1 + PC2 + ",
                           "PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + ", 
                           age_column, " + alt_count + egg_donor_factor", collapse = "")
  
  # Create a model using the glm 
  m1 <- glm(formula_string, family = family, data = gt) %>% 
    summary()
  
  # Return model for use in GWAS 
  return(m1)
  
}


# Calculate GWAS at each genomic site
gwas_per_site <- function(snp_index, bed, bim, pcs, phenotype, bed_dataset_indices, metadata, parent, phenotype_name) {
  
  # Get characteristics for each site
  snp_name <- colnames(bed)[snp_index]
  snp_chr <- bim[snp_index]$chr
  snp_pos <- bim[snp_index]$pos
  
  # Get genotype info for each site 
  gt <- get_gt(pcs, phenotype, bim, bed, bed_dataset_indices, snp_index, metadata)
  
  # Make GWAS model 
  m1 <- make_model(parent, phenotype_name, gt)
  
  # Get info for GWAS output 
  coef <- data.table(term = rownames(m1$coefficients), m1$coefficients)
  # Calculate alternate allele frequency 
  gt_filtered  <- gt %>% filter(!is.na(alt_count))
  alt_af <- sum(gt_filtered$alt_count) / (2*nrow(gt_filtered))
  
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
run_gwas <- function(dataset_type, discovery_test, metadata, bed, bim, pcs, phenotype, parent, phenotype_name, threads = 32) {
  
  # Subset bed file corresponding to correct dataset type 
  bed_dataset <- discovery_test_split(dataset_type, discovery_test, metadata, bed)
  
  # Get indices to execute function
  bed_dataset_indices <- 1:nrow(bed_dataset)
  
  # Calculate GWAS across each site
  gwas_results <- pbmclapply(1:ncol(bed_dataset),
                             function(x) gwas_per_site(x, bed_dataset, bim, pcs, phenotype, bed_dataset_indices, metadata, parent, phenotype_name),
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




