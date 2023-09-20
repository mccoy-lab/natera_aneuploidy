## GWAS for maternal meiotic errors + maternal genotypes 

# load libraries
library(data.table)
library(BEDMatrix)
library(dplyr)
library(pbmcapply)
library(purrr)

# take in arguments 
args <- commandArgs(trailingOnly = TRUE)
metadata <- args[1]
bed <- args[2]
discovery_test_f <- args[3]
eigenvec <- args[4]
pheno <- args[5]
bim <- args[6]
out_fname <- args[7]
dataset_type <- args[8] # discovery or test

# read in files from arguments 
metadata <- fread(metadata)
bed <- BEDMatrix(bed)
discovery_test_f <- fread(discovery_test_f)
eigenvec <- fread(eigenvec)
pheno <- fread(pheno, sep = ",")
bim <- fread(bim) %>% setnames(., c("chr", "snp_id", "drop", "pos", "ref", "alt"))


# separate into discovery or test 
discovery_test <- function(dataset_type, metadata, bed, discovery_test_f) {
  if (dataset_type == "discovery") {
    dataset <- discovery_test_f[discovery_test_f$is_discovery == TRUE,]
  } else if (data_type == "test") {
    dataset <- discovery_test_f[discovery_test_f$is_discovery == FALSE,]
  }
  # get caseIDs belonging to mothers in preferred set
  metadata_set <- metadata[metadata$array %in% dataset$array,]
  # add column for array_array to match rownames of bed
  metadata_set$array_array <- paste0(metadata_set$array, "_", metadata_set$array)
  # subset bed file to just the mothers in the relevant set 
  bed_indices <- which(rownames(bed) %in% metadata_set$array_array)
  bed_dataset <- bed[bed_indices,]
  # return relevant subset of bed file (discovery vs. test)
  return(bed_dataset)
}

# get relevant parts of bed file 
bed_dataset <- discovery_test(dataset_type, metadata, bed, discovery_test_f)
# get variables to call function 
bed_dataset_indices <- 1:nrow(bed_dataset)

# format pca 
pca_scores <- as.data.table(eigenvec)
# keep PCs 1 through 20
pca_scores <- pca_scores[,2:22]
colnames(pca_scores) <- c("array", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11", "PC12", "PC13", "PC14", "PC15", "PC16", "PC17", "PC18", "PC19", "PC20")

# function to run GWAS on each site
gwas <- function(snp_index, genotypes, phenotypes, metadata, locs, pcs, subject_indices) {
  
  gt <- data.table(names(genotypes[subject_indices, snp_index]), unname(genotypes[subject_indices, snp_index])) %>%
    setnames(., c("array", "alt_count")) %>%
    .[, array := sub("(.*?_.*?)_.*", "\\1", array)]
  
  snp_name <- colnames(genotypes)[snp_index]
  snp_chr <- locs[snp_index]$chr
  snp_pos <- locs[snp_index]$pos
  
  gt <- merge(gt, metadata, by = "array") %>%
    merge(pheno, by = "casefile_id") %>%
	  merge(pcs, by = "array") %>%
	  .[!duplicated(array)]
  
  m1 <- glm(cbind(num_aneuploid, num_euploid) ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + patient_age + alt_count, family = "quasibinomial", data = gt) %>% summary()
  
  coef <- data.table(term = rownames(m1$coefficients), m1$coefficients)

  return(data.table(snp = snp_name, pos = snp_pos, 
                    beta = unlist(coef[term == "alt_count", 2]),
                    se = unlist(coef[term == "alt_count", 3]),
                    t = unlist(coef[term == "alt_count", 4]),
                    p.value = unlist(coef[term == "alt_count", 5])))
  
}

# run GWAS
gwas_results_discovery <- pbmclapply(1:ncol(bed_dataset), function(x) gwas(x, bed_dataset, pheno, metadata, bim, pca_scores, bed_dataset_indices), mc.cores = 48L)

# bind output
gwas_results_dt <- rbindlist(gwas_results_discovery[unlist(map(gwas_results_discovery, is.data.table))]) %>%
  .[!is.na(p.value)] %>%
  .[, c("snp_id", "effect_allele") := tstrsplit(snp, "_", fixed = TRUE)] %>%
  merge(., bim[, -"snp_id"], by = "pos") %>%
  setorder(., p.value)

# write to file
write.table(gwas_results_dt, out_fname, append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE, quote = FALSE)

