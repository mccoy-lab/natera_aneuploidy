## GWAS for maternal meiotic errors + parental genotypes 

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
discovery_test <- args[3] # either maternal or paternal info read in 
eigenvec <- args[4]
pheno <- args[5]
bim <- args[6]
dataset_type <- args[7] # discovery or test
phenotype <- args[8] # aneuploidy type, embryo count, or age 
out_fname <- args[9]

# read in files from arguments 
metadata <- fread(metadata)
bed <- BEDMatrix(bed)
discovery_test <- fread(discovery_test)
eigenvec <- fread(eigenvec)
pheno <- fread(pheno, sep = ",")
bim <- fread(bim) %>% setnames(., c("chr", "snp_id", "drop", "pos", "ref", "alt"))

print("read in all files")

# separate into discovery or test 
discovery_test_split <- function(dataset_type, metadata, bed, discovery_test) {
  if (dataset_type == "discovery") {
    dataset <- discovery_test[discovery_test$is_discovery == TRUE,]
  } else if (data_type == "test") {
    dataset <- discovery_test[discovery_test$is_discovery == FALSE,]
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
bed_dataset <- discovery_test_split(dataset_type, metadata, bed, discovery_test)

print("split discovery/test")
# get variables to call function 
bed_dataset_indices <- 1:nrow(bed_dataset)

print("got bed indices")

# format pca 
pca_scores <- as.data.table(eigenvec)
# keep PCs 1 through 20
pca_scores <- pca_scores[,2:12]
colnames(pca_scores) <- c("array", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")
print("formatted pca")

# function to run GWAS on each site if it's a ploidy phenotype (triploid, haploid, maternal meiotic aneuploidy)
gwas_aneuploidy <- function(snp_index, genotypes, phenotypes, metadata, locs, pcs, subject_indices) {
  
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
  
  m1 <- glm(cbind(num_affected, num_euploid) ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + patient_age + alt_count, family = "quasibinomial", data = gt) %>% summary()
  
  coef <- data.table(term = rownames(m1$coefficients), m1$coefficients)

  return(data.table(snp = snp_name, pos = snp_pos, 
                    beta = unlist(coef[term == "alt_count", 2]),
                    se = unlist(coef[term == "alt_count", 3]),
                    t = unlist(coef[term == "alt_count", 4]),
                    p.value = unlist(coef[term == "alt_count", 5]))) 
}

# function to run GWAS on each site if it's embryo count 
gwas_embryo_count <- function(snp_index, genotypes, phenotypes, metadata, locs, pcs, subject_indices) {
  
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
  
  m1 <- glm(num_embryos ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + weighted_age + alt_count, family = "poisson", data = gt) %>% summary()
  
  coef <- data.table(term = rownames(m1$coefficients), m1$coefficients)

  return(data.table(snp = snp_name, pos = snp_pos, 
                    beta = unlist(coef[term == "alt_count", 2]),
                    se = unlist(coef[term == "alt_count", 3]),
                    t = unlist(coef[term == "alt_count", 4]),
                    p.value = unlist(coef[term == "alt_count", 5])))
}

# function to run GWAS on each site for maternal age 
gwas_maternal_age <- function(snp_index, genotypes, phenotypes, metadata, locs, pcs, subject_indices) {
  
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
  
  m1 <- glm(weighted_age ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + alt_count, family = "poisson", data = gt) %>% summary()
  
  coef <- data.table(term = rownames(m1$coefficients), m1$coefficients)

  return(data.table(snp = snp_name, pos = snp_pos, 
                    beta = unlist(coef[term == "alt_count", 2]),
                    se = unlist(coef[term == "alt_count", 3]),
                    t = unlist(coef[term == "alt_count", 4]),
                    p.value = unlist(coef[term == "alt_count", 5])))  
}

# run GWAS based on phenotype passed argument 
if (phenotype == "maternal_meiotic_aneuploidy") {
  # aneuploidy phenotypes 
  gwas_results <- pbmclapply(1:ncol(bed_dataset), function(x) gwas_aneuploidy(x, bed_dataset, pheno, metadata, bim, pca_scores, bed_dataset_indices), mc.cores = 48L)
} else if (phenotype == "embryo_count") {
  # embryo count
  gwas_results <- pbmclapply(1:ncol(bed_dataset), function(x) gwas_embryo_count(x, bed_dataset, pheno, metadata, bim, pca_scores, bed_dataset_indices), mc.cores = 48L)
} else if (phenotype == "maternal_age") {
  # maternal age 
  gwas_results <- pbmclapply(1:ncol(bed_dataset), function(x) gwas_maternal_age(x, bed_dataset, pheno, metadata, bim, pca_scores, bed_dataset_indices), mc.cores = 48L)
}

# bind output
gwas_results_dt <- rbindlist(gwas_results[unlist(map(gwas_results, is.data.table))]) %>%
  .[!is.na(p.value)] %>%
  .[, c("snp_id", "effect_allele") := tstrsplit(snp, "_", fixed = TRUE)] %>%
  merge(., bim[, -"snp_id"], by = "pos") %>%
  setorder(., p.value)

# write to file
write.table(gwas_results_dt, out_fname, append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE, quote = FALSE)

