## GWAS for maternal meiotic errors + maternal genotypes 

# load libraries
library(data.table)
library(BEDMatrix)
library(dplyr)
library(pbmcapply)
library(purrr)

# take in arguments 
args <- commandArgs(trailingOnly = TRUE)
eigenvec <- args[1]
bed <- args[2]
bim <- args[3]
pheno <- args[4]
metadata <- args[5]
discovery_test_f <- args[6]
out_fname <- args[7]

# read in files from arguments 
eigenvec <- fread(eigenvec)
bed <- BEDMatrix(bed)
bim <- fread(bim) %>% setnames(., c("chr", "snp_id", "drop", "pos", "ref", "alt"))
pheno <- fread(pheno, sep = ",")
metadata <- fread(metadata)
discovery_test_f <- fread(discovery_test_f)

# separate discovery and validation sets 
discovery_f <- discovery_test_f[discovery_test_f$is_discovery == TRUE,]
# get caseIDs belonging to mothers in discovery set
metadata_f_discovery <- metadata[metadata$array %in% discovery_f$array,]
# add column for array_array to match rownames of bed 
metadata_f_discovery$array_array<- paste0(metadata_f_discovery$array, "_", metadata_f_discovery$array)
# subset bed file to just the mothers in the discovery set 
bed_discovery_indices <- which(rownames(bed) %in% metadata_f_discovery$array_array)
bed_discovery_f <- bed[bed_discovery_indices,]
# get variables to call function 
bed_indices_discovery <- 1:nrow(bed_discovery_f)

# subset bed 
set.seed(1)
bed_subset <- bed_discovery_f[, sample(1:nrow(bim), 10000)]
bed_subset <- bed_discovery_f

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
gwas_results_discovery <- pbmclapply(1:ncol(bed_discovery_f), function(x) gwas(x, bed_discovery_f, pheno, metadata, bim, pca_scores, bed_indices_discovery), mc.cores = 48L)

# bind output
gwas_results_dt <- rbindlist(gwas_results_discovery[unlist(map(gwas_results_discovery, is.data.table))]) %>%
  .[!is.na(p.value)] %>%
  .[, c("snp_id", "effect_allele") := tstrsplit(snp, "_", fixed = TRUE)] %>%
  merge(., bim[, -"snp_id"], by = "pos") %>%
  setorder(., p.value)

# write to file
write.table(gwas_results_dt, out_fname, append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE, quote = FALSE)

