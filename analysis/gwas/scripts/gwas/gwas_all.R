## GWAS for phenotypes from Natera 

# load libraries
library(data.table)
library(BEDMatrix)
library(dplyr)
library(pbmcapply)
library(purrr)

# Usage: /scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/scripts/gwas/gwas_all.R \
# "/data/rmccoy22/natera_spectrum/data/summary_metadata/spectrum_metadata_merged.csv" \
# "/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/results/gwas/opticall_concat_9.norm.b38.bed" \
# "/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/results/discover_validate_split_mother.txt" \
# "/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/results/parental_genotypes_pcs/parental_genotypes.eigenvec" \
# "/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/results/phenotypes/haploidy_by_mother.csv" \
# "/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/results/gwas/opticall_concat_9.norm.b38.bim" \
# "discovery" \ # dataset type 
# "haploidy" \ # phentype 
# "mother" \ # parent 
# 16 \ # number of threads 
# "/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/results/gwas/gwas_haploidy_by_mother_discovery_9.tsv"

# get command line arguments 
args <- commandArgs(trailingOnly = TRUE)
# metadata for all years 
metadata <- args[1]
# parental genotypes as described https://github.com/mccoy-lab/natera_spectrum/tree/main/genotyping
bed <- args[2]
# discovery-test split for applicable parent
discovery_test <- args[3] 
# principal components for parents 
pcs <- args[4]
# phenotype file of interest
phenotype <- args[5]
# file containing snp name and map positions, corresponds with bed file 
bim <- args[6]
# string variable noting whether to GWAS using the discovery or test set for parents 
dataset_type <- args[7] 
# string variable naming the phenotype being considered (e.g., maternal_meiotic)
phenotype_name <- args[8] 
# string variable noting which parent to GWAS 
parent <- args[9]
# number of threads to use in mc.cores for GWAS 
threads <- as.numeric(args[10]) 
# output file name 
out_fname <- args[11]

# source Rscript with functions ``
source("/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/scripts/gwas/helper_functions/gwas_helper_functions.R")

# read in files 
metadata <- fread(metadata)
bed <- BEDMatrix(bed)
discovery_test <- fread(discovery_test)
pcs <- fread(pcs)
colnames(pcs)[1] <- "array"
pheno <- fread(phenotype)
bim <- fread(bim) %>% 
  setnames(., c("chr", "snp_id", "drop", "pos", "ref", "alt"))



# function to run GWAS on each site if it's a ploidy phenotype 
# (triploid, haploid, maternal meiotic)
gwas_aneuploidy <- function(snp_index, genotypes, phenotypes, metadata, locs, 
                            pcs, subject_indices, parent) {
    
    gt <- get_gt(snp_index, genotypes, phenotypes, metadata, locs, 
                 pcs, subject_indices)
    
    
    snp_name <- colnames(genotypes)[snp_index]
    snp_chr <- locs[snp_index]$chr
    snp_pos <- locs[snp_index]$pos
    
    
    
    coef <- data.table(term = rownames(m1$coefficients), m1$coefficients)
    
    gt_filtered  <- gt %>% filter(!is.na(alt_count))
    alt_af <- sum(gt_filtered$alt_count) / (2*nrow(gt_filtered))
    
    return(data.table(snp = snp_name, 
                      pos = snp_pos, 
                      beta = unlist(coef[term == "alt_count", 2]),
                      se = unlist(coef[term == "alt_count", 3]),
                      t = unlist(coef[term == "alt_count", 4]),
                      p.value = unlist(coef[term == "alt_count", 5]),
                      af = alt_af)) 
}

# function to run GWAS on each site if it's embryo count 
gwas_embryo_count <- function(snp_index, genotypes, phenotypes, metadata, locs, 
                              pcs, subject_indices) {
    
    gt <- get_gt(snp_index, genotypes, phenotypes, metadata, locs, 
                 pcs, subject_indices)
    
    m1 <- glm(num_embryos ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + 
                PC9 + PC10 + weighted_age + alt_count, family = "poisson", 
              data = gt) %>% summary()
    
    
    
    
    coef <- data.table(term = rownames(m1$coefficients), m1$coefficients)
    
    gt_filtered  <- gt %>% filter(!is.na(alt_count))
    alt_af <- sum(gt_filtered$alt_count) / (2*nrow(gt_filtered))
    
    return(data.table(snp = snp_name, 
                      pos = snp_pos, 
                      beta = unlist(coef[term == "alt_count", 2]),
                      se = unlist(coef[term == "alt_count", 3]),
                      t = unlist(coef[term == "alt_count", 4]),
                      p.value = unlist(coef[term == "alt_count", 5]),
                      af = alt_af))
}

# function to run GWAS on each site for maternal age 
gwas_maternal_age <- function(snp_index, genotypes, phenotypes, metadata, locs, 
                              pcs, subject_indices) {
    
    gt <- get_gt(snp_index, genotypes, phenotypes, metadata, locs, 
                 pcs, subject_indices)
    
    gt_filtered  <- gt %>% filter(!is.na(alt_count))
    alt_af <- sum(gt_filtered$alt_count) / (2*nrow(gt_filtered))
    
    m1 <- glm(weighted_age ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + 
                PC9 + PC10 + alt_count, family = "poisson", data = gt) %>%
      summary()
    
    coef <- data.table(term = rownames(m1$coefficients), m1$coefficients)
    
    return(data.table(snp = snp_name, 
                      pos = snp_pos, 
                      beta = unlist(coef[term == "alt_count", 2]),
                      se = unlist(coef[term == "alt_count", 3]),
                      t = unlist(coef[term == "alt_count", 4]),
                      p.value = unlist(coef[term == "alt_count", 5]),
                      af = alt_af))  
}

# run GWAS based on phenotype passed argument 
if (phenotype_name %in% c("maternal_meiotic_aneuploidy", 
                          "maternal_triploidy", "paternal_triploidy", 
                          "maternal_haploidy", "paternal_haploidy")) {
    # aneuploidy phenotypes 
    gwas_results <- 
      pbmclapply(1:ncol(bed_dataset), 
                 function(x) gwas_aneuploidy(x, bed_dataset, pheno, metadata, 
                                             bim, pca_scores, 
                                             bed_dataset_indices, parent), 
                 mc.cores = threads)
} else if (phenotype_name == "embryo_count") {
    # embryo count
    gwas_results <- 
      pbmclapply(1:ncol(bed_dataset), 
                 function(x) gwas_embryo_count(x, bed_dataset, pheno, metadata, 
                                               bim, pca_scores, 
                                               bed_dataset_indices), 
                 mc.cores = threads)
} else if (phenotype_name == "maternal_age") {
    # maternal age 
    gwas_results <- 
      pbmclapply(1:ncol(bed_dataset), 
                 function(x) gwas_maternal_age(x, bed_dataset, pheno, metadata, 
                                               bim, pca_scores, 
                                               bed_dataset_indices), 
                 mc.cores = threads)
}

# bind output
gwas_results_dt <- 
  rbindlist(gwas_results[unlist(map(gwas_results, is.data.table))]) %>%
    .[!is.na(p.value)] %>%
    .[, c("snp_id", "effect_allele") := tstrsplit(snp, "_", fixed = TRUE)] %>%
    merge(., bim[, -"snp_id"], by = "pos") %>%
    setorder(., p.value)

# write to file
write.table(gwas_results_dt, out_fname, append = FALSE, sep = "\t", dec = ".",
            row.names = TRUE, col.names = TRUE, quote = FALSE)

