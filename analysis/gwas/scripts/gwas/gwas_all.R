## GWAS for maternal meiotic errors + parental genotypes 

# load libraries
library(data.table)
library(BEDMatrix)
library(dplyr)
library(pbmcapply)
library(purrr)

# Usage: ./gwas_all.R \
# "/data/rmccoy22/natera_spectrum/data/summary_metadata/spectrum_metadata_merged.csv" \
# "/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/results/gwas/opticall_concat_9.norm.b38.bed" \
# "/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/results/discover_validate_split_mother.txt" \
# "/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/results/parental_genotypes_pcs/parental_genotypes.eigenvec" \
# "/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/results/phenotypes/haploidy_by_mother.csv" \
# "/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/results/gwas/opticall_concat_9.norm.b38.bim" \
# "discovery" \
# "haploidy" \
# "mother" \
# 16 / number of threads 
# "/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/results/gwas/gwas_haploidy_by_mother_discovery_9.tsv"

# take in arguments 
args <- commandArgs(trailingOnly = TRUE)
# metadata for all years 
metadata <- args[1]
# parental genotypes as described https://github.com/mccoy-lab/natera_spectrum/tree/main/genotyping
bed <- args[2]
# discovery-test split for applicable parent
discovery_test <- args[3] 
# PCs for parents, output from plink 
pca_scores <- args[4]
# phenotype file of interest
pheno <- args[5]
# file containing snp name and map positions, corresponds with bed file 
bim <- args[6]
# string variable noting whether to GWAS using the discovery or test set for parents 
dataset_type <- args[7] 
# string variable naming the phenotype being considered (e.g., aneuploidy type, embryo count)
phenotype_name <- args[8] 
# string variable noting which parent to GWAS 
parent <- args[9]
# number of threads to use in mc.cores for GWAS 
threads <- as.numeric(args[10]) 
# output file name 
out_fname <- args[11]


# Check if input files exist
stopifnot(file.exists(metadata), "Metadata file not found")
stopifnot(file.exists(bed), "BED file not found")
stopifnot(file.exists(discovery_test), "Discovery test file not found")
stopifnot(file.exists(pca_scores), "PCA scores file not found")
stopifnot(file.exists(pheno), "Phenotype file not found")
stopifnot(file.exists(bim), "BIM file not found")

# Check valid dataset_type, phenotype_name, and parent
stopifnot(dataset_type %in% c("discovery", "test"), "Invalid dataset_type")
stopifnot(phenotype_name %in% c("maternal_meiotic_aneuploidy", "maternal_triploidy", "paternal_triploidy", "maternal_haploidy", "paternal_haploidy", "embryo_count", "maternal_age"), "Invalid phenotype_name")
stopifnot(parent %in% c("mother", "father"), "Invalid parent")

# Check valid threads value
stopifnot(is.numeric(threads) && threads >= 1, "Invalid number of threads")

# Check valid out_fname
stopifnot(!is.null(out_fname) && nchar(out_fname) > 0, "Invalid output file name")


# read in files from arguments 
metadata <- fread(metadata)
bed <- BEDMatrix(bed)
discovery_test <- fread(discovery_test)
pca_scores <- fread(pca_scores) 
colnames(pca_scores)[1] <- "array"
pheno <- fread(pheno)
bim <- fread(bim) %>% 
  setnames(., c("chr", "snp_id", "drop", "pos", "ref", "alt"))

"
For testing
genotypes <- bed
phenotypes <- pheno
locs <- bim
pcs <- pca_scores
subject_indices <- bed_dataset_indices
"

# separate into discovery or test 
discovery_test_split <- function(dataset_type, metadata, bed, discovery_test) {
    if (dataset_type == "discovery") {
        dataset <- discovery_test[discovery_test$is_discovery == TRUE,]
    } else if (dataset_type == "test") {
        dataset <- discovery_test[discovery_test$is_discovery == FALSE,]
    }
    # get caseIDs belonging to mothers in preferred set
    metadata_set <- metadata[metadata$array %in% dataset$array,]
    # add column for array_array to match rownames of bed
    metadata_set$array_array <- paste0(metadata_set$array, "_", 
                                       metadata_set$array)
    # subset bed file to just the mothers in the relevant set 
    bed_indices <- which(rownames(bed) %in% metadata_set$array_array)
    bed_dataset <- bed[bed_indices,]
    # return relevant subset of bed file (discovery vs. test)
    return(bed_dataset)
}

# get bed file corresponding to correct dataset type 
bed_dataset <- discovery_test_split(dataset_type, metadata, bed, discovery_test)

# get variables to call function 
bed_dataset_indices <- 1:nrow(bed_dataset)

# function to pre-process gt info for each site 
get_gt <- function(snp_index, genotypes, phenotypes, metadata, locs, 
                   pcs, subject_indices) {
    gt <- data.table(names(genotypes[subject_indices, snp_index]), 
                     unname(genotypes[subject_indices, snp_index])) %>%
        setnames(., c("array", "alt_count")) %>%
        .[, array := sub("(.*?_.*?)_.*", "\\1", array)]
    
    gt <- merge(gt, metadata, by = "array") %>%
        merge(phenotypes, by = "array") %>%
        merge(pcs, by = "array") %>%
        .[!duplicated(array)]
    
    return(gt)
}
# the genotypes info has each arrayID pasted twice, 
# e.g., 10005770025_R06C01_10005770025_R06C01
# the reg ex keeps everything before the second underscore, 
# returning the array info to its standard format, e.g., 10005770025_R06C01

# function to run GWAS on each site if it's a ploidy phenotype 
# (triploid, haploid, maternal meiotic aneuploidy, parental_triploidy)
gwas_aneuploidy <- function(snp_index, genotypes, phenotypes, metadata, locs, 
                            pcs, subject_indices, parent) {
    
    gt <- get_gt(snp_index, genotypes, phenotypes, metadata, locs, 
                 pcs, subject_indices)
    
    # Determine the age covariate column based on the "parent" argument
    if (parent == "mother") {
        age_column <- "patient_age"
    } else if (parent == "father") {
        age_column <- "partner_age"
    } 
    
    snp_name <- colnames(genotypes)[snp_index]
    snp_chr <- locs[snp_index]$chr
    snp_pos <- locs[snp_index]$pos
    
    formula_string <- paste0("cbind(aneu_true, aneu_false) ~ PC1 + PC2 + ",
                             "PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + ", 
                             age_column, " + alt_count", collapse = "")
    m1 <- glm(formula_string, family = "quasibinomial", data = gt) %>% summary()
    
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

