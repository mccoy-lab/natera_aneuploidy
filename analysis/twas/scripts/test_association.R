#!/usr/bin/env Rscript

# test_association.R
# =================
# author: Margaret R. Starostik, Biology Dept., Johns Hopkins University
# email: mstaros1@jhu.edu
# last updated: October 21, 2024
# testing purposes: Rscript --vanilla /home/mstaros1/scr16_rmccoy22/mstarostik/natera_aneuploidy/code/test_association.R /home/mstaros1/scr16_rmccoy22/mstarostik/natera_aneuploidy/analysis/twas/predict_expression/Artery_Aorta_chr10_predict.txt /home/mstaros1/scr16_rmccoy22/mstarostik/natera_aneuploidy/analysis/twas/association/single_tissue/Artery_Aorta_chr10_association.txt
# =================

############################################################################
# AUXILIARY FUNCTIONS
############################################################################

# load packages
library(data.table)
library(dplyr)
library(tidyr)

# association testing
test_association <- function(merged, genes) {
      # initialize association dataframe
      association.df <- NULL
      
      # perform association test between each gene and phenotype
      for (gene in genes) {
            gene.of.interest <- filtered_expression_matrix %>% 
                                    select(all_of(c("FID", gene))) %>%
                                    rename(predicted_gene_expression = gene)

      # merge phenotype and expression data
      merged <- inner_join(phenotypes, principal.components, by = c("mother" = "#IID")) %>%
                  left_join(., gene.of.interest, by = c("mother" = "FID"))

      # of interest: maternal meiotic aneuploidy
      ## note: model is from what was used in GWAS (https://github.com/mccoy-lab/natera_aneuploidy/blob/main/analysis/gwas/scripts/gwas/gwas_lmm_chunks.R)
      pcs <- paste(paste0("PC", 1:20), collapse = " + ")
      model.formula <- paste0("cbind(aneu_true, aneu_false) ~ ", 
                              pcs,
                               " + (egg_donor == 'yes')",
                              " + (sperm_donor == 'yes')",
                              " + scale(patient_age_cycle)",
                              " + scale(partner_age_cycle)",
                              " + predicted_gene_expression")
    
      model <- glm(model.formula, 
                  family = binomial,
                  data = merged)

      # extract twas results
      extract.results <- c(summary(model)$coefficients["predicted_gene_expression", "Estimate"],
                          summary(model)$coefficients["predicted_gene_expression", "z value"],
                          summary(model)$coefficients["predicted_gene_expression", "Pr(>|z|)"], 
                          summary(model)$coefficients["predicted_gene_expression", "Std. Error"])
      association.results <- c(gene, extract.results)
      association.df <- as.data.frame(rbind(association.df, association.results))
  }

      # specify column names of association dataframe
      colnames(association.df) <- c("gene", "beta", "zval", "pval", "se")

       # merge association dataframe with reference annotations
      results <- right_join(gtf, 
                            association.df,
                            by = c("gene_id" = "gene"))

      # write twas results to file
      association.filename <- args[2]
      #association.filename <- "/home/mstaros1/scr16_rmccoy22/mstarostik/natera_aneuploidy/analysis/twas/association/singel_tissue/Artery_Aorta_chr10_association.txt"
      write.table(results,
                  file = association.filename,
                  sep = "\t",
                  row.names = FALSE,
                  col.names = TRUE,
                  quote = FALSE)
}

#############################################################################
# MAIN FUNCTIONS
#############################################################################

##############
# read in and prepare data
##############

# phenotypes file for maternal meiotic aneuploidy
## note: Sara assigned egg and sperm donor ages in her GWAS based on:
## (1) mean egg donor age from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7530253/
## (2) mean sperm donor age from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9118971/#:~:text=Donors%20were%20aged%2027%20years,aged%2030%20years%20and%20younger
phenotypes <- fread("/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/results/phenotypes/maternal_meiotic_aneuploidy_by_mother.csv") %>%
              mutate(patient_age_cycle = ifelse(egg_donor == "yes", as.numeric(25), patient_age_cycle),
                    partner_age_cycle = ifelse(sperm_donor == "yes", as.numeric(27), partner_age_cycle))

# principal components
principal.components <- fread("/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/results/parental_genotypes_pcs/parental_genotypes.eigenvec")

# Gencode v26 gene-level annotations that was used in GTEx v8
gtf <- fread("/scratch16/rmccoy22/mstarostik/natera_aneuploidy/reference/gencode.v26.annotation.gtf",
            skip = 8,
            sep = "\t",
            quote = "",
            col.names = c("seqname","source","feature","start","end","score","strand","frame","attribute"),
            data.table = FALSE) %>%
            filter(feature == "gene") %>%
        separate(attribute, 
                into = c("gene_id", "gene_type", "gene_name"), 
                sep = "; ", 
                extra = "drop") %>%
        mutate(gene_id = gsub('gene_id |"', '', gene_id),
              gene_type = gsub('gene_type |"', '', gene_type),
              gene_name = gsub('gene_name |"', '', gene_name)) %>%
        select(-c("source", "feature", "score", "frame"))

# expression matrix
args <- commandArgs(trailingOnly = TRUE)
expression.filename <- args[1]
#expression.filename <- "/home/mstaros1/scr16_rmccoy22/mstarostik/natera_aneuploidy/analysis/twas/predict_expression/Artery_Aorta_chr10_predict.txt"
expression_matrix <- fread(expression.filename)

## keep only protein-coding genes?
## remove genes with zero variance
filtered_expression_matrix <- expression_matrix %>%
select_if(function(x) is.character(x) | is.numeric(x) && var(x) != 0)
## remove genes not expressed in XX individuals?
## normalize gene expression matrix?
## scale gene expression matrix values?

##############
# perform association tests and write results to file
##############
genes <- colnames(filtered_expression_matrix)[c(-1, -2)]
test_association(merged, genes)