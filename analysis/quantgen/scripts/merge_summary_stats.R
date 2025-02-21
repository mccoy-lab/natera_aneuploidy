# =================
# author: Sara A. Carioscia, Biology Dept., Johns Hopkins University
# email: scarios1@jhu.edu
# last update: February 13, 2025
# aim: for each lead SNP from the aneuploidy and recombination GWAS results, 
#       get the beta and p-value from all other traits and create table.
# =================

# load libraries
library(dplyr)
library(data.table)

# Read in arguments
args <- commandArgs(trailingOnly = TRUE)
primary_table <- args[1]
out_filename <- args[2]
summary_stat_files <- args[3:length(args)]

# Load table with lead SNPs from aneuploidy and recombination traits
primary_table <- read.table(primary_table, header = FALSE)
colnames(primary_table) <- c("Primary_Trait", "SNP", "Primary_Trait_EffectAllele", 
                             "Primary_Trait_Beta", "Primary_Trait_SE", "Primary_Trait_P-value")

# Iterate through each file of summary stats
for (file_path in summary_stat_files) {
  
  # Extract the trait name from the file name
  trait_name <- gsub("_summary_stats_cpra|\\.tsv$", "", basename(file_path))
  
  # Load the summary statistic file
  #summary_stats <- read.table(file_path, header = TRUE, fill = TRUE)
  summary_stats <- fread(file_path, fill = TRUE)
  
  # Remove any rows for which SNP was NA 
  summary_stats <- summary_stats[!is.na(summary_stats$SNP),]
  
  # Merge with the primary table on SNP
  primary_table <- left_join(primary_table, summary_stats %>% 
                               dplyr::select(SNP, BETA, SE, P, A1), 
                             by = "SNP", relationship = "many-to-many")
  
  # Rename the newly added columns
  colnames(primary_table)[(ncol(primary_table) - 3):ncol(primary_table)] <- c(
    paste(trait_name, "Beta", sep = "_"),
    paste(trait_name, "SE", sep = "_"),
    paste(trait_name, "P-value", sep = "_"),
    paste(trait_name, "EffectAllele", sep = "_")
    )

  # Free up memory after using each summary stats file
  rm(summary_stats)  
  gc()  
}

# Save the final table
write.table(primary_table, out_filename, col.names = TRUE, row.names = FALSE, quote = FALSE)
