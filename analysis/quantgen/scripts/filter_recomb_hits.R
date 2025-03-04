# =================
# author: Sara A. Carioscia, Biology Dept., Johns Hopkins University
# email: scarios1@jhu.edu
# last update: January 22, 2025
# aim: starting from lead SNP in each recombination hit, ensure that no two SNPs for the same trait are within 1MB of one another.
#      If so, keep only the most significant
# =================

# load libraries
library(dplyr)

# Read in arguments 
args <- commandArgs(trailingOnly = TRUE)
input_data <- args[1]
out_fname <- args[2]

# Load table 
recomb_hits <- read.table(input_data, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Filter to only genome-wide significant 
recomb_hits <- recomb_hits[recomb_hits$Bonferroni == "true",]

# Filter to only female traits
recomb_hits_female <- recomb_hits[grepl("Female", recomb_hits$PHENO),]

# Function to process the dataset
filter_significant_entries <- function(data) {
  # Split chromosome and position from CPRA 
  data <- data %>%
    mutate(
      chr = sub(":.*", "", ID),
      position = as.numeric(sub(".*:(\\d+):.*", "\\1", ID))
    )
  
  # Group by trait and chromosome
  data <- data %>%
    group_by(trait = PHENO, chr) %>%
    arrange(position, pvalue = as.numeric(P)) %>% # Sort by position and p-value
    mutate(
      keep = c(TRUE, diff(position) > 1e6 | lead(diff(position), default = 1e6) > 1e6) # Check if within 1MB
    ) %>%
    ungroup()
  
  # Filter to keep only the most significant entries (lowest p-value)
  data <- data %>%
    filter(keep) %>%
    select(-keep, -chr, -position, -trait)
  
  return(data)
}

# Apply function
filtered_data <- filter_significant_entries(recomb_hits_female)

# Save filtered dataset
write.table(filtered_data, out_fname, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
