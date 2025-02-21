#!/bin/bash

# =================
# author: Sara A. Carioscia, Biology Dept., Johns Hopkins University
# email: scarios1@jhu.edu
# last updated: February 10, 2025
# aim: Extract lead SNP from each peak in aneuploidy and recombination GWAS results
# =================

# Usage: extract_significant_snps.sh <aneuploidy_summary_stats> <lead_variants_recombination> <output_file>

aneuploidy_summary_stats=$1
lead_variants_recombination=$2
output_file=$3

# Ensure the output directory exists
mkdir -p "$(dirname "$output_file")"

# Extract lead SNP from aneuploidy GWAS
tail -n +2 "$aneuploidy_summary_stats" | awk 'length($2) == 1 && length($3) == 1' | sort -k6,6g | head -n 1 | awk 'BEGIN {OFS="\t"} {print "MaternalMeioticAneu", $1, $2, $4, $5, $6}' > "$output_file"

# Extract lead SNPs from recombination GWAS results
awk 'BEGIN{OFS="\t"} {print $1, $NF, $10, $8, $9, $3}' "$lead_variants_recombination" >> "$output_file"
