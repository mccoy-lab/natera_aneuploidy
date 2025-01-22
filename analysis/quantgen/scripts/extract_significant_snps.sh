#!/bin/bash
# Usage: extract_significant_snps.sh <aneuploidy_summary_stats> <lead_variants_recombination> <output_file>

aneuploidy_summary_stats=$1
lead_variants_recombination=$2
output_file=$3

# Ensure the output directory exists
mkdir -p "$(dirname "$output_file")"

# Extract lead SNP from aneuploidy GWAS
tail -n +2 "$aneuploidy_summary_stats" | sort -k8,8g | head -n 1 | awk '{print "MaternalMeioticAneu", $9, $6, $8}' > "$output_file"

# Extract lead SNPs from recombination GWAS results
awk 'BEGIN{OFS="\t"} NR > 1 {print $1, $NF, $8, $3}' "$lead_variants_recombination" >> "$output_file"
