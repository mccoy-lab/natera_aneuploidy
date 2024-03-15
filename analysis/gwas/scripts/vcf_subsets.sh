#!/bin/bash

#SBATCH --job-name=make_vcf_maps
#SBATCH -N 1
#SBATCH --time=1:00:00
#SBATCH --partition=shared
#SBATCH --mem=1G
#SBATCH --mail-type=fail
#SBATCH --mail-user=scarios1@jhu.edu

# Usage: sbatch vcf_subsets.sh spectrum_imputed_chr21_rehead_filter_cpra.vcf.gz /scratch16/rmccoy22/scarios1/sandbox/spectrum_imputed_chr21_subsets

module load bcftools

# load vcf from argument 
input_vcf=$1
# load outdir from argument 
outdir=$2

# get outfile prefix 
prefix=$(basename "$input_vcf" .vcf.gz)

# create map files for each chunk of the input vcf 
bcftools query -f'%CHROM\t%POS\n' $input_vcf | split -l 5000 -d --additional-suffix=".txt" - "${outdir}/${prefix}"

# make list of all resulting map files 
ls "${outdir}/${prefix}"*".txt" > "${outdir}/mapfiles_${prefix}.txt"

# count number of map files 
n=$(wc -l < "${outdir}/mapfiles_${prefix}.txt")

# submit array job to run bcftools and plink on each map file 
sbatch --array=1-$n plink_subsets.sh "${outdir}/mapfiles_${prefix}.txt" "$input_vcf" "$outdir"