#!/bin/bash

#SBATCH --job-name=bed_split_vcfs
#SBATCH -N 1
#SBATCH --time=1:00:00
#SBATCH --partition=shared
#SBATCH --mem=1G
#SBATCH --mail-type=fail
#SBATCH --mail-user=scarios1@jhu.edu

module load bcftools 
module load plink

echo ${SLURM_ARRAY_TASK_ID}

# read in list of map files
map_file=$1
# get input vcf to split 
input_vcf=$2
# get output dir to write plink files 
outdir=$3

# get name of map file to use as vcf subsetter
line_num=${SLURM_ARRAY_TASK_ID}
map_name=$(sed "${line_num}q;d" "$map_file")

# get file prefix for newly output plink files
prefix=$(basename "$map_name" .txt)

# use map to make a newly subsetted bcf
bcftools view -T $map_name -Ob $input_vcf > "${outdir}/${prefix}.bcf"

# run plink to convert to bed/bim/fam
plink --bcf "${outdir}/${prefix}.bcf" --double-id --allow-extra-chr --make-bed --out "${outdir}/${prefix}"
