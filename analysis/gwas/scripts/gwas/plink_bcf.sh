#!/bin/bash

# Load PLINK 1.9 module
module load plink/1.90b6.4

bcf_file=$1
outfix=$2

plink --memory 3000 --bcf $bcf_file --double-id --allow-extra-chr --make-bed --out $outfix
