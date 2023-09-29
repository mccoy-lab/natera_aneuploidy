## Make phenotype file for embryos affected by maternal meiotic aneuploidy (cn = 1p or 3m for a number of chromosomes between the thresholds)

# load libraries
library(data.table)
library(tidyr)
library(dplyr)

# Usage: 
# /scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/scripts/phenotypes/maternal_meiotic_aneuploidy.R \ 
# /scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/results/phenotypes/maternal_meiotic_aneuploidy_by_mother.csv \
# mother \
# /data/rmccoy22/natera_spectrum/karyohmm_outputs/compiled_output/natera_embryos_v2.karyohmm_v14.bph_sph_trisomy.071023.tsv.gz \ 
# 2 \
# 5 \ # 5 or more chromosomes at cn=0 is considered failed amplification 
# 3 # 3 or more aneuploid chromosomes is not considered "maternal aneuploidy" but rather another ploidy (number of chromosomes greater than which the embryo is not just "aneuploid" but rather has an entire ploidy) 

# get command line arguments
args <- commandArgs(trailingOnly = TRUE)
# output file name
out_fname <- args[1]
# parent to measure phenotype
parent <- args[2]
# ploidy calls from karyohmm
embryos <- args[3]
# minimum bayes factor for filtering
bayes_factor_cutoff <- as.numeric(args[4])
# maximum number of chromosomes that are allowed nullisomies; anything more is considered failed amplification
nullisomy_threshold <- as.numeric(args[5])
# maximum number of chromosomes that are allowed aneuploidy; anything more is a whole-genome ploidy
ploidy_threshold <- as.numeric(args[6])

# source Rscript with functions `filter_data` and `count_ploidy_by_parent`
source("helper_functions/get_ploidy.R")

# read in and filter data
embryos <- fread(embryos)
embryos <- filter_data(embryos, bayes_factor_cutoff)

# remove embryos with failed amplification (5 or more nullisomies), triploidies, or haploidies 
count_nullisomies <- embryos %>% 
    group_by({{parent}}, child) %>% 
    summarise(num_nullisomies = sum(bf_max_cat == "0"))
successful_amp <- count_nullisomies[count_nullisomies$num_nullisomies < nullisomy_threshold,]
# grab triploid embryos 
count_triploidies <- embryos %>% 
    group_by({{parent}}, child) %>% 
    summarise(num_trisomies = sum(bf_max_cat == "3m" | bf_max_cat == "3p")) 
non_trip <- count_triploidies[count_triploidies$num_trisomies < ploidy_threshold,]              
# grab haploid embryos (this would also catch isoUPD embryos)
count_haploidies <- embryos %>% 
    group_by({{parent}}, child) %>% 
    summarise(num_monosomies = sum(bf_max_cat == "1m" | bf_max_cat == "1p"))
non_hap <- count_haploidies[count_haploidies$num_monosomies < ploidy_threshold,]
# remove failed amp, triploid, haploid embryos 
embryos_filtered <- embryos[embryos$child %in% successful_amp$child & embryos$child %in% non_trip$child & embryos$child %in% non_hap$child]


# count maternal meiotic aneuploidies per embryo, based on parent
# group ploidy by respective parent 
ploidy_counts_by_parent <- count_ploidy_by_parent(embryos_filtered, !!as.name(parent), phenotype, ploidy_threshold, parent)
colnames(ploidy_counts_by_parent)[1] <- "array"

# write to file 
write.csv(ploidy_counts_by_parent, out_fname, row.names = FALSE)
