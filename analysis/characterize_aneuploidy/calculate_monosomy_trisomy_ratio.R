# Number of monosomies and trisomies 

# =================
# author: Sara A. Carioscia, Biology Dept., Johns Hopkins University
# email: scarios1@jhu.edu
# last update: February 24, 2025
# aim: plot proportion of each chromosome that is aneuploid, by parent of origin 
# =================

# load libraries
library(data.table)
library(dplyr)
library(ggplot2)

# Filter data exactly as in phenotyping script

ploidy_calls <- fread("/data/rmccoy22/natera_spectrum/karyohmm_outputs/compiled_output/natera_embryos.karyohmm_v30a.bph_sph_trisomy.full_annotation.031624.tsv.gz")
#sex_chr <- fread("/data/rmccoy22/natera_spectrum/karyohmm_outputs/compiled_output/natera_embryos.karyohmm_v30a.sex_embryos.012425.tsv.gz")

# Read in segmentals for filtering
segmental_calls <- fread("/scratch16/rmccoy22/abiddan1/natera_segmental/analysis/segmental_qc/results/tables/segmental_calls_postqc_refined.tsv.gz")

# Save raw data for comparisons 
raw_data <- ploidy_calls

# Keep each embryo only once (one call for each chromosome)
ploidy_calls <- ploidy_calls %>%
  distinct(child, chrom, .keep_all = TRUE)

# Remove embryos that have noise more than 3sd from mean
ploidy_calls <- ploidy_calls[ploidy_calls$embryo_noise_3sd == FALSE, ]

# Remove day 3 embryos 
filter_day_5 == TRUE
if (filter_day_5 == TRUE) {
  ploidy_calls <- ploidy_calls[ploidy_calls$day3_embryo == FALSE,]
}

# Remove embryos with failed amplification
# Count number of chromosomes called as nullisomies for each embryo
count_nullisomies <- ploidy_calls %>%
  group_by(mother, child) %>%
  summarise(num_nullisomies = sum(bf_max_cat == "0"))
# Identify embryos with fewer nullisomies than the threshold
successful_amp <- count_nullisomies[count_nullisomies$num_nullisomies < 5, ]

# Keep only embryos without failed amplification
ploidy_calls <- ploidy_calls[ploidy_calls$child %in% successful_amp$child, ]


# Keep only chrom that have probabilities for all 6 cn states
ploidy_calls <- ploidy_calls[complete.cases(
  ploidy_calls[,c("0", "1m", "1p", "2", "3m", "3p")]), ]

# Keep only rows with bayes factor greater than the threshold for 
# bayes factor qc
ploidy_calls <- ploidy_calls[ploidy_calls$bf_max > 2, ]

# Keep only chromosomes with sufficiently high probability cn call
ploidy_calls <- ploidy_calls[ploidy_calls$post_max > 0.9, ]

# Keep only chromosomes that are not affected by post-QC segmental aneu
# Add column that makes ID of mother-father-child-chrom in whole chr 
ploidy_calls$uid <- paste0(ploidy_calls$mother, "+", ploidy_calls$father, "+",
                           ploidy_calls$child, "+", ploidy_calls$chrom)
# Add column that makes ID of mother-father-child-chrom in segmental
segmental_calls$uid <- paste0(segmental_calls$mother, "+", 
                              segmental_calls$father, "+", 
                              segmental_calls$child, "+", 
                              segmental_calls$chrom)
# Remove from ploidy calls any chromosomes in segmentals 
ploidy_calls <- ploidy_calls[!ploidy_calls$uid %in% segmental_calls$uid,]

# Remove trisomies that are likely mosaic 
ploidy_calls <- ploidy_calls[!(ploidy_calls$bf_max_cat %in% c("3m", "3p") & 
                                 (ploidy_calls$post_bph_centro < 0.340 & 
                                    ploidy_calls$post_bph_noncentro < 0.340)), ]


# Hold data as equivalent after filter_data function from phenotyping 
ploidy_calls_filtered <- ploidy_calls 


# address triploidies and haploidies 
count_trisomies <- ploidy_calls %>%
  group_by(mother, child) %>%
  summarise(num_trisomies = sum(bf_max_cat %in% c("3m", "3p")))
count_monosomies <- ploidy_calls %>%
  group_by(mother, child) %>%
  summarise(num_monosomies = sum(bf_max_cat %in% c("1m", "1p")))
non_trip <- count_trisomies[count_trisomies$num_trisomies < 15,]
non_hap <- count_monosomies[count_monosomies$num_monosomies < 15,]
non_wg <- rbind(non_trip, non_hap)
# keep only the ones in non trip and non hap 
ploidy_calls_no_wg <- ploidy_calls[ploidy_calls$child %in% non_wg$child]


# calculate ratio of monosomy to trisomy 
num_monosomy <- length(which(ploidy_calls_no_wg$bf_max_cat %in% c("1p", "1m")))
num_trisomy <- length(which(ploidy_calls_no_wg$bf_max_cat %in% c("3p", "3m")))
binom.test(num_monosomy, num_monosomy + num_trisomy)

"""
data:  num_monosomy and num_monosomy + num_trisomy
number of successes = 34511, number of trials = 92485, p-value < 2.2e-16
alternative hypothesis: true probability of success is not equal to 0.5
95 percent confidence interval:
 0.3700340 0.3762787
sample estimates:
probability of success 
             0.3731524 

"""

