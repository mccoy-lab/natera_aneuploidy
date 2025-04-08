# Aneuploidy by chromosome, by parent of origin

# =================
# author: Sara A. Carioscia, Biology Dept., Johns Hopkins University
# email: scarios1@jhu.edu
# last update: March 4, 2025
# aim: plot and count number of aneuploid autosomes per embryo 
# =================

# load libraries
library(data.table)
library(dplyr)
library(ggplot2)

# Read in data (aneuploidy calls for autosomes)
ploidy_calls <- fread("/data/rmccoy22/natera_spectrum/karyohmm_outputs/compiled_output/natera_embryos.karyohmm_v30a.bph_sph_trisomy.full_annotation.031624.tsv.gz")

# Read in segmentals for filtering
segmental_calls <- fread("/scratch16/rmccoy22/abiddan1/natera_segmental/analysis/segmental_qc/results/tables/segmental_calls_postqc_refined.tsv.gz")

# Filter data to match chromosome calls used in aneuploidy GWAS

# Keep each embryo only once (one call for each chromosome)
ploidy_calls <- ploidy_calls %>%
  distinct(child, chrom, .keep_all = TRUE)

# Remove embryos that have noise more than 3sd from mean
ploidy_calls <- ploidy_calls[ploidy_calls$embryo_noise_3sd == FALSE, ]

# Remove day 3 embryos 
ploidy_calls <- ploidy_calls[ploidy_calls$day3_embryo == FALSE,]

# Remove embryos with failed amplification
# Count number of chromosomes called as nullisomies for each embryo
count_nullisomies <- ploidy_calls %>%
  group_by(mother, child) %>%
  summarise(num_nullisomies = sum(bf_max_cat == "0"))
# Identify embryos with fewer nullisomies than the threshold
successful_amp <- count_nullisomies[count_nullisomies$num_nullisomies < 5, ]
# Keep only embryos without failed amplification
ploidy_calls <- ploidy_calls[ploidy_calls$child %in% successful_amp$child, ]


# Keep only rows with bayes factor greater than the threshold for
# bayes factor qc
bayes_factor_cutoff <- 2
ploidy_calls <- ploidy_calls[ploidy_calls$bf_max > bayes_factor_cutoff, ]

# Keep only chromosomes with sufficiently high probability cn call
min_prob <- 0.9
ploidy_calls <- ploidy_calls[ploidy_calls$post_max > min_prob, ]

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


# Get number of aneuploid and disomic chromosomes per embryo
count_aneuploidies <- ploidy_calls %>% 
  group_by(mother, father, child) %>% 
  count(bf_max_cat %in% c("0", "1p", "1m", "3p", "3m"))

aneuploidies_by_embryo <- count_aneuploidies %>% 
  pivot_wider(names_from = `bf_max_cat %in% c("0", "1p", "1m", "3p", "3m")`, values_from = n)
colnames(aneuploidies_by_embryo) <- c("mother", "father", "child", "num_disomy", "num_aneuploid")
# Fill in NA with 0
aneuploidies_by_embryo$num_aneuploid[is.na(aneuploidies_by_embryo$num_aneuploid)] <- 0
aneuploidies_by_embryo$num_disomy[is.na(aneuploidies_by_embryo$num_disomy)] <- 0

# Create PDF
pdf("/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/characterize_aneuploidy/aneuploid_chr_by_embryo.pdf", height = 8, width = 12)

# Plot
ggplot(data = aneuploidies_by_embryo, aes(x = num_aneuploid)) + 
  geom_bar() + 
  labs(
    #title = "Number of aneuploid chromosomes per embryo",
    x = "Number of Aneuploid Chromosomes",
    y = "Number of Embryos Affected"
  ) +
  theme_minimal() +
  theme(legend.position = "none") + 
  theme(axis.text=element_text(size = 20),
        axis.title=element_text(size = 20))

# Disconnect
dev.off()


# Get numbers for paper 
percent_euploid <- length(which(aneuploidies_by_embryo$num_aneuploid == 0)) / nrow(aneuploidies_by_embryo)
percent_single <- length(which(aneuploidies_by_embryo$num_aneuploid == 1)) / nrow(aneuploidies_by_embryo)
percent_complex <- length(which(aneuploidies_by_embryo$num_aneuploid > 1)) / nrow(aneuploidies_by_embryo)

# Get breakdown of whole genome gain vs. loss
count_trisomies <- ploidy_calls %>%
  group_by(mother, child) %>%
  summarise(num_mat_trisomies = sum(bf_max_cat == "3m"), 
            num_pat_trisomies = sum(bf_max_cat == "3p"), 
            total_instances = n(),
            percent_mat_trisomy = num_mat_trisomies / total_instances,
            percent_pat_trisomy = num_pat_trisomies / total_instances)
count_monosomies <- ploidy_calls %>%
  group_by(mother, child) %>%
  summarise(num_mat_monosomies = sum(bf_max_cat == "1p"), 
            num_pat_monosomies = sum(bf_max_cat == "1m"), 
            total_instances = n(),
            percent_mat_monosomy = num_mat_monosomies / total_instances,
            percent_pat_monosomy = num_pat_monosomies / total_instances)
# Calculate percentages for gain/loss by parent of origin 
percent_maternal_wg_gain <- length(which(count_trisomies$percent_mat_trisomy >= 0.9)) / nrow(aneuploidies_by_embryo)
percent_paternal_wg_gain <- length(which(count_trisomies$percent_pat_trisomy >= 0.9)) / nrow(aneuploidies_by_embryo)
percent_maternal_wg_loss <- length(which(count_monosomies$percent_mat_monosomy >= 0.9)) / nrow(aneuploidies_by_embryo)
percent_paternal_wg_loss <- length(which(count_monosomies$percent_pat_monosomy >= 0.9)) / nrow(aneuploidies_by_embryo)
# Get overall percent for gain/loss 
percent_gain <- sum(percent_maternal_wg_gain + percent_paternal_wg_gain)
percent_loss <- sum(percent_maternal_wg_loss + percent_paternal_wg_loss)
# Get overall percent for whole genome
percent_wg <- percent_gain + percent_loss
