# Aneuploidy by chromosome, by parent of origin

# =================
# author: Sara A. Carioscia, Biology Dept., Johns Hopkins University
# email: scarios1@jhu.edu
# last update: January 21, 2025
# aim: plot proportion of each chromosome that is aneuploid, by parent of origin 
# =================

# load libraries
library(data.table)
library(dplyr)
library(ggplot2)

# Read in data (aneuploidy calls for autosomes, sex chromosomes)
ploidy_calls <- fread("/data/rmccoy22/natera_spectrum/karyohmm_outputs/compiled_output/natera_embryos.karyohmm_v30a.bph_sph_trisomy.full_annotation.031624.tsv.gz")
sex_chr <- fread("/data/rmccoy22/natera_spectrum/karyohmm_outputs/compiled_output/natera_embryos.karyohmm_v30a.sex_embryos.012425.tsv.gz")

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

# Remove embryos with whole-genome gain or loss
count_trisomies <- ploidy_calls %>%
  group_by(mother, child) %>%
  summarise(num_trisomies = sum(bf_max_cat %in% c("3m", "3p")))
count_monosomies <- ploidy_calls %>%
  group_by(mother, child) %>%
  summarise(num_monosomies = sum(bf_max_cat %in% c("1m", "1p")))
# Identify embryos with fewer mono/trisomies than the threshold for 
# maternal meiotic aneuploidy 
non_triploid <- count_monosomies[count_monosomies$num_monosomies <= 5, ]
non_haploid <- count_trisomies[count_trisomies$num_trisomies <= 5, ]
# Keep only embryos without whole-genome gain or loss
ploidy_calls <- ploidy_calls[ploidy_calls$child %in% non_triploid$child, ]
ploidy_calls <- ploidy_calls[ploidy_calls$child %in% non_haploid$child, ]

# Keep only chrom that have probabilities for all 6 cn states
ploidy_calls <- ploidy_calls[complete.cases(
  ploidy_calls[,c("0", "1m", "1p", "2", "3m", "3p")]), ]

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


# Get number of children present after filtering, for determining proportion
num_embryos <- length(unique(ploidy_calls$child))

# Write to autosomes for downstream use 
autosomes <- ploidy_calls


# Categorize calls for sex chromosomes

# Remove embryos that were already removed from autosomal calls 
# (e.g., trip/haploid, failed amplification, day 3)
sex_chr <- sex_chr[sex_chr$child %in% ploidy_calls$child]

# Create column with max probability for filtering 
sex_chr$max_prob <- apply(sex_chr[, 6:11], 1, max)
sex_chr <- sex_chr[sex_chr$max_prob > min_prob]
sex_chr$max_prob_Y <- apply(sex_chr[, 16:17], 1, max)
sex_chr <- sex_chr[sex_chr$max_prob_Y > min_prob]

# Combine x and y max BF cats to get all possible unique combinations
sex_chr$karyotype <- paste0(sex_chr$x_maxBFcat, "_", sex_chr$y_maxBFcat)

# Reformat to show maternal/paternal and gain/loss
sex_chr$status <- ifelse(sex_chr$karyotype == "x0_y1", "loss_maternal",
                         ifelse(sex_chr$karyotype == "x1m_y0", "loss_paternalY",
                                ifelse(sex_chr$karyotype == "x1p_y0", "loss_maternal",
                                       ifelse(sex_chr$karyotype == "x2_y1", "gain_paternal",
                                              ifelse(sex_chr$karyotype == "x3_y0", "gain_maternal",
                                                     sex_chr$karyotype)))))


# Format calls for plotting 

# Sort chrom numerically by removing 'chr' prefix
autosomes$chrom <- factor(autosomes$chrom, 
                          levels = unique(autosomes$chrom[order(as.numeric(gsub("chr", "", autosomes$chrom)))]))

# Set the order of bf_max_cat to control stacking
autosomes$bf_max_cat <- factor(autosomes$bf_max_cat, levels = c("3m", "1p", "3p", "1m"))

# Reformat autosomal data to allow combination with sex chr
autosomes_clean <- autosomes %>%
  filter(bf_max_cat %in% c("1m", "1p", "3m", "3p")) %>%
  mutate(
    chrom = gsub("chr", "", chrom),  # Remove 'chr' prefix
    parent_origin = case_when(
      bf_max_cat %in% c("1p", "3m") ~ "maternal",
      bf_max_cat %in% c("1m", "3p") ~ "paternal"
    )
  )

# Consolidate maternal/paternal in sex chromosomes across different call states
sex_chr_clean <- sex_chr %>%
  filter(status %in% c("loss_maternal", "loss_paternal", "gain_maternal", "gain_paternal", "loss_paternalY")) %>%
  mutate(
    chrom = case_when(
      status %in% c("gain_maternal", "loss_maternal", "gain_paternal") ~ "X",
      status %in% c("loss_paternalY") ~ "Y",
      TRUE ~ NA_character_
    ),
    parent_origin = case_when(
      status %in% c("gain_maternal", "loss_maternal") ~ "maternal",
      status %in% c("gain_paternal", "loss_paternalY") ~ "paternal"
    )
  )

# Combine autosomes and sex chromosomes
combined_data <- bind_rows(autosomes_clean, sex_chr_clean)

# Ensure chrom order is numeric for autosomes, then X and Y
combined_data <- combined_data %>%
  mutate(chrom = factor(chrom, levels = c(as.character(1:22), "X", "Y")))

# Aggregate data to get counts and proportions
combined_data <- combined_data %>%
  group_by(chrom, parent_origin) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(proportion = count / num_embryos)  # Calculate the proportion


# Plot 

# Create PDF
pdf("/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/characterize_aneuploidy/parent_aneuploidy_by_chr.pdf", height = 8, width = 12)

# Plot combined data as proportion of embryos affected
ggplot(combined_data, aes(x = chrom, y = proportion, fill = parent_origin)) +
  geom_bar(stat = "identity", position = "stack") + # Use 'identity' for precomputed values
  scale_fill_manual(
    values = c("maternal" = "#A020F0", "paternal" = "#008B8B")
  ) +
  labs(
    title = "Aneuploidy per chromosome by parent of origin",
    x = "Chromosome",
    y = "Proportion of Embryos Affected"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

# Disconnect
dev.off()



