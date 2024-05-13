# Number of affected chr per embryo 
# Aneuploidy incidence (by type) per chr

# load libraries
library(data.table)
library(dplyr)
library(ggplot2) 

# Usage: /aneuploidy_post/utils/embryo_ploidy.R \ 
# "/data/rmccoy22/natera_spectrum/karyohmm_outputs/compiled_output/natera_embryos.karyohmm_v30a.bph_sph_trisomy.full_annotation.031624.tsv.gz" \
# "mother"
# "/data/rmccoy22/natera_spectrum/data/summary_metadata/spectrum_metadata_merged.csv" 
# 15 # threshold number of chr to count as whole-genome 
# "/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/aneuploidy_post/utils/test.pdf"

# get command line arguments 
args <- commandArgs(trailingOnly = TRUE)
# ploidy calls 
ploidy_calls <- args[1]
# parent to get gain/loss by 
parent <- args[2]
# metadata to remove day3 embryos 
metadata <- args[3]
# threshold number of affected chr for whole genome ploidy 
ploidy_threshold <- args[4]
# out plot for whole genome gain/loss counts 
out_ploidy_by_chr <- args[5]


# read in data 
ploidy_calls <- fread(ploidy_calls)
metadata <- fread(metadata)
max_meiotic = 3


# get functions `filter_data` and `day5_only`
setwd(".")
source("../../gwas/scripts/phenotypes/helper_functions/phenotyping_helper_functions.R")

# filter data 
ploidy_calls_filtered <- filter_data(ploidy_calls, parent, 
                                     bayes_factor_cutoff = 2, 
                                     nullisomy_threshold = 5, min_prob = 0.9) 
ploidy_calls_filtered <- day5_only(ploidy_calls_filtered, metadata)


# get count of affected chr per embryo
affected_chr_per_embryo <- ploidy_calls_filtered %>%
  group_by(child) %>%
  summarise(num_affected_chromosomes = sum(bf_max_cat != 2)) %>%
  ungroup()

# count number of embryos with single-chr aneuploidies 
num_meiotic <- length(which(affected_chr_per_embryo$num_affected_chromosomes > 0 
                            & affected_chr_per_embryo$num_affected_chromosomes <= max_meiotic))
      
# count number of embryos with whole-genome gain/loss
num_ploidy <- length(which(affected_chr_per_embryo$num_affected_chromosomes 
                           >= ploidy_threshold))

# count number of embryos total
num_embryos <- length(unique(affected_chr_per_embryo$child))

# plot count of affected chromosomes per embryo
ggplot(affected_chr_per_embryo, aes(x = num_affected_chromosomes)) +
  geom_bar() +
  labs(title = "Number of Affected Chromosomes per Embryo",
       x = "Number of Affected Chromosomes",
       y = "Count of Embryos") + 
  theme_minimal()



# explore error type by chromosome 
# create new column with just numeric value of chromosome to improve plotting
# remove the prefix "chr" from chromosome names
ploidy_calls_filtered$chr <- gsub("chr", "", ploidy_calls_filtered$chrom)
# convert chromosome names to factors and reorder them numerically
ploidy_calls_filtered$chr <- factor(ploidy_calls_filtered$chr, levels = as.character(1:22))

# remove embryos experiencing whole genome gain/loss 
if (remove_whole_genome == TRUE) {
  affected_chr_per_embryo <- ploidy_calls_filtered %>%
    group_by(child) %>%
    summarise(num_affected_chromosomes = sum(bf_max_cat != 2)) %>%
    ungroup()
  whole_genome <- affected_chr_per_embryo[affected_chr_per_embryo$num_affected_chromosomes < ploidy_threshold,]
  ploidy_calls_filtered <- ploidy_calls_filtered[ploidy_calls_filtered$child %in% whole_genome$child,]
}

# statistics for different amounts of aneuploidy on each chromosome 

# test 1: does incidence of aneuploidy differ across chromosomes? yes, p-value < 2.2e-16 
aneu_vs_eu <- ploidy_calls_filtered
aneu_vs_eu$bf_max_cat <- factor(aneu_vs_eu$bf_max_cat)
aneu_vs_eu$ploidy <- ifelse(aneu_vs_eu$bf_max_cat %in% c("3m", "3p", "1m", "1p", "0"), "aneuploid", "euploid")
aneu_vs_eu_table <- table(aneu_vs_eu$chrom, aneu_vs_eu$ploidy)
results_aneu_vs_eu <- chisq.test(aneu_vs_eu_table)

# test 2: does ratio of maternal to paternal aneuploidy differ across chromosomes? yes, p-value < 2.2e-16
mat_vs_pat <- ploidy_calls_filtered[ploidy_calls_filtered$bf_max_cat %in% c("3m", "3p", "1m", "1p"),]
mat_vs_pat$bf_max_cat <- factor(mat_vs_pat$bf_max_cat)
mat_vs_pat$parent <- ifelse(mat_vs_pat$bf_max_cat %in% c("3m", "1p"), "maternal", "paternal") 
mat_vs_pat_table <- table(mat_vs_pat$chrom, mat_vs_pat$parent)
result_mat_vs_pat <- chisq.test(mat_vs_pat_table)

# test 3: does ratio of monosomy to trisomy differ across chromosomes? yes, p-value < 2.2e-16
mono_tri <- ploidy_calls_filtered[ploidy_calls_filtered$bf_max_cat %in% c("3m", "3p", "1m", "1p"),]
mono_tri$bf_max_cat <- factor(mono_tri$bf_max_cat)
mono_tri$cn <- ifelse(mono_tri$bf_max_cat %in% c("3m", "3p"), "trisomy", "monosomy") 
mono_tri_table <- table(mono_tri$chrom, mono_tri$cn)
result_mono_tri <- chisq.test(mono_tri_table)

# test 4: does ratio of aneuploidy type differ between chromosomes? yes, p-value < 2.2e-16 
incidence <- ploidy_calls_filtered
incidence$bf_max_cat <- factor(incidence$bf_max_cat)
incidence_table <- table(incidence$chrom, incidence$bf_max_cat)
results_incidence_table <- chisq.test(incidence_table)

# test 5: does rate of paternal-origin aneuploidy differ between chromosomes? 
pat_all_chr <- ploidy_calls_filtered[ploidy_calls_filtered$bf_max_cat %in% c("3p", "1m"),]
mono_tri$bf_max_cat <- factor(mono_tri$bf_max_cat)

# reorder bf_max_cat for plotting 
plot_aneuploidy_order <- c("3m", "1p", "3p", "1m", "0", "2")
# Reorder bf_max_cat factor levels
ploidy_calls_filtered_sorted <- ploidy_calls_filtered
ploidy_calls_filtered_sorted$bf_max_cat <- factor(ploidy_calls_filtered$bf_max_cat, 
                                                  levels = plot_aneuploidy_order)

# set colors for each aneuploidy type 
colors <- c("3m" = "#E69F00", "1p" = "#CC79A7", "3p" = "#0072B2", "1m" = "#56B4E9", "0" = "#009E73", "2" = "#D55E00")


# plot aneuploidies from both parents
ggplot(ploidy_calls_filtered_sorted[ploidy_calls_filtered_sorted$bf_max_cat != "2",], 
       aes(x = chr, fill = factor(bf_max_cat))) + 
  geom_bar(stat = "count") +
  labs(title = "Ploidy Status by Chromosome",
       x = "Chromosome",
       y = "Count") +
  theme_minimal() + 
  scale_fill_manual(name = "Aneuploidy Type", 
                    values = colors,
                    labels = c("Gain of Maternal", "Loss of Maternal", 
                                 "Gain of Paternal", "Loss of Paternal", 
                                 "Nullisomy")) 

# plot paternal gain/loss 
ggplot(ploidy_calls_filtered_sorted[ploidy_calls_filtered_sorted$bf_max_cat %in% c("1m", "3p"),], 
       aes(x = chr, fill = factor(bf_max_cat))) + 
  geom_bar(stat = "count") +
  labs(title = "Paternal Aneuploidies by Chromosome",
       x = "Chromosome",
       y = "Count") +
  theme_minimal() + 
  scale_fill_manual(name = "Aneuploidy Type", 
                    values = colors,
                    labels = c("Loss of Paternal", "Gain of Paternal"))

# plot maternal gain/loss 
ggplot(ploidy_calls_filtered_sorted[ploidy_calls_filtered_sorted$bf_max_cat %in% c("1p", "3m"),], 
       aes(x = chr, fill = factor(bf_max_cat))) + 
  geom_bar(stat = "count") +
  labs(title = "Maternal Aneuploidies by Chromosome",
       x = "Chromosome",
       y = "Count") +
  theme_minimal() + 
  scale_fill_manual(name = "Aneuploidy Type",
                    values = colors,
                    labels = c("Loss of Maternal", "Gain of Maternal"))

# plot aneuploidies as percentage 
ggplot(ploidy_calls_filtered_sorted[ploidy_calls_filtered_sorted$bf_max_cat != "2",], 
       aes(x = chr, fill = factor(bf_max_cat))) + 
  geom_bar(position = "fill") +
  labs(title = "Proportion of Aneuploidy Types",
       x = "Chromosome",
       y = "Count") +
  theme_minimal() + 
  scale_fill_manual(name = "Aneuploidy Type",
                    values = colors,
                    labels = c("Gain of Maternal", "Loss of Maternal", 
                                 "Gain of Paternal", "Loss of Paternal", 
                                 "Nullisomy"))

# plot all outcomes (including disomy) as percentage
ggplot(ploidy_calls_filtered_sorted, 
       aes(x = chr, fill = factor(bf_max_cat))) + 
  geom_bar(position = "fill") +
  theme_minimal() + 
  scale_fill_manual(name = "Aneuploidy Type",
                    values = colors, 
                    labels = c("Gain of Maternal", "Loss of Maternal", 
                               "Gain of Paternal", "Loss of Paternal", 
                               "Nullisomy", "Disomy (Euploid)"))


# plot all aneuploidies with no colors 
ggplot(ploidy_calls_filtered_sorted[!(ploidy_calls_filtered_sorted$bf_max_cat %in% c("2", "0")),], 
       aes(x = chr)) + 
  geom_bar(stat = "count") +
  theme_minimal() + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) 


ggplot(ploidy_calls_filtered_sorted[!(ploidy_calls_filtered_sorted$bf_max_cat %in% c("2", "0")),], 
       aes(x = chr, fill = bf_max_cat)) + 
  geom_bar(stat = "count") +
  theme_minimal() + 
  scale_fill_manual(values = colors_parents) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none") 


colors_parents <- c("mother" = "#CC79A7", "father" = "#0072B2")


# classify aneuploidy as maternal or paternal 
ploidy_calls_filtered_sorted$parent <- ifelse(
  ploidy_calls_filtered_sorted$bf_max_cat %in% c("3m", "1p"), "mother",
  ifelse(
    ploidy_calls_filtered_sorted$bf_max_cat %in% c("3p", "1m"), "father",
    "neither"
  )
)

# plot aneuploidies as proportion, with no color 
ggplot(ploidy_calls_filtered_sorted[!(ploidy_calls_filtered_sorted$bf_max_cat %in% c("2", "0")),], 
       aes(x = chr, y = ..count../73910)) + 
  geom_bar(stat = "count", position = position_stack(reverse = TRUE)) +
  theme_minimal() + 
  scale_fill_manual(values = colors_parents) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        legend.position = "none")

# plot aneuploidies, color by parent of origin and show as proportion 
ggplot(ploidy_calls_filtered_sorted[!(ploidy_calls_filtered_sorted$bf_max_cat %in% c("2", "0")),], 
       aes(x = chr, y = ..count../73910, fill = factor(parent))) + 
  geom_bar(stat = "count", position = position_stack(reverse = TRUE)) +
  theme_minimal() + 
  scale_fill_manual(values = colors_parents) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none") 

# chi sq for paternal error across chromosomes 
chi_sq_test <- chisq.test(table(ploidy_calls_filtered_sorted$chr, ploidy_calls_filtered_sorted$father))



# format to write plot to pdf 
pdf(out_ploidy_by_chr)
ggsave(out_ploidy_by_chr, ploidy_by_chr)
