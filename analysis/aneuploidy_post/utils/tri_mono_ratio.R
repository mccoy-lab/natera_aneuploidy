# Comparison of trisomy:monosomy ratio by number of genes per chrom

# load libraries
library(data.table)
library(dplyr)
library(rtracklayer)
library(ggplot2) 

# Usage: /aneuploidy_post/utils/tri_mono_ratio.R \ 
# "/data/rmccoy22/natera_spectrum/karyohmm_outputs/compiled_output/natera_embryos.karyohmm_v30a.bph_sph_trisomy.full_annotation.031624.tsv.gz" \
# "/data/rmccoy22/natera_spectrum/data/summary_metadata/spectrum_metadata_merged.csv"  
# mother 

# get command line arguments 
args <- commandArgs(trailingOnly = TRUE)
# ploidy calls 
ploidy_calls <- args[1]
# metadata to filter to day 5s 
metadata <- args[2]
# parent for filtering 
parent <- args[3]

# get functions `filter_data` and `day5_only`
setwd(".")
source("../../gwas/scripts/phenotypes/helper_functions/phenotyping_helper_functions.R")


# NATERA DATA 

ploidy_calls <- fread(ploidy_calls)
metadata <- fread(metadata)

# filter data 
ploidy_calls_filtered <- filter_data(ploidy_calls, parent, 
                                     bayes_factor_cutoff = 2, 
                                     nullisomy_threshold = 5, min_prob = 0.9) 
ploidy_calls_filtered <- day5_only(ploidy_calls_filtered, metadata)

# make table with chromosome, trisomies, monosomies, ratio
copy_numbers <- ploidy_calls_filtered %>% 
  group_by(chrom) %>%
  summarise(trisomies = sum(bf_max_cat %in% c("3m", "3p")), 
            monosomies = sum(bf_max_cat %in% c("1m", "1p")), 
            aneuploidies = sum(bf_max_cat != "2")) 
# make column for ratio 
# copy_numbers$tri_mono_ratio = copy_numbers$trisomies / copy_numbers$monosomies


# CHROMOSOME DATA 

# get number of genes per chromosome 
gencode <- readGFF("/data/rmccoy22/resources/GENCODE/v38/gencode.v38.annotation.sorted.gtf")
genes_by_chr <- gencode[gencode$type == "gene" & gencode$gene_type == "protein_coding",] %>%
  # keep just genes and remove pseudogenes
  filter(!grepl("pseudogene", gene_type, ignore.case = TRUE)) %>% 
  group_by(seqid) %>%
  count() 
colnames(genes_by_chr) <- c("chrom", "num_genes")


# get size of each chromosome 
# https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&chromInfoPage=&pix=1752
chrom_sizes <- fread("/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/aneuploidy_post/utils/hg38.chrom.sizes.txt") 
colnames(chrom_sizes) <- c("chrom", "bp")


# COMBINE DATA 

# add number of genes to aneuploidy table
ploidy_chr <- left_join(copy_numbers, genes_by_chr, by = "chrom")
# add size of each chrom to aneuploidy table 
ploidy_chr <- left_join(ploidy_chr, chrom_sizes, by = "chrom")

# convert chrom to factor in numeric order
chrom_order <- c(paste0("chr", 1:22), "chrM", "chrX", "chrY")
ploidy_chr$chrom <- factor(ploidy_chr$chrom, levels = chrom_order)


# PLOTS 

# number of genes per chromosome 
ggplot(ploidy_chr, aes(x = chrom, y = num_genes)) + 
  geom_point() + 
  labs(
    x = "Chromosome",
    y = "Number of Genes",
    title = "Genes per Chromosome",
  ) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# number of genes per chromosome, corrected for length of chromosome 
ggplot(ploidy_chr, aes(x = chrom, y = num_genes/bp)) + 
  geom_point() + 
  labs(
    x = "Chromosome",
    y = "Number of Genes / Size of Chromosome",
    title = "Genes per Chromosome",
  ) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# proportion of aneuploidy vs. rate of monosomy
ggplot(ploidy_chr, aes(x = aneuploidies / nrow(ploidy_calls_filtered)/22, 
                         y = monosomies / aneuploidies, label = chrom)) +
  geom_point() +
  geom_text(size = 3, vjust = -0.5, hjust = 1) +  
  labs(
    x = "Proportion of Aneuploid Chromosomes",
    y = "Proportion of Monosomies",
    #title = "Enrichment for Monosomies",
  ) +
  theme_minimal()

# ratio of monosomies vs. number of genes (uncorrected for chr length)
ggplot(ploidy_chr, aes(x = num_genes, 
                       y = monosomies / aneuploidies, label = chrom)) +
  geom_point(size = 3) +
  geom_text(size = 4, vjust = -1, hjust = 1) +
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 13), 
        axis.text.y = element_text(size = 13), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) + 
  ylim(0.1, 0.5) + 
  xlim(200, 2100)

# ratio of monosomies vs. number of genes (corrected for chr length)
ggplot(ploidy_chr, aes(x = num_genes / (bp/1000), 
                         y = monosomies / aneuploidies, label = chrom)) +
  geom_point(size = 3) +
  geom_text(size = 4, vjust = -1, hjust = 1) +
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 13), 
        axis.text.y = element_text(size = 13), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) + 
  ylim(0.1, 0.5) + 
  xlim(0.002, 0.0255)

# ratio of monosomies vs. size of chromosome 
ggplot(ploidy_chr, aes(x = bp / 1000, 
                       y = monosomies / aneuploidies, label = chrom)) +
  geom_point(size = 3) + 
  geom_text(size = 4, vjust = -1, hjust = 1) + 
  labs(
    x = "Size of Chromosome (kbp)",
    y = expression(paste(n[monosomy], "/", n[monosomy] + n[trisomy])),
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 13), 
        axis.text.y = element_text(size = 13), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) + 
  ylim(0.1, 0.5) + 
  xlim(45000, 260000)


# reorganize data to get ratio of monosomy to trisomy by chromosome 
chroms <- ploidy_chr$chrom
monosomy_ratios <- (ploidy_chr$monosomies / (ploidy_chr$monosomies + ploidy_chr$trisomies))
monosomy_type <- rep("monosomy", times = 22)
trisomy_ratios <- (ploidy_chr$trisomies / (ploidy_chr$monosomies + ploidy_chr$trisomies))
trisomy_type <- rep("trisomy", times = 22)
mono <- data.frame(chroms, monosomy_ratios, monosomy_type)
colnames(mono) <- c("chromosomes", "ratio", "type")
tri <- data.frame(chroms, trisomy_ratios, trisomy_type)
colnames(tri) <- c("chromosomes", "ratio", "type")
mono_tri_ratios <- rbind(mono, tri)
mono_tri_ratios$chromosomes <- gsub("chr", "", mono_tri_ratios$chromosomes)
mono_tri_ratios$chromosomes <- factor(mono_tri_ratios$chromosomes, levels = as.character(1:22))
  

# color bar for plotting 
colors_mono_tri <- c("trisomy" = "#E69F00", "monosomy" = "#0072B2")
# adjust order for plotting 
mono_tri_ratios$type <- factor(mono_tri_ratios$type, levels = c("monosomy", "trisomy"))

# plot actual rates of monosomy and trisomy 
ggplot(mono_tri_ratios, aes(x = factor(chromosomes), y = ratio, fill = type)) +
  geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) +
  geom_hline(yintercept = 0.5, linetype = "dotted") + 
  scale_y_continuous(labels = scales::percent) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = colors_mono_tri)


# Simulate equal ratio of maternal monosomy to trisomy 
# Create a sample dataframe
chromosomes <- rep(1:22, each = 2)  
type <- rep(c("monosomy", "trisomy"), times = 22)  # Equal proportions of "monosomy" and "trisomy"
ratio <- c(rep(0.5, 44))  # 50% for each type
sim_mono_tri <- data.frame(chromosomes, types, ratio)

# plot simulation of equal rates of monosomy and trisomy 
ggplot(sim_mono_tri, aes(x = factor(chromosomes), y = ratio, fill = type)) +
  geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) +
  scale_y_continuous(labels = scales::percent) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = colors_mono_tri)


