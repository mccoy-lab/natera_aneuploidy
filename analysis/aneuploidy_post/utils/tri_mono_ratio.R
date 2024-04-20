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
genes_by_chr <- gencode[gencode$type == "gene",] %>%
  # keep just genes and remove pseudogenes
  filter(!grepl("pseudogene", gene_type, ignore.case = TRUE)) %>% 
  group_by(seqid) %>%
  count() %>% 
  rename(num_genes = n) %>% 
  rename(chrom = seqid)

# get size of each chromosome 
# https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&chromInfoPage=&pix=1752
chrom_sizes <- fread("hg38.chrom.sizes.txt") 
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
  geom_point() +
  geom_text(size = 3, vjust = -0.5, hjust = 1) +  
  labs(
    x = "Number of Genes",
    y = "Proportion of Monosomies",
    title = "Monosomies by Number of Genes",
  ) +
  theme_minimal()

# ratio of monosomies vs. number of genes (corrected for chr length)
ggplot(ploidy_chr, aes(x = num_genes / bp, 
                         y = monosomies / aneuploidies, label = chrom)) +
  geom_point() +
  geom_text(size = 3, vjust = -0.5, hjust = 1) +  
  labs(
    x = "Number of Genes / Size of Chr",
    y = expression(paste(n[monosomy], "/", n[aneuploidy])),
    #title = "Enrichment for Monosomies",
  ) +
  theme_minimal()

# ratio of monosomies vs. size of chromosome 
ggplot(ploidy_chr, aes(x = bp, 
                       y = monosomies / aneuploidies, label = chrom)) +
  geom_point() +
  geom_text(size = 3, vjust = -0.5, hjust = 1) +  
  labs(
    x = "Size of Chromosome (bp)",
    y = expression(paste(n[monosomy], "/", n[aneuploidy])),
    title = "Proportion of Monosomies by Chromosome Size",
  ) +
  theme_minimal()

