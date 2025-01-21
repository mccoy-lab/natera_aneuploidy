# Crossovers per embryo, by parent of origin 

# =================
# author: Sara A. Carioscia, Biology Dept., Johns Hopkins University
# email: scarios1@jhu.edu
# last update: January 15, 2025
# aim: plot number of maternal and paternal crossovers per euploid embryo 
# =================

# load libraries
library(data.table)
library(dplyr)
library(ggplot2)

# Read in data 
crossovers <- fread("/scratch16/rmccoy22/abiddan1/natera_recomb/analysis/co_post_process/results/v30b_heuristic_90_nsib_qual.crossover_filt.deCode_haldorsson19.merged.meta.tsv.gz")

# Subset to only euploid embryos 
crossovers_euploid <- crossovers[crossovers$euploid == TRUE,]

# Group by embryo and count number of maternal and paternal crossovers each 
crossovers_euploid$identifier <- paste0(crossovers_euploid$mother, "_", crossovers_euploid$father, "_", crossovers_euploid$child)
crossovers_by_parent <- crossovers_euploid %>%
  group_by(identifier) %>%
  summarize(
    num_maternal = sum(crossover_sex == "maternal"),
    num_paternal = sum(crossover_sex == "paternal"),
    .groups = "drop"
  )

# Filter to keep only those within 3 sd of the mean number of crossovers 
crossovers_by_parent_filtered <- crossovers_by_parent %>%
  filter(
    abs(num_maternal - mean(num_maternal)) <= 3 * sd(num_maternal),
    abs(num_paternal - mean(num_paternal)) <= 3 * sd(num_paternal)
  )


# Create PDF
pdf("/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/characterize_aneuploidy/crossovers_by_parent.pdf", height = 8, width = 12)

# Plot number of maternal crossovers and number of paternal crossovers on the x axis 
ggplot() + 
  geom_bar(data = crossovers_by_parent_filtered, aes(x = num_maternal), fill = "#A020F0") + 
  geom_bar(data = crossovers_by_parent_filtered, aes(x = num_paternal), fill = "#008B8B") + 
  xlab("Number of Autosomal Crossovers") + 
  ylab("Number of Euploid Embryos") + 
  theme_minimal()

# Disconnect
dev.off()
  