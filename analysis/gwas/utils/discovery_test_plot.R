# Plot cumulative distribution of the covariates across the discovery and 
# test sets 

library(ggplot2)
library(ggpubr)

# Usage: 
# "./discovery_test_plot.R" \ 
# "/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/results/discovery_test_split_mother.txt" \
# "/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/results/discovery_test_covariates.pdf"

# get command line arguments
args <- commandArgs(trailingOnly = TRUE)
discovery_test_mothers <- args[1]
out_pdf <- args[2]

# read in discovery test data 
metadata_merged_array_ages_mothers <- fread(discovery_test_mothers)

# plot weighted maternal age 
p1 <- ggplot(data = metadata_merged_array_ages_mothers, 
             aes(x = weighted_age, color = is_discovery)) +
  geom_density() +
  theme_bw() +
  theme(axis.line = element_line(), panel.grid = element_blank(), 
        panel.border = element_blank()) +
  xlab("Patient Age") +
  ylab("Density") + 
  scale_color_manual(labels = c("Test", "Discovery"), 
                     values = c("blue", "red")) + 
  guides(color = guide_legend("Data Split Assignment"))

# plot weighted maternal age 
p2 <- ggplot(data = metadata_merged_array_ages_mothers, 
             aes(x = weighted_partner_age, color = is_discovery)) +
  geom_density() +
  theme_bw() +
  theme(axis.line = element_line(), panel.grid = element_blank(), 
        panel.border = element_blank()) +
  xlab("Partner Age") +
  ylab("Density") + 
  scale_color_manual(labels = c("Test", "Discovery"), 
                     values = c("blue", "red")) + 
  guides(color = guide_legend("Data Split Assignment"))


# plot number of embryos per mother
p3 <- ggplot(data = metadata_merged_array_ages_mothers, 
             aes(x = child_count, color = is_discovery)) +
  geom_density() +
  theme_bw() +
  theme(axis.line = element_line(), panel.grid = element_blank(), 
        panel.border = element_blank()) +
  xlab("Number of Embryos") +
  ylab("Density") + 
  scale_color_manual(labels = c("Test", "Discovery"), 
                     values = c("blue", "red")) + 
  guides(color = guide_legend("Data Split Assignment"))

# plot all three as one figure 
pdf(covariate_figs)
ggpubr::ggarrange(p1, p2, p3, 
                  labels = "AUTO", 
                  common.legend = T, 
                  legend = "bottom", 
                  align = "hv", 
                  nrow = 1)
dev.off()
