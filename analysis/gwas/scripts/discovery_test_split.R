library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
#library(ggpubr)

# Usage: 
# ./discovery_test_split.R \
# "/data/rmccoy22/natera_spectrum/data/summary_metadata/spectrum_metadata_merged.csv" \
# "/data/rmccoy22/natera_spectrum/genotypes/opticall_parents_100423/genotypes/opticall_concat_total.norm.b38.fam" \
# "/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/king_result.king.cutoff.out.id" \ 
# "/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/results/discover_test_split_mother.txt" \
# "/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/results/discover_test_split_father.txt" 

# get command line arguments
args <- commandArgs(trailingOnly = TRUE)
metadata <- args[1]
fam <- args[2]
king_related_arrays <- args[3]
output_maternal <- args[4]
output_paternal <- args[5]

# read files from args
metadata <- fread(metadata)
fam <- fread(fam) 
king_related_arrays <- fread(king_related_arrays)

# Aggregate Across Families
# remove duplicate individuals  
metadata <- metadata %>%
  distinct(array, .keep_all = TRUE)

# apply same array ID to all individuals (partner, child) associated with each 
# mother, across different casefiles 
metadata_merged_array <- metadata %>%
  mutate(array_id_merged = ifelse(family_position == "mother", 
                                  paste0(array, "_merged"), NA_character_)) %>%
  group_by(casefile_id) %>%
  fill(array_id_merged) %>%
  ungroup() 


# create dataframe with 1) number of embryos total per mother and 
# 2) weighted age of mother by embryos 
weighted_ages <- metadata_merged_array %>%
  filter(family_position == "child") %>%
  group_by(array_id_merged) %>%
  summarise(weighted_age = sum(patient_age) / n(),
            weighted_partner_age = sum(partner_age) / n(),
            child_count = n()) %>% 
  as.data.frame()

# add the embryo count and weighted age columns to the metadata table
metadata_merged_array_ages <- merge(weighted_ages, metadata_merged_array, 
                                    by = "array_id_merged") %>%
  as.data.table()


# Keep only parents, and only those that were successfully genotyped 
# identify which parents were genotyped (present in the genotyped .fam file)
metadata_merged_array_ages[, is_genotyped := array %in% fam$V1] %>%
  setorder(array)
# keep only individuals who are genotyped 
metadata_merged_array_ages <- metadata_merged_array_ages[is_genotyped == TRUE]


# Remove individuals that were related 
# create table of full metadata for any related individuals 
related_metadata <- merge(metadata_merged_array_ages, king_related_arrays, 
                          by.x = "array", by.y = "IID")
# add indicator in metadata for whether individual is related 
metadata_merged_array_ages[, related_samples_to_drop := 
                             casefile_id %in% related_metadata$casefile_id]
# remove families with related individuals from metadata 
metadata_merged_array_ages <- 
  metadata_merged_array_ages[related_samples_to_drop == FALSE]

# Filter parent age range 
# keep only parents with ages between 18-90 or NA 
metadata_merged_array_ages <- metadata_merged_array_ages[
  ((patient_age > 18 & patient_age < 90) | is.na(patient_age)) & 
    ((partner_age > 18 & partner_age < 90) | is.na(partner_age)),]


# Assign egg and sperm donor ages 
# get average egg donor age from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7530253/
# get average sperm donor age from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9118971/#:~:text=Donors%20were%20aged%2027%20years,aged%2030%20years%20and%20younger.
metadata_merged_array_ages[egg_donor == "yes", weighted_age := 25]
metadata_merged_array_ages[sperm_donor == "yes", weighted_partner_age := 27]


# Separate into discovery and test sets while maintaining split on 
# key covariates (maternal age, embryo count)
# get just mothers to split into test and discovery set 
metadata_merged_array_ages_mothers <- metadata_merged_array_ages[
  metadata_merged_array_ages$family_position == "mother",]

# Adapted from https://gettinggeneticsdone.blogspot.com/2011/03/splitting-dataset-revisited-keeping.html
# splitdf splits a data frame into a discovery and a test set
splitdf <- function(dataframe, trainfrac, seed = NULL) {
  if (trainfrac <= 0 | trainfrac >= 1) 
    stop("Training fraction must be between 0 and 1, not inclusive")
  if (!is.null(seed)) set.seed(seed)
  index <- 1:nrow(dataframe)
  trainindex <- sample(index, trunc(length(index)/(1/trainfrac)))
  trainset <- dataframe[trainindex, ]
  testset <- dataframe[-trainindex, ]
  list(trainset = trainset, testset = testset)
}

# Use splitdf to split the dataframe, keeping discovery and test sets equivalent
# across the input character covariates 
splitdf.randomize <- function(dataframe, min_p, trainfrac, 
                              ttestcolnames = c("cols","to","test"), ...) {
  d <- dataframe
  if (!all(ttestcolnames %in% names(d))) 
    stop(paste(ttestcolnames,"not in dataframe"))
  ps <- NULL
  while (is.null(ps) | any(ps < min_p)) {
    sets <- splitdf(d, trainfrac)
    trainset <- sets$trainset
    testset <- sets$testset
    #ttestcols <- which(names(d) %in% ttestcolnames)
    ps <- NULL
    p <- t.test(trainset[ ,2], testset[ ,2])$p.value
    ps = c(ps, p)
    p <- t.test(trainset[ ,3], testset[ ,3])$p.value
    ps = c(ps, p)
    print(paste(ttestcolnames," t-test p-value =",ps))
    cat("\n")
  }
  list(trainset=trainset,testset=testset)
}

# get the column names containing the covariates of interest (col #2 and 3)
cols <- c("weighted_age", "child_count")

# split dataset, keeping even distribution of those covariates
set.seed(5)
evensplit <- splitdf.randomize(metadata_merged_array_ages_mothers, min_p=0.05, trainfrac=0.85, cols)

# add column to metadata of mothers noting whether they're discovery 
metadata_merged_array_ages_mothers[, is_discovery := array %in% evensplit$trainset$array]


# write mothers to file 
fwrite(metadata_merged_array_ages_mothers[,c("array", "family_position", "is_discovery")], 
       file = output_maternal, 
       sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# propagate those assignments to the fathers
case_assign <- metadata_merged_array_ages_mothers[, c("casefile_id", "is_discovery")]
metadata_fathers <- metadata_merged_array_ages %>%
  .[family_position == "father"] %>%
  merge(., case_assign, by = "casefile_id")
# remove duplicates
metadata_fathers <- metadata_fathers[!duplicated(metadata_fathers$array)]

# write fathers to file 
fwrite(metadata_fathers[, c("array", "family_position", "is_discovery")], 
       file = output_paternal, 
       sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


# # plot cumulative dist of the covariates 
# p1 <- ggplot(data = metadata_merged_array_ages_mothers, aes(x = weighted_age, color = is_discovery)) +
#   geom_density() +
#   theme_bw() +
#   theme(axis.line = element_line(), panel.grid = element_blank(), panel.border = element_blank()) +
#   xlab("Patient Age") +
#   ylab("Density") + 
#   scale_color_manual(labels = c("Validation", "Discovery"), values = c("blue", "red")) + 
#   guides(color=guide_legend("Data Split Assignment"))

# p2 <- ggplot(data = metadata_merged_array_ages_mothers, aes(x = partner_age, color = is_discovery)) +
#   geom_density() +
#   #xlim(20, 83) + 
#   theme_bw() +
#   theme(axis.line = element_line(), panel.grid = element_blank(), panel.border = element_blank()) +
#   xlab("Partner Age") +
#   ylab("Density") + 
#   scale_color_manual(labels = c("Validation", "Discovery"), values = c("blue", "red")) + 
#   guides(color=guide_legend("Data Split Assignment"))

# p3 <- ggplot(data = metadata_merged_array_ages_mothers, aes(x = child_count, color = is_discovery)) +
#   geom_density() +
#   theme_bw() +
#   theme(axis.line = element_line(), panel.grid = element_blank(), panel.border = element_blank()) +
#   xlab("Number of Embryos") +
#   ylab("Density") + 
#   scale_color_manual(labels = c("Validation", "Discovery"), values = c("blue", "red")) + 
#   guides(color=guide_legend("Data Split Assignment"))


# # plot all three as one figure 
# pdf(plot_fpn)
# ggpubr::ggarrange(p1, p2, p3, # list of plots
#                   labels = "AUTO", # labels
#                   common.legend = T, # COMMON LEGEND
#                   legend = "bottom", # legend position
#                   align = "hv", # Align them both, horizontal and vertical
#                   nrow = 1)  # number of rows
# dev.off()

