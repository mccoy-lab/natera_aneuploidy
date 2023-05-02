library(data.table)
library(ggplot2)
library(ggpubr)

library(dplyr)

# input files necessary 
metadata <- fread("/data/rmccoy22/natera_spectrum/data/summary_metadata/spectrum_metadata_merged.csv")
fam <- fread("/data/rmccoy22/natera_spectrum/genotypes/opticall_parents_031423/genotypes/opticall_concat_total.norm.b38.fam")
king_related_arrays <- fread("/scratch16/rmccoy22/scarios1/natera_spectrum/gwas/background/king_outputs/kingunrelated_toberemoved.txt", header = FALSE)

# output filepaths 
plot_fp <- "/scratch16/rmccoy22/scarios1/natera_spectrum/gwas/background/discovery_validation_split_covariates.pdf"
output_maternal <- "/scratch16/rmccoy22/scarios1/natera_spectrum/gwas/background/discover_validate_split_f.txt"
output_paternal <- "/scratch16/rmccoy22/scarios1/natera_spectrum/gwas/background/discover_validate_split_m.txt"

# remove families affected by donors 
metadata <- replace(metadata, metadata=='', NA)
donors <- rbind(metadata[!is.na(egg_donor)], metadata[!is.na(sperm_donor)])
metadata[, donors_to_drop := casefile_id %in% donors$casefile_id]
metadata <- metadata[donors_to_drop == FALSE]

# create new column that tags every individual affiliated with each mother even if in different caseIDs
metadata_merged_array <- metadata %>%
  mutate(array_id_merged = ifelse(family_position == "mother", paste0(array, "_merged"), NA_character_)) %>%
  group_by(casefile_id) %>%
  fill(array_id_merged) %>%
  ungroup() 
# create dataframe that calculates the weighted age (based on age at each embryo) and number of embryos 
weighted_ages <- metadata_merged_array %>%
  filter(family_position == "child") %>%
  group_by(array_id_merged) %>%
  summarise(weighted_age = sum(patient_age) / n(),
            child_count = n()) %>% 
  as.data.frame()

# add the weighted age and embryo count columns to the main metadata table
metadata_merged_array_ages <- merge(weighted_ages, metadata_merged_array, by = "array_id_merged")
metadata_merged_array_ages <- metadata_merged_array_ages %>% as.data.table()


# get individuals who were genotyped successfully and in the bed output .fam file 
metadata_merged_array_ages[, is_genotyped := array %in% fam$V1] %>%
  setorder(array)
# keep only individuals who are genotyped 
metadata_merged_array_ages <- metadata_merged_array_ages[is_genotyped == TRUE]

# remove all families affected by a related individual 
related_metadata <- merge(metadata_merged_array_ages, king_related_arrays, by.x = "array", by.y = "V1")
# keep only individuals in families that did not have a related individual 
metadata_merged_array_ages[, related_samples_to_drop := casefile_id %in% related_metadata$casefile_id]
metadata_merged_array_ages <- metadata_merged_array_ages[related_samples_to_drop == FALSE]

# filter parent age range - keep individuals for which the age is NA 
metadata_merged_array_ages <- metadata_merged_array_ages[((patient_age > 18 & patient_age < 90) | is.na(patient_age)) & ((partner_age > 18 & partner_age < 90) | is.na(partner_age)),]

# get just mothers to split 
metadata_merged_array_ages_mothers <- metadata_merged_array_ages[metadata_merged_array_ages$family_position == "mother",]
# keep each mother only once 
metadata_merged_array_ages_mothers <- metadata_merged_array_ages_mothers[!duplicated(metadata_merged_array_ages_mothers$array)]


## Try this from https://gettinggeneticsdone.blogspot.com/2011/03/splitting-dataset-revisited-keeping.html
#splitdf splits a data frame into a training and testing set. returns a list of two data frames: trainset and testset. you can optionally apply a random seed.
splitdf <- function(dataframe, trainfrac, seed=NULL) {
  if (trainfrac<=0 | trainfrac>=1) stop("Training fraction must be between 0 and 1, not inclusive")
  if (!is.null(seed)) set.seed(seed)
  index <- 1:nrow(dataframe)
  trainindex <- sample(index, trunc(length(index)/(1/trainfrac)))
  trainset <- dataframe[trainindex, ]
  testset <- dataframe[-trainindex, ]
  list(trainset=trainset,testset=testset)
}

#this function utilizes the function above. you give it a data frame you want to randomize, and a character vector with column names you want to be sure are 
#equally distributed among the two different sets. these columns must be continuous variables. 
splitdf.randomize <- function(dataframe, min_p, trainfrac, ttestcolnames=c("cols","to","test"), ...) {
  d <- dataframe
  if (!all(ttestcolnames %in% names(d))) stop(paste(ttestcolnames,"not in dataframe"))
  ps <- NULL
  while (is.null(ps) | any(ps<min_p)) {
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

# get the column names of the covariates I'm testing (col #2 and 3)
cols <- c("weighted_age", "child_count")

# split up the dataset, keeping even distribution of those covariates
set.seed(5)
evensplit <- splitdf.randomize(metadata_merged_array_ages_mothers, min_p=0.05, trainfrac=0.85, cols)
# weighted_age  t-test p-value = 0.825539076325965" "child_count  t-test p-value = 0.649560315142923

# add column to metadata of mothers noting whether they're discovery 
metadata_merged_array_ages_mothers[, is_discovery := array %in% evensplit$trainset$array]


# plot cumulative dist of the covariates 
p1 <- ggplot(data = metadata_merged_array_ages_mothers, aes(x = weighted_age, color = is_discovery)) +
  geom_density() +
  theme_bw() +
  theme(axis.line = element_line(), panel.grid = element_blank(), panel.border = element_blank()) +
  xlab("Patient Age") +
  ylab("Density") + 
  scale_color_manual(labels = c("Validation", "Discovery"), values = c("blue", "red")) + 
  guides(color=guide_legend("Data Split Assignment"))


p2 <- ggplot(data = metadata_merged_array_ages_mothers, aes(x = partner_age, color = is_discovery)) +
  geom_density() +
  #xlim(20, 83) + 
  theme_bw() +
  theme(axis.line = element_line(), panel.grid = element_blank(), panel.border = element_blank()) +
  xlab("Partner Age") +
  ylab("Density") + 
  scale_color_manual(labels = c("Validation", "Discovery"), values = c("blue", "red")) + 
  guides(color=guide_legend("Data Split Assignment"))


p3 <- ggplot(data = metadata_merged_array_ages_mothers, aes(x = child_count, color = is_discovery)) +
  geom_density() +
  theme_bw() +
  theme(axis.line = element_line(), panel.grid = element_blank(), panel.border = element_blank()) +
  xlab("Number of Embryos") +
  ylab("Density") + 
  scale_color_manual(labels = c("Validation", "Discovery"), values = c("blue", "red")) + 
  guides(color=guide_legend("Data Split Assignment"))


# Plot all three variables as one figure 
pdf(plot_fp)
ggpubr::ggarrange(p1, p2, p3, # list of plots
                  labels = "AUTO", # labels
                  common.legend = T, # COMMON LEGEND
                  legend = "bottom", # legend position
                  align = "hv", # Align them both, horizontal and vertical
                  nrow = 1)  # number of rows
dev.off()

# write mothers to file 
fwrite(metadata_merged_array_ages_mothers[,c("array", "family_position", "is_discovery")], 
       file = output_maternal, 
       sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

case_assign <- metadata_merged_array_ages_mothers[, c("casefile_id", "is_discovery")]

# then propogate those assignments to the fathers
metadata_fathers <- metadata_merged_array_ages %>%
  .[family_position == "father"] %>%
  merge(., case_assign, by = "casefile_id")
# remove duplicates
metadata_fathers <- metadata_fathers[!duplicated(metadata_fathers$array)]
# write fathers to file 
fwrite(metadata_fathers[, c("array", "family_position", "is_discovery")], 
       file = output_paternal, 
       sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)




# Check if a random split would yield even distribution of covariates 
#discovery <- sample(unique(metadata_mothers_child_count$array), round(length(unique(metadata_mothers_child_count$array)) * 0.7))
#metadata_mothers_child_count[, is_discovery := array %in% discovery] %>%
#  setorder(., -is_discovery)

#t.test(metadata_mothers_child_count$patient_age ~ metadata_mothers_child_count$is_discovery) # p-value = 0.8304
#t.test(metadata_mothers_child_count$count ~ metadata_mothers_child_count$is_discovery) # p-value = 0.2757


