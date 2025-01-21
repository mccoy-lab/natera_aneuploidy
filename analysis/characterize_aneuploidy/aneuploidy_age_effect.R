# Plot errors affecting maternal and paternal errors by maternal or paternal age

# =================
# author: Sara A. Carioscia, Biology Dept., Johns Hopkins University
# email: scarios1@jhu.edu
# last update: January 15, 2025
# aim: plot age association with aneuploidy, for each parent 
# =================

# load libraries
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2) 

# Read in data 
maternal_phenotype <- fread("/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/results/phenotypes/maternal_meiotic_aneuploidy_by_mother.csv")
paternal_phenotype <- fread("/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/results/phenotypes/paternal_meiotic_aneuploidy_by_father.csv")

# Function that reformats data for use in plotting
reformat_data <- function(phenotype, parent) {
  
  # Remove donors (age doesn't show the age of the donor) 
  phenotype_non_donor <- phenotype[which(phenotype$egg_donor != "yes" & phenotype$sperm_donor != "yes"),]
  
  # Set age column based on parent of input
  if (parent == "mother") {
    phenotype_non_donor$age <- phenotype_non_donor$patient_age_cycle
  } else if (parent == "father") {
    phenotype_non_donor$age <- phenotype_non_donor$partner_age_cycle
  }
  
  # Round parental ages for plotting
  phenotype_non_donor$rounded_age <- round(phenotype_non_donor$age)
  
  # Keep only ages with at least 25 individuals in the dataset 
  num_each_age <- phenotype_non_donor %>% count(rounded_age)
  num_each_age_sufficient <- num_each_age[num_each_age$n >= 25]
  phenotype_non_donor <- phenotype_non_donor %>% 
    filter(rounded_age %in% num_each_age_sufficient$rounded_age)
  
  # Calculate proportion aneuploid embryos 
  phenotype_non_donor$aneuploidy_ratio <- 
    phenotype_non_donor$aneu_true / phenotype_non_donor$total_embryos
  
  
  # Calculate model for plotting
  
  # Calculate mean proportion of aneuploidy for each age
  age_tranched <- phenotype_non_donor %>%
    group_by(rounded_age) %>%
    summarize(
      avg_proportion = mean(aneuploidy_ratio),
      count = n(),  # Count of observations for each age
      .groups = "drop"
    )
  
  # Compute logistic model 
  model <- glm(
    aneuploidy_ratio ~ rounded_age + I(rounded_age^2),
    data = phenotype_non_donor,
    family = binomial(link = "logit")
  )
  
  # Add the model predictions to the aggregated data
  age_tranched$predicted <- predict(
    model, 
    newdata = age_tranched, 
    type = "response"
  )
  
  # Return data with the age info and model 
  return(age_tranched)
  
}

# Call function 
maternal_data <- reformat_data(maternal_phenotype, "mother")
paternal_data <- reformat_data(paternal_phenotype, "father")

# Plot for both parents 
ggplot() +
  geom_point(data = maternal_data, aes(x = rounded_age, y = avg_proportion), color = "#A020F0") +  
  geom_line(data = maternal_data, aes(x = rounded_age, y = predicted), color = "#A020F0", size = 1) +
  geom_point(data = paternal_data, aes(x = rounded_age, y = avg_proportion), color = "#008B8B") + 
  geom_line(data = paternal_data, aes(x = rounded_age, y = predicted), color = "#008B8B", size = 1) +
  labs(
    title = "Aneuploidy by Age",
    x = "Parental Age",
    y = "Proportion of Embryos with Aneuploidy"
  ) +
  theme_minimal() + 
  ylim(0, 1)

