# Number of embryos per mother (overall and per visit)
# Proportion of aneuploid embryos by maternal age (with regression) 

# load libraries
library(data.table)
library(dplyr)
library(ggplot2)

# load phenotype data
embryo_count_by_mother <- fread("/scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/results/phenotypes/embryo_count_by_mother.csv")

# create column with rate of aneuploid embryos 
embryo_count_by_mother <- embryo_count_by_mother %>%
  mutate(aneuploid_euploid_ratio = aneuploid / (euploid + aneuploid))


## EMBRYO COUNT PER MOTHER 

# plot number of overall embryos per mother 
ggplot(embryo_count_by_mother, aes(x = num_embryos)) +
  geom_histogram(binwidth = 1) +
  labs(title = "Number of Embryos per Mother",
       x = "Number of Embryos",
       y = "Number of Mothers") + 
  theme_minimal()

# plot number of overall embryos per mother per visit
ggplot(embryo_count_by_mother, aes(x = num_embryos/num_visits)) +
  geom_histogram(binwidth = 1) +
  labs(title = "Number of Embryos per Mother per Visit",
       x = "Number of Embryos",
       y = "Number of Mothers") + 
  theme_minimal()


## ANEUPLOIDY RATIO BY MATERNAL AGE 

# plot proportion of aneuploid embryos by maternal age with logistic regression
ggplot(embryo_count_by_mother, aes(x = weighted_age, y = aneuploid_euploid_ratio)) +
  geom_point() +  
  geom_smooth(method = "glm", method.args = list(family = binomial), se = FALSE) +  # Logistic regression line
  labs(title = "Proportion of Aneuploid Embryos by Maternal Age",
       x = "Maternal Age",
       y = "Proportion of Aneuploid Embryos") + 
  theme_minimal()

# plot proportion of aneuploid embryos by maternal age with logistic regression - no olds 
ggplot(embryo_count_by_mother[embryo_count_by_mother$weighted_age <= 60], aes(x = weighted_age, y = aneuploid_euploid_ratio)) +
  geom_point() +  
  geom_smooth(method = "glm", method.args = list(family = binomial), se = FALSE) +  # Logistic regression line
  labs(title = "Proportion of Aneuploid Embryos by Maternal Age",
       x = "Maternal Age",
       y = "Proportion of Aneuploid Embryos") + 
  theme_minimal()

# plot proportion of aneuploid embryos by maternal age with stats
ggplot(embryo_count_by_mother[embryo_count_by_mother$weighted_age <= 60], aes(x = weighted_age)) +
  geom_point(aes(y = aneuploid_euploid_ratio), color = "red") +  # Scatter plot for proportion of aneuploid embryos
  geom_smooth(aes(y = aneuploid_euploid_ratio), method = "glm", se = FALSE, color = "blue") +  # Smooth line for proportion of aneuploid embryos
  stat_ecdf(aes(y = aneuploid_euploid_ratio * 100), color = "green", geom = "point") +  # ECDF for proportion of aneuploid embryos
  labs(title = "Proportion of Aneuploid Embryos and its Cumulative Distribution by Maternal Age",
       x = "Maternal Age",
       y = "Proportion of Aneuploid Embryos / Cumulative Proportion (%)") +
  scale_y_continuous(sec.axis = sec_axis(~./100, name = "Cumulative Proportion (%)")) +  # Secondary y-axis for ECDF
  theme_minimal()