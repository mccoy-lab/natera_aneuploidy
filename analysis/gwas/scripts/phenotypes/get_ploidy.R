## Script called in each large-scale ploidy phenotype R file

# function to filter embryo data
filter_data <- function(embryos, bayes_factor_cutoff) {
	# keep only rows that have probabilities for all 6 cn states
	embryos <- embryos[complete.cases(embryos[,c("0", "1m", "1p", "2", "3m", "3p")]),]
	# keep only rows that met the threshold for bayes factor qc
	embryos <- embryos[embryos$bf_max > bayes_factor_cutoff,]

	return(embryos)
}


count_ploidy_by_parent <- function(data, parent, phenotype, ploidy_threshold) {
  	if (phenotype == "maternal_triploidy") {
		cn <- "3m"
	} else if (phenotype == "paternal_triploidy") {
		cn <- "3p"
	} else if (phenotype == "maternal_haploidy") {
		cn <- "1p"
	} else if (phenotype == "paternal_haploidy") {
		cn <- "1m"
	} else {
		stop("Invalid 'phenotype' argument.")
	}

  result <- data %>%
    group_by({{parent}}, child) %>%
    summarise(num_affected = sum(bf_max_cat == cn)) %>%
		mutate(is_ploidy = if_else(num_affected >= ploidy_threshold, "aneu_true", "aneu_false")) %>%
		count(is_ploidy) %>%
		pivot_wider(names_from = is_ploidy, values_from = n, values_fill = 0) %>%
		replace(is.na(.), 0)
  
  return(result)
}

count_complex_ploidy_by_parent <- function(data, parent) {
	# group ploidy by respective parent 
	result <- embryos %>%
		group_by({{parent}}, child) %>%
		summarise(parents_affected = sum(any(bf_max_cat %in% c("3m", "1p")) & any(bf_max_cat %in% c("1m", "3p")))) %>% 
		mutate(is_ploidy = if_else(parents_affected > 1, "aneu_true", "aneu_false")) %>%
		count(is_ploidy) %>%
		pivot_wider(names_from = is_ploidy, values_from = n, values_fill = 0) %>%
		replace(is.na(.), 0)

	return(result)
}
