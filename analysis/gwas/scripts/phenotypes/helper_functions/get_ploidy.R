## Script called in each large-scale ploidy phenotype R file

# function to filter embryo data
filter_data <- function(embryos, bayes_factor_cutoff) {

	# Check if bayes_factor_cutoff is numeric and positive
  	if (!is.numeric(bayes_factor_cutoff) || bayes_factor_cutoff <= 0) {
    	stop("Invalid 'bayes_factor_cutoff' argument. It should be a positive numeric value.")
  	}
    
    # keep only rows that have probabilities for all 6 cn states
    embryos <- embryos[complete.cases(embryos[,c("0", "1m", "1p", "2", "3m", "3p")]),]
    # keep only rows that met the threshold for bayes factor qc
    embryos <- embryos[embryos$bf_max > bayes_factor_cutoff,]
    
    return(embryos)
}


count_ploidy_by_parent <- function(data, parent, phenotype, ploidy_threshold) {
	# check input 
	if (!(parent %in% c("mother", "father"))) {
    	stop("Invalid 'parent' argument. Use 'mother' or 'father'.")
	}

    # choose relevant aneuploidies based on phenotype 
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
    
    # group ploidy by respective parent
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
	# check input 
	if (!(parent %in% c("mother", "father"))) {
    	stop("Invalid 'parent' argument. Use 'mother' or 'father'.")
	}

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

run_phenotype <- function(embryos, bayes_factor_cutoff, parent, phenotype, ploidy_threshold) {
	# filter embryo data 
    embryos_filtered <- filter_data(embryos, bayes_factor_cutoff)

	# group ploidy by respective parent 
	ploidy_counts_by_parent <- count_ploidy_by_parent(embryos, !!as.name(parent), phenotype, ploidy_threshold)
	colnames(ploidy_counts_by_parent)[1] <- "array"

	return(ploidy_counts_by_parent)
}
