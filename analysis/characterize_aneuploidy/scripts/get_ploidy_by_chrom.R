library(ggplot2)
library(data.table)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)
input_data <- args[1]
out_plot_fp <- args[2]
out_file_fp <- args[3]

#merged_input <- fread("/scratch16/rmccoy22/abiddan1/natera_aneuploidy/sandbox/natera_embryos_v1.karyohmm_v11.051223.tsv")

get_putative_cn <- function(input_data) {

  input_data <- fread(input_data)

  selected_columns <- c("0", "1m", "1p", "2", "3m", "3p")
  highest_values <- apply(input_data[, selected_columns, with = FALSE], 1, function(x) max(x, na.rm = TRUE))
  second_highest_values <- apply(input_data[, selected_columns, with = FALSE], 1, function(x) {
    sorted_x <- sort(x, decreasing = TRUE, na.last = TRUE)
    sorted_x[2]
  })

  input_data[,highest := highest_values]
  input_data[, second_highest := second_highest_values]
  # remove any rows that have NA for all ploidy states (these had baf all OF 0 or 1)
  input_data_hold <- input_data
  input_data_filtered <- input_data[!apply(X = is.na(input_data[,11:16]), MARGIN = 1, FUN = sum) > 0,]
  input_data <- input_data_filtered
  # add column that says what the putative cn is 
  input_data[, putative_cn := colnames(input_data[, 11:16])[apply(input_data[, 11:16], 1, which.max)]]

  # remove "chr" from chromosome for ease of plotting
  input_data_hold <- input_data
  input_data$chromosome <- gsub("chr", "", input_data$chrom) %>% as.numeric()
  
  # return dataframe
  return(input_data)
}

# call function 
putative_cn_dt <- get_putative_cn(input_data)
# write to file
fwrite(putative_cn_dt, out_file_fp, sep = "\t") 

# plot aneuploidies by chromosome 
g <- ggplot(putative_cn_dt[putative_cn != "2",], aes(fill = putative_cn, chromosome))
p <- g + geom_bar()
ggsave(out_plot_fp, p)
