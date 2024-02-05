# load libraries
library(data.table)
library(ggplot2)

# Usage: 
# ./plot_ploidy_by_chrom.R # this script 
# /data/rmccoy22/natera_spectrum/karyohmm_outputs/compiled_output/natera_embryos.karyohmm_v18.bph_sph_trisomy.full_annotation.112023.filter_bad_trios.tsv.gz
# /scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/characterize_aneuploidy/ploidy_by_chrom.pdf

args = commandArgs(trailingOnly=TRUE)
input_data <- args[1]
out_fp <- args[2]

# read in karyohmm table 
input_data <- fread(input_data)

# plot aneuploidies by chromosome 
g <- ggplot(input_data[bf_max_cat != "2",], aes(fill = bf_max_cat, chrom))
p <- g + geom_bar()
ggsave(out_fp, p)
