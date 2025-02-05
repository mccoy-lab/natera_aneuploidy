# =================
# author: Sara A. Carioscia, Biology Dept., Johns Hopkins University
# email: scarios1@jhu.edu
# last update: February 3, 2025
# aim: create file of CPRA:RSID dictionary for use in renaming summary stats 
# =================

# load libraries
library(parallel)

# load arguments
args <- commandArgs(trailingOnly = TRUE)

# check if the correct number of arguments is provided
if (length(args) < 3) {
  stop("Usage: Rscript build_rsid_dict.R <dbsnp_file> <out_fname>")
}

dbsnp_file <- args[1] # path and fname to the dbSNP file (gzipped)
num_cores <- as.numeric(args[2])  # number of cores for use in parallel processing 
out_fname <- args[3]  # path and fname to the output file

if (is.na(num_cores) || num_cores <= 0) {
  stop("num_cores is not valid. It must be a positive numeric value.")
}

# function to build an RSID dictionary from a gzipped dbSNP file
build_rsid_dict <- function(dbsnp_fp) {
  
  cpra2rsid <- list()
  
  # read the file into a character vector
  lines <- readLines(gzfile(dbsnp_fp))
  
  # filter out header lines
  lines <- lines[!startsWith(lines, "#")]
  
  # function to process a chunk of lines
  process_chunk <- function(chunk) {
    local_cpra2rsid <- list()  # To store results for this chunk
    
    for (line in chunk) {
      fields <- strsplit(line, "\\s+")[[1]]
      chrom <- fields[1]
      pos <- fields[2]
      ref <- fields[4]
      alts <- strsplit(fields[5], ",")[[1]]
      rsid <- fields[3]
      
      for (a in alts) {
        cpra <- paste0("chr", chrom, ":", pos, ":", ref, ":", a)
        local_cpra2rsid[[cpra]] <- rsid
      }
    }
    return(local_cpra2rsid)  
  }
  
  # Ensure that num_lines is numeric
  num_lines <- length(lines)
  if (!is.numeric(num_lines)) {
    stop("num_lines is not numeric")
  }
  
  # group lines into chunks for parallel processing
  num_lines <- length(lines)
  chunk_size <- ceiling(num_lines / num_cores)
  line_chunks <- split(lines, ceiling(seq_along(lines) / chunk_size))
  
  # use mclapply to process chunks in parallel
  results <- mclapply(line_chunks, process_chunk, mc.cores = num_cores)
  
  # Combine the results from all chunks
  cpra2rsid <- do.call(c, results)  # Combine the list from each chunk
  
  return(cpra2rsid)
}

# build RSID dictionary using parallel processing
dict_dbsnp <- build_rsid_dict(dbsnp_file)

# remove automatic names from split() function 
dict_dbsnp_clean <- data.frame(CPRA = unname(names(dict_dbsnp)), SNP = unlist(dict_dbsnp))

# write dictionary to tab-separated file
write.table(
  data.frame(CPRA = names(dict_dbsnp_clean), SNP = unlist(dict_dbsnp_clean)),
  file = out_fname,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE
)

