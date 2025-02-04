# =================
# author: Sara A. Carioscia, Biology Dept., Johns Hopkins University
# email: scarios1@jhu.edu
# last update: February 3, 2025
# aim: create file of CPRA:RSID dictionary for use in renaming summary stats 
# =================

# load arguments
args <- commandArgs(trailingOnly = TRUE)

# check if the correct number of arguments is provided
if (length(args) < 2) {
  stop("Usage: Rscript build_rsid_dict.R <dbsnp_file> <out_fname>")
}
dbsnp_file <- args[1]  # Path to the dbSNP file (gzipped)
out_fname <- args[2]  # Path to the output file

# function to build an RSID dictionary from a gzipped dbSNP file
build_rsid_dict <- function(dbsnp_fp) {
  cpra2rsid <- list()
  
  # read the file into a character vector
  lines <- readLines(gzfile(dbsnp_fp))
  
  for (line in lines) {
    # skip VCF header lines starting with ##, and colnames line starting with #
    if (startsWith(line, "#")) next  
    
    fields <- strsplit(line, "\\s+")[[1]]
    chrom <- fields[1]
    pos <- fields[2]
    ref <- fields[4]
    # create one CPRA for each comma-separated ALT allele option 
    alts <- strsplit(fields[5], ",")[[1]]
    rsid <- fields[3]
    
    for (a in alts) {
      cpra <- paste0("chr", chrom, ":", pos, ":", ref, ":", a)
      cpra2rsid[[cpra]] <- rsid
    }
  }
  
  return(cpra2rsid)
}

# build RSID dictionary
dict_dbsnp <- build_rsid_dict(dbsnp_file)

# write dictionary to tab-separated file
write.table(
  data.frame(CPRA = names(dict_dbsnp), SNP = unlist(dict_dbsnp)),
  file = out_fname,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE
)

