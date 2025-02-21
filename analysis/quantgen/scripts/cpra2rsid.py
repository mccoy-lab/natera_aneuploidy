""" 
original author: Arjun Biddanda, Biology Dept., Johns Hopkins University 
modified by: Sara A. Carioscia, Biology Dept., Johns Hopkins University 
email: scarios1@jhu.edu 
last update: February 8, 2025
aim: Convert CPRA to rsids in Natera GWAS summary stats (aneuploidy and recombination)
  for use in ldsc pipeline.  
"""

import sys
import gzip

def parse_gzipped_vcf(file_path):
    cpra_rsid_dict = {}
    
    with gzip.open(file_path, 'rt') as f:
        for line in f:
            chrom, pos, ref, alt, rsid = line.strip().split('\t')
            
            # Handle multiple alternative alleles
            for alt_allele in alt.split(','):
                cpra = f"chr{chrom}:{pos}:{ref}:{alt_allele}"
                cpra_rsid_dict[cpra] = rsid
    
    return cpra_rsid_dict

def annotate_cpra_file(input_file, output_file, cpra_rsid_dict, file_type):
    open_file = gzip.open if file_type == "ld_scores" else open
    with open_file(input_file, 'rt') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if file_type == "summary_stats":    
                cpra = line.strip().split('\t')[0]  # First column is CPRA
            elif file_type == "ld_scores":
                cpra = line.strip().split('\t')[1] # Second column is CPRA 
            rsid = cpra_rsid_dict.get(cpra, 'NA')  # Lookup RSID or assign 'NA'
            outfile.write(f"{rsid}\t{line.strip()}\n")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python3 script.py <gzipped_input_file> <input_file> <output_file> <filetype>")
        sys.exit(1)
    
    file_path = sys.argv[1]
    input_file = sys.argv[2]
    output_file = sys.argv[3]
    file_type = sys.argv[4]
    
    cpra_rsid = parse_gzipped_vcf(file_path)

    print("First few CPRA:RSID mappings:")
    for i, (cpra, rsid) in enumerate(cpra_rsid.items()):
        if i >= 10:
            break
        print(cpra, "->", rsid)

    annotate_cpra_file(input_file, output_file, cpra_rsid, file_type)

