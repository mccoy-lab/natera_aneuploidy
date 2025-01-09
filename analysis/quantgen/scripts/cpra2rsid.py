#!/usr/bin/env python3

import gzip
import argparse


def build_rsid_dict(af_tsv_fp):
    """Function to build a dictionary of rsids."""
    cpra2rsid = {}
    with gzip.open(af_tsv_fp, "rt") as fp:
        #for line in tqdm(fp, desc="Building RSID dictionary"):
        for line in fp:
            [chrom, pos, ref, alts, rsid] = line.rstrip().split()
            for a in alts.split(","):
                cpra = f"chr{chrom}:{pos}:{ref}:{a}"
                cpra2rsid[cpra] = rsid
    return cpra2rsid

def load_rsid_dict(dictionary_fp):
    """Load RSID dictionary from a tab-separated file."""
    cpra2rsid = {}
    with gzip.open(dictionary_fp, "rt") as f:
        for line in f:
            chrom, pos, ref, alt, rsid = line.rstrip().split("\t")  # Split by tab
            cpra = f"chr{chrom}:{pos}:{ref}:{alt}"
            cpra2rsid[cpra] = rsid
    return cpra2rsid


def cpra2rsid(cpra, cpra_alt, rsid_dict):
    """CPRA to rsid direct conversion."""
    if cpra in rsid_dict:
        return rsid_dict[cpra]
    elif cpra_alt in rsid_dict:
        return rsid_dict[cpra_alt]
    else:
        return "NA"  # Use "NA" if no RSID is found


def main():
    # Set up the argument parser
    parser = argparse.ArgumentParser(
        description="Convert CPRA IDs to RSIDs in summary statistics."
    )
    parser.add_argument(
        "--sumstats",
        required=True,
        help="Path to the input summary statistics file.",
    )
    parser.add_argument(
        "--dbsnp",
        required=True,
        help="Path to the dbSNP file (gzipped).",
    )
    parser.add_argument(
        "--dictionary",
        required=True,
        help="Path to the prebuilt RSID dictionary file or specify 'build' to generate a new dictionary.",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Path to the output file with RSIDs.",
    )

    # Parse command-line arguments
    args = parser.parse_args()

    # Load or build the RSID dictionary
    if args.dictionary.lower() == 'build':
        if not args.dbsnp:
            raise ValueError(
                "--dbsnp must be specified to build a dictionary if --dictionary is 'build'."
            )
        print("Building RSID dictionary...")
        rsid_dict = build_rsid_dict(args.dbsnp)
    else:
        print("Loading RSID dictionary from file...")
        rsid_dict = load_rsid_dict(args.dictionary)

    # # Build the RSID dictionary
    # rsid_dict = build_rsid_dict(args.dbsnp)

    # Open the output file for writing
    with open(args.output, "w+") as out:
        # Open the input summary statistics file
        with open(args.sumstats, "r") as sumstats:
            header = sumstats.readline().rstrip()  # Read and clean the header
            out.write(header + "\tRSID\n")  # Add RSID column header
            # Process each line in the summary statistics file
            for line in sumstats:
            #for line in tqdm(sumstats, desc="Processing summary stats"):
                columns = line.split()
                cpra = columns[0].rstrip()  # Get the CPRA from the first column
                cpra_splt = cpra.split(":")  # Split the CPRA into components
                # Construct the alternate CPRA format
                cpra_alt = f"{cpra_splt[0]}:{cpra_splt[1]}:{cpra_splt[3]}:{cpra_splt[2]}"
                
                # Debugging: Print cpra and cpra_alt to check if they're formatted correctly
                print(f"CPRA: {cpra}, CPRA_ALT: {cpra_alt}")
                
                # Get RSID from dictionary
                rsid = cpra2rsid(cpra, cpra_alt, rsid_dict)
                
                # Debugging: Check if the RSID is being correctly retrieved
                print(f"RSID for {cpra}: {rsid}")
                
                # Write the original line and the corresponding RSID
                out.write(line.rstrip() + "\t" + rsid + "\n")


if __name__ == "__main__":
    main()

