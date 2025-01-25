#!python3

# =================
# author: Arjun Biddanda, Biology Dept., Johns Hopkins University
# email: abiddan1@jhu.edu
# last update: Nov 10, 2024
# aim: Estimate aneuploidy on sex-chromosomes.
# =================


import numpy as np
import pandas as pd
import pickle, gzip
from tqdm import tqdm
from pathlib import Path
from io import StringIO

# ---- Parameters for inference in Natera Data ---- #
metadata_file = "../../data/spectrum_metadata_merged.csv"
alleles_file = "/data/rmccoy22/natera_spectrum/data/illumina_files/humancytosnp-12v2-1_h.update_alleles.txt"
cluster_file = "/scratch16/rmccoy22/abiddan1/natera_spectrum/r_expected/HumanCytoSNP-12v2-1_NS550.cluster.tsv.gz"
strand_file = "/data/rmccoy22/natera_spectrum/data/illumina_files/humancytosnp-12v2-1_h-b37.strand"
strand_refalt = "/data/rmccoy22/natera_spectrum/data/illumina_files/humancytosnp-12v2-1_h-b37.strand.RefAlt"
cytosnp_map_v12 = (
    "/data/rmccoy22/natera_spectrum/data/illumina_files/snp_map_cyto12b_f004.txt"
)
lrrs = ["none"]

# Create the VCF data dictionary for the sex-chromosomes
vcf_dict = {}
chroms = ["chrX", "chrY"]
# vcf_dict["chrX"] = (
# "/data/rmccoy22/natera_spectrum/genotypes/opticall_parents_100423/genotypes/opticall_concat_23.norm.b38.vcf.gz"
# )
vcf_dict["chrX"] = (
    "/data/rmccoy22/natera_spectrum/genotypes/opticall_parents_100423/genotypes/eagle_phased_hg38/natera_parents.b38.beagle.chr23.annot.vcf.gz"
)


vcf_dict["chrY"] = (
    "/data/rmccoy22/natera_spectrum/genotypes/opticall_parents_100423/genotypes/opticall_concat_24.norm.b38.vcf.gz"
)

# Read in the aggregate metadata file
meta_df = pd.read_csv(metadata_file)


def find_child_data(
    child_id, meta_df=meta_df, raw_data_path="/data/rmccoy22/natera_spectrum/data/"
):
    """Find the child csv file based on the provided meta_df."""
    child_year = meta_df[meta_df.array == child_id].year.values
    for year in child_year:
        child_fp = Path(f"{raw_data_path}{year}/{child_id}.csv.gz")
        if child_fp.is_file():
            return child_fp, True
        else:
            continue
    return None, False


def create_trios(
    meta_df, sample_file, raw_data_path="/data/rmccoy22/natera_spectrum/data/"
):
    """Create a list of valid trio datasets."""
    assert "family_position" in meta_df.columns
    assert "casefile_id" in meta_df.columns
    assert "array" in meta_df.columns
    grouped_df = (
        meta_df.groupby(["casefile_id", "family_position"])["array"]
        .agg(lambda x: list(x))
        .reset_index()
    )
    valid_trios = []
    for case in tqdm(np.unique(grouped_df.casefile_id)):
        cur_df = grouped_df[grouped_df.casefile_id == case]
        for m in cur_df[cur_df.family_position == "mother"].array.values[0]:
            for f in cur_df[cur_df.family_position == "father"].array.values[0]:
                for c in cur_df[cur_df.family_position == "child"].array.values[0]:
                    valid_trios.append((case, m, f, c))

    valid_df = pd.DataFrame(
        valid_trios, columns=["casefile_id", "mother", "father", "child"]
    )
    parents = [line.rstrip() for line in open(sample_file, "r")]
    valid_df["parents_in"] = valid_df.mother.isin(parents) & valid_df.father.isin(
        parents
    )
    valid_df["child_found"] = [
        find_child_data(c)[1] for c in tqdm(valid_df.child.values)
    ]
    valid_filt_df = valid_df[
        valid_df.parents_in & valid_df.child_found
    ].drop_duplicates()[["mother", "father", "child"]]
    valid_filt_trios = [
        (m, f, c)
        for (m, f, c) in zip(
            valid_filt_df.mother.values,
            valid_filt_df.father.values,
            valid_filt_df.child.values,
        )
    ]
    return valid_filt_trios


total_data = []
if Path("results/natera_inference/valid_trios.txt").is_file():
    with open("results/natera_inference/valid_trios.txt", "r") as fp:
        for i, line in enumerate(fp):
            [m, f, c] = line.rstrip().split()
            total_data.append(
                f"results/natera_inference/{m}+{f}/{c}.sex_chroms.total_ploidy.tsv"
            )


# ------- Rules Section ------- #
localrules:
    all,


rule all:
    input:
        "results/natera_inference/valid_trios.txt",
        total_data,


rule generate_parent_sample_list:
    """Generate the parental samples list."""
    input:
        vcf=[vcf_dict[c] for c in chroms],
    output:
        "results/natera_inference/parent_samples.txt",
    resources:
        time="0:30:00",
        mem_mb="1G",
    shell:
        "bcftools query -l {input.vcf} | sort | uniq > {output}"


rule obtain_valid_trios:
    """Obtain the valid trios here."""
    input:
        metadata_tsv=metadata_file,
        parent_samples="results/natera_inference/parent_samples.txt",
    output:
        valid_trios="results/natera_inference/valid_trios.txt",
    resources:
        time="0:30:00",
        mem_mb="1G",
    run:
        valid_trios = create_trios(
            meta_df,
            input.parent_samples,
            raw_data_path="/data/rmccoy22/natera_spectrum/data/",
        )
        with open(output.valid_trios, "w") as out:
            for m, f, c in valid_trios:
                out.write(f"{m}\t{f}\t{c}\n")


rule preprocess_baf_data_sex_chrom:
    """Preprocess the BAF data for a trio estimating posteriors downstream

    This rule performs the following steps:

    1. Assess that the csv.gz corresponding to the proper embryo has the correct columns

    2. Obtain the two diploid (phased) genomes for the mother & father at variants with MAF > 0.01 for speed

    Then for each variant on the chromosome perform:

    *  Assess if an allele is a valid SNP nucleotide (ACGT)

    * Assign BAF to be reflective of the alternative allele frequency in the parental VCF (accounting for the complement strand)

    * Filter out positions where either parent has a missing genotype
    
    The specific steps can be found in `preprocess_sex_chroms.py` in greater detail with code examples as well.
    """
    input:
        metadata_csv=metadata_file,
        cytosnp_map=cytosnp_map_v12,
        alleles_file=alleles_file,
        egt_cluster=cluster_file,
        vcf_file=[vcf_dict[c] for c in chroms],
        child_data=lambda wildcards: find_child_data(wildcards.child_id)[0],
    output:
        baf_pkl="results/natera_inference/{mother_id}+{father_id}/{child_id}.sex_chroms.bafs.pkl.gz",
    resources:
        partition="parallel",
        time="1:00:00",
        mem_mb="5G",
    wildcard_constraints:
        chrom="|".join(chroms),
    params:
        chroms=chroms,
    script:
        "scripts/preprocess_sex_chroms.py"


rule assign_sex_chrom_copy:
    """Apply the revised sex-chromosome ploidy HMM to the pre-processed BAF data for this embryo and output a TSV."""
    input:
        baf_pkl="results/natera_inference/{mother_id}+{father_id}/{child_id}.sex_chroms.bafs.pkl.gz",
    output:
        karyo_tsv="results/natera_inference/{mother_id}+{father_id}/{child_id}.sex_chroms.total_ploidy.tsv",
    resources:
        partition="parallel",
        time="0:10:00",
        mem_mb="2G",
    params:
        mother_id=lambda wildcards: f"{wildcards.mother_id}",
        father_id=lambda wildcards: f"{wildcards.father_id}",
        child_id=lambda wildcards: f"{wildcards.child_id}",
        chroms=chroms,
    script:
        "scripts/sex_chrom_hmm.py"
