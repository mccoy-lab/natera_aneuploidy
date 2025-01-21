#!python3

# =================
# author: Arjun Biddanda, Biology Dept., Johns Hopkins University
# email: abiddan1@jhu.edu
# last update: March 14, 2024
# aim: Estimate aneuploidy across autosomes using karyoHMM.
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
meanr_file = (
    "/scratch16/rmccoy22/abiddan1/natera_spectrum/deprecated/r_expected/meanR_child.txt"
)
strand_file = "/data/rmccoy22/natera_spectrum/data/illumina_files/humancytosnp-12v2-1_h-b37.strand"
strand_refalt = "/data/rmccoy22/natera_spectrum/data/illumina_files/humancytosnp-12v2-1_h-b37.strand.RefAlt"
cytosnp_map_v12 = (
    "/data/rmccoy22/natera_spectrum/data/illumina_files/snp_map_cyto12b_f004.txt"
)

# Create the VCF data dictionary for each chromosome ...
vcf_dict = {}
chroms = [f"chr{i}" for i in range(1, 23)]
for i, c in enumerate(range(1, 23)):
    vcf_dict[chroms[i]] = (
        f"/data/rmccoy22/natera_spectrum/genotypes/opticall_parents_100423/genotypes/eagle_phased_hg38/natera_parents.b38.chr{c}.vcf.gz"
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
            total_data.append(f"results/natera_inference/{m}+{f}/{c}.total.ploidy.tsv")

# ------- Rules Section ------- #


localrules:
    all,
    hmm_model_chromosomes,
    generate_posterior_table,


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


rule preprocess_baf_data:
    """Preprocess the BAF data for a trio estimating posteriors downstream

    This rule performs the following steps:

    1. Assess that the csv.gz corresponding to the proper embryo has the correct columns

    2. Obtain the two diploid (phased) genomes for the mother & father at variants with MAF > 0.01 for speed

    Then for each variant on the chromosome perform:

    *  Assess if an allele is a valid SNP nucleotide (ACGT)

    * Assign BAF to be reflective of the alternative allele frequency in the parental VCF (accounting for the complement strand)

    * Filter out positions where either parent has a missing genotype
    
    The specific steps can be found in `preprocess_natera.py` in greater detail with code examples as well.
    """
    input:
        metadata_csv=metadata_file,
        cytosnp_map=cytosnp_map_v12,
        alleles_file=alleles_file,
        egt_cluster=cluster_file,
        meanr_file=meanr_file,
        vcf_file=[vcf_dict[c] for c in chroms],
        child_data=lambda wildcards: find_child_data(wildcards.child_id)[0],
    output:
        baf_pkl="results/natera_inference/{mother_id}+{father_id}/{child_id}.bafs.pkl.gz",
    resources:
        time="1:00:00",
        mem_mb="5G",
    wildcard_constraints:
        chrom="|".join(chroms),
    params:
        chroms=chroms,
    script:
        "scripts/preprocess_natera.py"


rule hmm_model_comparison:
    """Apply the ploidy HMM to the pre-processed BAF data for this embryo."""
    input:
        baf_pkl="results/natera_inference/{mother_id}+{father_id}/{child_id}.bafs.pkl.gz",
    output:
        hmm_pkl="results/natera_inference/{mother_id}+{father_id}/{child_id}.hmm_model.pkl.gz",
    resources:
        time="1:00:00",
        mem_mb="4G",
    params:
        unphased=False,
        mother_id=lambda wildcards: f"{wildcards.mother_id}",
        father_id=lambda wildcards: f"{wildcards.father_id}",
        child_id=lambda wildcards: f"{wildcards.child_id}",
        chroms=chroms,
    script:
        "scripts/baf_hmm_bulk.py"


def bayes_factor(posteriors, priors=None):
    """Compute Bayes Factors for evidence of specific aneuploidy states."""
    if priors is None:
        priors = np.ones(posteriors.size) / posteriors.size

    assert posteriors.size == priors.size
    assert np.isclose(np.sum(priors), 1.0)
    bfs = np.zeros(posteriors.size)
    for i in range(posteriors.size):
        denom = np.sum(
            [posteriors[j] * priors[i] for j in range(posteriors.size) if j != i]
        )
        bfs[i] = posteriors[i] * (1 - priors[i]) / denom
    return bfs


rule hmm_model_chromosomes:
    """Local rule that collapses all ploidy assignments into a single table."""
    input:
        hmm_models="results/natera_inference/{mother_id}+{father_id}/{child_id}.hmm_model.pkl.gz",
    output:
        ploidy="results/natera_inference/{mother_id}+{father_id}/{child_id}.total.ploidy.tsv",
    resources:
        time="0:30:00",
        mem_mb="1G",
    run:
        with open(output.ploidy, "w") as out:
            full_hmm_output = pickle.load(gzip.open(input.hmm_models, "r"))
            out.write(
                "mother\tfather\tchild\tchrom\tsigma_baf\tpi0_baf\t0\t1m\t1p\t2\t3m\t3p\tbf_max\tbf_max_cat\tprop01_exp_het\n"
            )
            cats = np.array(["0", "1m", "1p", "2", "3m", "3p"])
            for c in chroms:
                data = full_hmm_output[c]
                post_vals = np.array([data[x] for x in cats])
                bayes_factor_chrom = bayes_factor(post_vals)
                max_bf = np.max(bayes_factor_chrom)
                max_cat = cats[np.argmax(bayes_factor_chrom)]
                out.write(
                    f"{data['mother_id']}\t{data['father_id']}\t{data['child_id']}\t{c}\t{data['sigma_baf']}\t{data['pi0_baf']}\t{data['0']}\t{data['1m']}\t{data['1p']}\t{data['2']}\t{data['3m']}\t{data['3p']}\t{max_bf}\t{max_cat}\t{data['prop01_exp_het']}\n"
                )


rule generate_posterior_table:
    """Generates a full TSV with posterior probabilities for each embryo across ploidy states.
    
    The columns of the TSV contain the specific karyotype states like 0, 1m, 1p, 2m, 2p, 2, 3m, 3p.
    
    Each row corresponds to a specific SNP position.
    """
    input:
        baf_data="results/natera_inference/{mother_id}+{father_id}/{child_id}.bafs.pkl.gz",
        hmm_models="results/natera_inference/{mother_id}+{father_id}/{child_id}.hmm_model.pkl.gz",
        ploidy="results/natera_inference/{mother_id}+{father_id}/{child_id}.total.ploidy.tsv",
    output:
        posterior="results/natera_inference/{mother_id}+{father_id}/{child_id}.total.posterior.tsv.gz",
    resources:
        time="0:30:00",
        mem_mb="1G",
    run:
        full_hmm_output = pickle.load(gzip.open(input.hmm_models, "r"))
        full_baf = pickle.load(gzip.open(input.baf_data, "r"))
        for c in tqdm(chroms):
            data = full_hmm_output[c]
            baf_data = full_baf[c]
            try:
                gammas = data["gammas"]
                cur_df = pd.DataFrame(gammas.T)
                cur_df.columns = data["states"]
                cur_df["chrom"] = c
                cur_df["pos"] = baf_data["pos"]
                cur_df["rsid"] = baf_data["rsids"]
                tot_dfs.append(cur_df)
            except:
                pass
        if tot_dfs == []:
            df = pd.DataFrame(
                columns=["chrom", "pos", "rsid", "0", "1m", "1p", "2", "3m", "3p"]
            )
        else:
            df = pd.concat(tot_dfs)
            cols_to_move = ["chrom", "pos", "rsid"]
            df = df[
                cols_to_move + [col for col in df.columns if col not in cols_to_move]
            ]
        df.to_csv(output.posterior, sep="\t", index=None)
