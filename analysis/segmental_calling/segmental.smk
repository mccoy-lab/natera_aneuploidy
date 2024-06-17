#!python3

import numpy as np
import pandas as pd

import pickle, gzip
from tqdm import tqdm
from pathlib import Path
from io import StringIO


# ---- Parameters for inference in Natera Data ---- #
metadata_file = "../../data/spectrum_metadata_merged.csv"
aneuploidy_calls = "/data/rmccoy22/natera_spectrum/karyohmm_outputs/compiled_output/natera_embryos.karyohmm_v30a.031624.tsv.gz"

# Create the VCF data dictionary for each chromosome ...
vcf_dict = {}
chroms = [f"chr{i}" for i in range(1, 23)]
for i, c in enumerate(range(1, 23)):
    vcf_dict[
        chroms[i]
    ] = f"/data/rmccoy22/natera_spectrum/genotypes/opticall_parents_100423/genotypes/eagle_phased_hg38/natera_parents.b38.chr{c}.vcf.gz"

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


# Run inference for all valid triplets (restricting to euploid chromosomes)
total_data = []
if Path("results/segmental_inference/valid_trios.txt").is_file():
    with open("results/segmental_inference/valid_trios.txt", "r") as fp:
        for i, line in enumerate(fp):
            [m, f, c] = line.rstrip().split()
            if m != "mother":
                total_data.append(
                    f"results/segmental_inference/{m}+{f}+{c}.segmental.tsv"
                )
        total_data = np.unique(total_data).tolist()


# ------- Rules Section ------- #
localrules:
    all,
    generate_parent_sample_list,
    obtain_valid_trios,


rule all:
    input:
        "results/segmental_inference/valid_trios.txt",
        total_data,


rule generate_parent_sample_list:
    """Generate the parental samples list."""
    input:
        vcf=[vcf_dict[c] for c in chroms],
    output:
        "results/segmental_inference/parent_samples.txt",
    resources:
        time="0:30:00",
        mem_mb="1G",
    shell:
        "bcftools query -l {input.vcf} | sort | uniq > {output}"


rule obtain_valid_trios:
    """Obtain the valid trios here."""
    input:
        metadata_tsv=metadata_file,
        parent_samples="results/segmental_inference/parent_samples.txt",
    output:
        valid_trios="results/segmental_inference/valid_trios.txt",
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
            out.write("mother\tfather\tchild\n")
            for m, f, c in valid_trios:
                out.write(f"{m}\t{f}\t{c}\n")


def define_triplets(
    mother_id,
    father_id,
    child_id,
    trio_file="results/segmental_inference/valid_trios.txt",
    base_path="/home/abiddan1/scratch16/natera_aneuploidy/analysis/aneuploidy/results/natera_inference",
):
    trio_df = pd.read_csv(trio_file, sep="\t")
    filt_df = trio_df[
        (trio_df.mother == mother_id)
        & (trio_df.father == father_id)
        & (trio_df.child == child_id)
    ]
    assert filt_df.shape[0] > 0
    return f"{base_path}/{mother_id}+{father_id}/{child_id}.hmm_model.pkl.gz"


rule est_segmental_changepoints:
    """Estimate segmental aneuploidy using changepoint detection on the HMM trace."""
    input:
        triplets="results/segmental_inference/valid_trios.txt",
        hmm_trace=lambda wildcards: define_triplets(
            mother_id=wildcards.mother,
            father_id=wildcards.father,
            child_id=wildcards.child,
        ),
    output:
        changept_calls="results/segmental_inference/{mother}+{father}+{child}.segmental.tsv",
    resources:
        time="1:00:00",
        mem_mb="5G",
    params:
        chroms=chroms,
        min_size=20,
        penalty=10,
    script:
        "scripts/changept_segmental.py"
