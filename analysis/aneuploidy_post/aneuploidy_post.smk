#!python3

import numpy as np
import pandas as pd

import pickle, gzip
from tqdm import tqdm
from pathlib import Path
from io import StringIO


# ---- Parameters for post whole-chromosome aneuploidy inference in Natera Data ---- #
metadata_file = "../../data/spectrum_metadata_merged.csv"
centromeres_file = "../../data/gaps/centromeres_grch38.bed"
aneuploidy_calls = "/data/rmccoy22/natera_spectrum/karyohmm_outputs/compiled_output/natera_embryos.karyohmm_v18.102523.tsv.gz"
results_dir = "../aneuploidy/results/natera_inference"


# ------- Rules Section ------- #
localrules:
    all,


def expand_bph_sph(fp="results/bph_sph/valid_trisomies.tsv"):
    """Expand the BPH vs. SPH results to aggregate."""
    trisomy_df = pd.read_csv(fp, sep="\t")
    res_files = []
    for m, p, c, chrom in zip(
        trisomy_df.mother.values,
        trisomy_df.father.values,
        trisomy_df.child.values,
        trisomy_df.chrom.values,
    ):
        res_files.append(f"results/bph_sph/{m}+{p}+{c}.{chrom}.tsv")
    return res_files


# ---- Target definition ---- #
TARGETS = ["results/bph_sph/valid_trisomies.tsv"]
if Path("results/bph_sph/valid_trisomies.tsv").is_file():
    TARGETS.append(expand_bph_sph())


rule all:
    input:
        TARGETS,


# -------- 1. Isolate BPH vs. SPH signature of trisomies ---------- #
rule isolate_trisomies:
    input:
        aneuploidy_calls=aneuploidy_calls,
    output:
        trisomy_tsv="results/bph_sph/valid_trisomies.tsv",
    run:
        aneuploidy_df = pd.read_csv(input.aneuploidy_calls, sep="\t")
        assert "bf_max_cat" in aneuploidy_df.columns
        assert "mother" in aneuploidy_df.columns
        assert "father" in aneuploidy_df.columns
        assert "child" in aneuploidy_df.columns
        assert "chrom" in aneuploidy_df.columns
        trisomy_df = aneuploidy_df[
            (aneuploidy_df.bf_max_cat == "3m") | (aneuploidy_df.bf_max_cat == "3p")
        ][["mother", "father", "child", "chrom"]]
        trisomy_df.to_csv(output.trisomy_tsv, sep="\t", index=None)


rule trisomy_bph_sph:
    """Obtain evidence for BPH vs. SPH for an observed trisomy in centromere-proximal vs distal regions."""
    input:
        baf_pkl=lambda wildcards: f"{results_dir}/{wildcards.mother}+{wildcards.father}/{wildcards.child}.bafs.pkl.gz",
        hmm_traceback=lambda wildcards: f"{results_dir}/{wildcards.mother}+{wildcards.father}/{wildcards.child}.hmm_model.pkl.gz",
        centromere_bed=centromeres_file,
        trisomy_tsv="results/bph_sph/valid_trisomies.tsv",
    output:
        bph_tsv="results/bph_sph/{mother}+{father}+{child}.{chrom}.tsv",
    resources:
        time="0:10:00",
        mem_mb="5G",
    params:
        bp_padding=10e6,
    script:
        "scripts/bph_vs_sph.py"
