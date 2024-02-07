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
aneuploidy_bph_sph_calls_merged = "/data/rmccoy22/natera_spectrum/karyohmm_outputs/compiled_output/natera_embryos.karyohmm_v18.bph_sph_trisomy.102523.tsv.gz"
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
        res_files.append(f"results/bph_sph/inferred/{m}+{p}+{c}.{chrom}.tsv")
    return res_files


def expand_mosaic_est(fp="results/mosaic_est/valid_mosaics.tsv"):
    mosaic_df = pd.read_csv(fp, sep="\t")
    res_files = []
    for m, p, c, chrom in zip(
        trisomy_df.mother.values,
        trisomy_df.father.values,
        trisomy_df.child.values,
        trisomy_df.chrom.values,
    ):
        res_files.append(f"results/mosaic_est/inferred/{m}+{p}+{c}.{chrom}.tsv")
    return res_files


# ---- Target definition ---- #
# TARGETS = ["results/bph_sph/valid_trisomies.tsv"]
# if Path("results/bph_sph/valid_trisomies.tsv").is_file():
# TARGETS.append(expand_bph_sph())

# if Path(aneuploidy_bph_sph_calls_merged).is_file():
# TARGETS.append("results/filt_aneuploidy.tsv.gz")


rule all:
    input:
        "results/filt_aneuploidy.tsv.gz",


# -------- 1. Isolate BPH vs. SPH signature of trisomies ---------- #
rule isolate_trisomies:
    input:
        aneuploidy_calls=aneuploidy_calls,
    output:
        trisomy_tsv="results/bph_sph/valid_trisomies.tsv",
    params:
        postThreshold=0.90,
    run:
        ppTrisomy = float(params.postThreshold)
        aneuploidy_df = pd.read_csv(input.aneuploidy_calls, sep="\t")
        assert "bf_max_cat" in aneuploidy_df.columns
        assert "mother" in aneuploidy_df.columns
        assert "father" in aneuploidy_df.columns
        assert "child" in aneuploidy_df.columns
        assert "chrom" in aneuploidy_df.columns
        assert "3m" in aneuploidy_df.columns
        assert "3p" in aneuploidy_df.columns
        trisomy_df = aneuploidy_df[
            (aneuploidy_df["3m"] >= ppTrisomy) | (aneuploidy_df["3p"] >= ppTrisomy)
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
        bph_tsv="results/bph_sph/inferred/{mother}+{father}+{child}.{chrom}.tsv",
    resources:
        time="0:10:00",
        mem_mb="5G",
    params:
        bp_padding=10e6,
    script:
        "scripts/bph_vs_sph.py"


rule aggregate_bph_sph:
    """Aggregate the BPH vs. SPH signature."""
    input:
        trisomy_tsv="results/bph_sph/valid_trisomies.tsv",
        bph_sph_results=expand_bph_sph(),
    output:
        aggregate_bph_sph="results/bph_sph/natera.total.bph_sph.tsv.gz",
    shell:
        "find results/bph_sph/inferred/ -name \"*.tsv\" | while read line; do cat $line; done | awk '!visited[$0]++' | gzip > {output.aggregate_bph_sph}"


# ----------- 2. Mosaic Estimation routines -------------------- #
rule isolate_putative_mosaics:
    input:
        aneuploidy_calls=aneuploidy_calls,
    output:
        mosaic_tsv="results/mosaic_est/valid_mosaics.tsv",
    params:
        postThreshold=0.90,
    run:
        ppThresh = float(params.postThreshold)
        aneuploidy_df = pd.read_csv(input.aneuploidy_calls, sep="\t")
        assert "bf_max_cat" in aneuploidy_df.columns
        assert "mother" in aneuploidy_df.columns
        assert "father" in aneuploidy_df.columns
        assert "child" in aneuploidy_df.columns
        assert "chrom" in aneuploidy_df.columns
        assert "3m" in aneuploidy_df.columns
        assert "3p" in aneuploidy_df.columns
        assert "1m" in aneuploidy_df.columns
        assert "1p" in aneuploidy_df.columns
        assert "2" in aneuploidy_df.columns
        assert "0" in aneuploidy_df.columns
        aneu_cats = ["0", "1m", "1p", "2", "3m", "3p"]
        mosaic_df = aneuploidy_df[
            np.max(aneuploidy_df[aneu_cats].values, axis=0) <= ppThresh
        ][["mother", "father", "child", "chrom"]]
        mosaic_df.to_csv(output.mosaic_tsv, sep="\t", index=None)


rule mosaic_est:
    """Estimate Mosaic Cell Fraction."""
    input:
        baf_pkl=lambda wildcards: f"{results_dir}/{wildcards.mother}+{wildcards.father}/{wildcards.child}.bafs.pkl.gz",
        mosaic_tsv="results/mosaic_est/valid_mosaics.tsv",
    output:
        mosaic_tsv="results/mosaic_est/{mother}+{father}+{child}.{chrom}.tsv",
    resources:
        time="0:10:00",
        mem_mb="5G",
    script:
        "scripts/mosaic_est.py"


rule aggregate_mosaic_est:
    """Aggregate the Mosaic estimation signature."""
    input:
        mosaic_tsv="results/mosaic_est/valid_mosaics.tsv",
        bph_sph_results=expand_mosaic_est(),
    output:
        aggregate_mosaic="results/mosaic_est/natera.total.mosaic_est.tsv.gz",
    shell:
        "find results/mosaic_est/inferred/ -name \"*.tsv\" | while read line; do cat $line; done | awk '!visited[$0]++' | gzip > {output.aggregate_mosaic_est}"


# ----------- 2a. Merging aneuploidy calls + BPH SPH + Mosaic Estimates -------- #


# ----------- 3. Aneuploidy full filtering mechanism ----------- #
rule run_aneuploidy_filtering:
    """Run the script to generate a filtered set of calls to be used in downstream analyses."""
    input:
        aneuploidy_tsv=aneuploidy_bph_sph_calls_merged,
        meta_csv=metadata_file,
    output:
        filt_aneuploidy_tsv="results/filt_aneuploidy.tsv.gz",
    params:
        sd=3,
        k=3,
        q=10,
    script:
        "scripts/aneuploidy_filtering.py"
