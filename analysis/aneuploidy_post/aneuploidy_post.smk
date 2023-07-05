#!python3

import numpy as np
import pandas as pd

import pickle, gzip
from tqdm import tqdm
from pathlib import Path
from io import StringIO

# ---- Parameters for post whole-chromosome aneuploidy inference in Natera Data ---- #
metadata_file = "../../data/spectrum_metadata_merged.csv"
centromeres_file = ""
aneuploidy_calls=""
results_dir=""

# Read in the aggregate metadata file
meta_df = pd.read_csv(metadata_file)

# ------- Rules Section ------- #
localrules:
    all,


rule all:
    input:


rule trisomy_bph_sph:
    """Obtain evidence for BPH vs. SPH for observed trisomies."""
    input:
        hmm_traceback = lambda wildcards: f'{results_dir}/{wildcards.mother}+{wildcards.father}/{wildcards.child}.'
    output:
        "results/bph_sph/{mother}+{father}+{child}.{chrom}.tsv",
    resources:
        time="0:30:00",
        mem_mb="1G",
    script:
        "scripts/bph_vs_sph.py"
