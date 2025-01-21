#!python3


# =================
# author: Arjun Biddanda, Biology Dept., Johns Hopkins University
# email: abiddan1@jhu.edu
# last update: Nov 10, 2024
# aim: evaluate accuracy of aneuploidy detection using simulations.
# =================


import os
import subprocess
import glob
import pandas as pd
from snakemake.utils import validate
from pathlib import Path
import yaml
import warnings
import datetime
import numpy as np

# Mark the datestamp for any date-specific files
DATESTAMP = datetime.datetime.now().strftime("%Y%m%d")

# configure shell behavior for all rules
shell.executable("/bin/bash")
shell.prefix("set -euo pipefail;")

# record commit-ish ID
label = subprocess.check_output(["git", "describe", "--always"]).strip()
print(f"natera-aneuploid simulation workflow {label}")

# create the log dir
Path("logs/").mkdir(parents=True, exist_ok=True)


# Setup and validate a configfile
configfile: "config.yaml"


# Creating the associated targets for simulation & evaluation...
TARGETS = []
if config["hmm_sims"]["model_comp"]:
    TARGETS.append("results/total_hmm_ploidy.tsv.gz")
if config["hmm_sims"]["mixed_ploidy"]:
    TARGETS.append("results/mixed_ploidy_sims.tsv.gz")
    # TARGETS.append("results/mixed_ploidy_sims.mosaic_est.tsv.gz")


localrules:
    all,
    sim_mixed_ploidy,
    hmm_model_chromosomes_mixed,


rule all:
    input:
        TARGETS,


# ---- Simulation 1: Simulate True Whole-chromosome Aneuploidy --------- #
rule sim_baf_lrr_ploidy:
    """Simulate parental haplotypes and BAF for an individual."""
    output:
        baf="results/hmm_simulations/ploidy{k}/sim{rep}_m{m}_pi{pi0}_sigma{sigma}_skew{skew}.a{a}.npz",
    wildcard_constraints:
        k="\d+",
        m="\d+",
        pi0="\d+",
        sigma="\d+",
        skew="\d+",
        a="\d+",
    resources:
        time="0:10:00",
        mem_mb="1G",
    params:
        sfs=config["afs"],
        k=lambda wildcards: int(wildcards.k),
        m=lambda wildcards: int(wildcards.m),
        sigma=lambda wildcards: float(wildcards.sigma) / 100,
        pi0=lambda wildcards: float(wildcards.pi0) / 100,
        seed=lambda wildcards: int(wildcards.rep)
        + int(wildcards.pi0)
        + int(wildcards.sigma),
        mat_skew=lambda wildcards: float(wildcards.skew) / 100,
        mother_id=lambda wildcards: f"k{wildcards.k}_m{wildcards.rep}",
        father_id=lambda wildcards: f"k{wildcards.k}_m{wildcards.rep}",
        child_id=lambda wildcards: f"k{wildcards.k}_m{wildcards.rep}",
        mixed_ploidy=False,
        alpha=lambda wildcards: int(wildcards.a) / 100,
    script:
        "scripts/sim_data.py"


rule hmm_baf_lrr_ploidy:
    input:
        baf=rules.sim_baf_lrr_ploidy.output.baf,
    output:
        hmm_out="results/hmm_ploidy_comp/ploidy{k}/sim{rep}_m{m}_pi{pi0}_sigma{sigma}_skew{skew}.a{a}.phase_error{p}.model_comp.npz",
    wildcard_constraints:
        p="0|1",
    resources:
        time="0:30:00",
        mem_mb="2G",
    params:
        model_comp=True,
        unphased=False,
        phase_error=lambda wildcards: wildcards.p == "1",
        mother_id=lambda wildcards: f"k{wildcards.k}_p{wildcards.p}_m{wildcards.rep}",
        father_id=lambda wildcards: f"k{wildcards.k}_p{wildcards.p}_m{wildcards.rep}",
        child_id=lambda wildcards: f"k{wildcards.k}_p{wildcards.p}_m{wildcards.rep}",
    script:
        "scripts/baf_hmm.py"


rule hmm_combine_baf_lrr:
    """Local rule that collapses all ploidy assignments to a single estimand for mixed-ploidy biopsies."""
    input:
        hmm_model="results/hmm_ploidy_comp/ploidy{k}/sim{rep}_m{m}_pi{pi0}_sigma{sigma}_skew{skew}.a{a}.phase_error{p}.model_comp.npz",
    output:
        ploidy=temp(
            "results/hmm_ploidy_comp/ploidy{k}/sim{rep}_m{m}_pi{pi0}_sigma{sigma}_skew{skew}.a{a}.phase_error{p}.model_comp.tsv"
        ),
    run:
        sigma = int(wildcards.sigma) / 100
        pi0 = int(wildcards.pi0) / 100
        a = int(wildcards.a) / 100
        phase_err = int(wildcards.p) == 1
        with open(output.ploidy, "w") as out:
            data = np.load(input.hmm_model)
            cats = np.array(["0", "1m", "1p", "2", "3m", "3p"])
            posteriors = [data[x] for x in cats]
            max_cat_full = cats[np.argmax(posteriors)]
            out.write(
                "mother\tfather\tchild\taploid\trep\tm\ta\tsigma\tpi0\tsigma_est\tpi0_est\t0\t1m\t1p\t2\t3m\t3p\tmax_cat\n"
            )
            out.write(
                f"{data['mother_id']}\t{data['father_id']}\t{data['child_id']}\t{data['aploid']}\t{wildcards.rep}\t{wildcards.m}\t{a}\t{sigma}\t{pi0}\t{data['sigma_est']}\t{data['pi0_est']}\t{data['0']}\t{data['1m']}\t{data['1p']}\t{data['2']}\t{data['3m']}\t{data['3p']}\t{max_cat_full}\n"
            )


rule collect_hmm_model_baf_lrr:
    input:
        hmm_tsvs=expand(
            "results/hmm_ploidy_comp/ploidy{k}/sim{rep}_m{m}_pi{pi0}_sigma{sigma}_skew{skew}.a{a}.phase_error{p}.model_comp.tsv",
            k=config["hmm_sims"]["simple_sims"]["ploidy"],
            rep=range(1, config["hmm_sims"]["simple_sims"]["reps"] + 1),
            m=config["hmm_sims"]["simple_sims"]["m"],
            pi0=config["hmm_sims"]["simple_sims"]["pi0"],
            sigma=config["hmm_sims"]["simple_sims"]["std_dev"],
            skew=config["hmm_sims"]["simple_sims"]["skew"],
            a=[30],
            p=[1],
        ),
    output:
        tot_hmm_tsv="results/total_hmm_ploidy.tsv.gz",
    run:
        dfs = []
        for p in input.hmm_tsvs:
            df = pd.read_csv(p, sep="\t")
            dfs.append(df)
        tot_df = pd.concat(dfs)
        tot_df.to_csv(output.tot_hmm_tsv, sep="\t", index=None)


# ----- Simulation 2: Simulate Mosaic Whole-Chromosome Aneuploidy -------- #
rule sim_mixed_ploidy:
    output:
        baf="results/hmm_simulations/mix_ploidy_{p_mono}_{p_tri}/sim{rep}_m{m}_n{n}_pi{pi0}_sigma{sigma}_skew{skew}.npz",
    resources:
        time="0:10:00",
        mem_mb="1G",
    params:
        sfs=config["afs"],
        n=lambda wildcards: int(wildcards.n),
        p_mono=lambda wildcards: int(wildcards.p_mono) / 100,
        p_tri=lambda wildcards: int(wildcards.p_tri) / 100,
        m=lambda wildcards: int(wildcards.m),
        sigma=lambda wildcards: float(wildcards.sigma) / 100,
        pi0=lambda wildcards: float(wildcards.pi0) / 100,
        seed=lambda wildcards: int(wildcards.rep),
        mat_skew=lambda wildcards: float(wildcards.skew) / 100,
        mother_id=lambda wildcards: f"n{wildcards.n}_p{wildcards.p_mono}_{wildcards.p_tri}_m{wildcards.rep}",
        father_id=lambda wildcards: f"n{wildcards.n}_p{wildcards.p_mono}_{wildcards.p_tri}_m{wildcards.rep}",
        child_id=lambda wildcards: f"n{wildcards.n}_p{wildcards.p_mono}_{wildcards.p_tri}_m{wildcards.rep}",
        mixed_ploidy=True,
        alpha=1.0,
    script:
        "scripts/sim_data.py"


rule hmm_model_comparison_mixed:
    input:
        baf=rules.sim_mixed_ploidy.output.baf,
    output:
        hmm_out="results/hmm_ploidy_comp/mix_ploidy_{p_mono}_{p_tri}/sim{rep}_m{m}_n{n}_pi{pi0}_sigma{sigma}_skew{skew}.hmm.model_comp.npz",
    resources:
        time="0:10:00",
        mem_mb="2G",
    params:
        model_comp=True,
        unphased=False,
        lrr=False,
        phase_error=True,
        mother_id=lambda wildcards: f"n{wildcards.n}_p{wildcards.p_mono}_{wildcards.p_tri}_m{wildcards.rep}",
        father_id=lambda wildcards: f"n{wildcards.n}_p{wildcards.p_mono}_{wildcards.p_tri}_m{wildcards.rep}",
        child_id=lambda wildcards: f"n{wildcards.n}_p{wildcards.p_mono}_{wildcards.p_tri}_m{wildcards.rep}",
    script:
        "scripts/baf_hmm.py"


rule hmm_model_chromosomes_mixed:
    """Local rule that collapses all ploidy assignments to a single estimand for mixed-ploidy biopsies."""
    input:
        hmm_model="results/hmm_ploidy_comp/mix_ploidy_{p_mono}_{p_tri}/sim{rep}_m{m}_n{n}_pi{pi0}_sigma{sigma}_skew{skew}.hmm.model_comp.npz",
        sims="results/hmm_simulations/mix_ploidy_{p_mono}_{p_tri}/sim{rep}_m{m}_n{n}_pi{pi0}_sigma{sigma}_skew{skew}.npz",
    output:
        ploidy=temp(
            "results/hmm_ploidy_comp/mix_ploidy_{p_mono}_{p_tri}/sim{rep}_m{m}_n{n}_pi{pi0}_sigma{sigma}_skew{skew}.hmm.ploidy.tsv"
        ),
    run:
        sigma = int(wildcards.sigma) / 100
        pi0 = int(wildcards.pi0) / 100
        prop_mono = int(wildcards.p_mono) / 100
        prop_tri = int(wildcards.p_tri) / 100
        n = int(wildcards.n)
        with open(output.ploidy, "w") as out:
            out.write(
                "mother\tfather\tchild\trep\tm\tn_cells\tn_cells_disomy\tn_cells_monosomy\tn_cells_trisomy\tp_mono\tp_tri\tsigma\tpi0\tsigma_est\tpi0_est\t0\t1m\t1p\t2\t3m\t3p\tmax_cat\n"
            )
            data = np.load(input.hmm_model)
            sim_data = np.load(input.sims)
            cell_ids = sim_data["ploidies"]
            n_cells_disomy = np.sum(cell_ids == 2)
            n_cells_monosomy = np.sum(cell_ids == 1)
            n_cells_trisomy = np.sum(cell_ids == 3)
            cats = np.array(["0", "1m", "1p", "2", "3m", "3p"])
            posteriors = [data[x] for x in cats]
            max_cat_full = cats[np.argmax(posteriors)]
            out.write(
                f"{data['mother_id']}\t{data['father_id']}\t{data['child_id']}\t{wildcards.rep}\t{wildcards.m}\t{n}\t{n_cells_disomy}\t{n_cells_monosomy}\t{n_cells_trisomy}\t{prop_mono}\t{prop_tri}\t{sigma}\t{pi0}\t{data['sigma_est']}\t{data['pi0_est']}\t{data['0']}\t{data['1m']}\t{data['1p']}\t{data['2']}\t{data['3m']}\t{data['3p']}\t{max_cat_full}\n"
            )


rule collect_mixed_ploidy:
    input:
        mixed_ploidy_tsvs_mono=expand(
            "results/hmm_ploidy_comp/mix_ploidy_{p_mono}_{p_tri}/sim{rep}_m{m}_n{n}_pi{pi0}_sigma{sigma}_skew{skew}.hmm.ploidy.tsv",
            p_mono=config["hmm_sims"]["mixed_sims"]["p_mono"],
            p_tri=0,
            rep=range(1, config["hmm_sims"]["mixed_sims"]["reps"] + 1),
            m=config["hmm_sims"]["mixed_sims"]["m"],
            pi0=config["hmm_sims"]["mixed_sims"]["pi0"],
            sigma=config["hmm_sims"]["mixed_sims"]["std_dev"],
            skew=config["hmm_sims"]["mixed_sims"]["skew"],
            n=config["hmm_sims"]["mixed_sims"]["ncells"],
        ),
        mixed_ploidy_tsvs_tri=expand(
            "results/hmm_ploidy_comp/mix_ploidy_{p_mono}_{p_tri}/sim{rep}_m{m}_n{n}_pi{pi0}_sigma{sigma}_skew{skew}.hmm.ploidy.tsv",
            p_mono=0,
            p_tri=config["hmm_sims"]["mixed_sims"]["p_tri"],
            rep=range(1, config["hmm_sims"]["mixed_sims"]["reps"] + 1),
            m=config["hmm_sims"]["mixed_sims"]["m"],
            pi0=config["hmm_sims"]["mixed_sims"]["pi0"],
            sigma=config["hmm_sims"]["mixed_sims"]["std_dev"],
            skew=config["hmm_sims"]["mixed_sims"]["skew"],
            n=config["hmm_sims"]["mixed_sims"]["ncells"],
        ),
    output:
        tot_mixed_tsv="results/mixed_ploidy_sims.tsv.gz",
    run:
        dfs = []
        for p in input.mixed_ploidy_tsvs_mono:
            df = pd.read_csv(p, sep="\t")
            dfs.append(df)
        for p in input.mixed_ploidy_tsvs_tri:
            df = pd.read_csv(p, sep="\t")
            dfs.append(df)
        tot_df = pd.concat(dfs)
        tot_df.to_csv(output.tot_mixed_tsv, sep="\t", index=None)


# ------------ 2a. Estimation of Mosaic-Cell Fraction --------------- #
rule estimate_mosaic_cell_fraction:
    input:
        baf=rules.sim_mixed_ploidy.output.baf,
    output:
        mosaic_tsv="results/mosaic_est/mix_ploidy_{p_mono}_{p_tri}/sim{rep}_m{m}_n{n}_pi{pi0}_sigma{sigma}_skew{skew}.tsv",
    resources:
        time="0:10:00",
        mem_mb="2G",
    params:
        model_comp=True,
        unphased=False,
        lrr=False,
        phase_error=True,
        mother_id=lambda wildcards: f"n{wildcards.n}_p{wildcards.p_mono}_{wildcards.p_tri}_m{wildcards.rep}",
        father_id=lambda wildcards: f"n{wildcards.n}_p{wildcards.p_mono}_{wildcards.p_tri}_m{wildcards.rep}",
        child_id=lambda wildcards: f"n{wildcards.n}_p{wildcards.p_mono}_{wildcards.p_tri}_m{wildcards.rep}",
    script:
        "scripts/mosaic_est.py"


rule collect_mosaic_est:
    input:
        mixed_ploidy_tsvs_mono=expand(
            "results/mosaic_est/mix_ploidy_{p_mono}_{p_tri}/sim{rep}_m{m}_n{n}_pi{pi0}_sigma{sigma}_skew{skew}.tsv",
            p_mono=config["hmm_sims"]["mixed_sims"]["p_mono"],
            p_tri=0,
            rep=range(1, config["hmm_sims"]["mixed_sims"]["reps"] + 1),
            m=config["hmm_sims"]["mixed_sims"]["m"],
            pi0=config["hmm_sims"]["mixed_sims"]["pi0"],
            sigma=config["hmm_sims"]["mixed_sims"]["std_dev"],
            skew=config["hmm_sims"]["mixed_sims"]["skew"],
            n=config["hmm_sims"]["mixed_sims"]["ncells"],
        ),
        mixed_ploidy_tsvs_tri=expand(
            "results/mosaic_est/mix_ploidy_{p_mono}_{p_tri}/sim{rep}_m{m}_n{n}_pi{pi0}_sigma{sigma}_skew{skew}.tsv",
            p_mono=0,
            p_tri=config["hmm_sims"]["mixed_sims"]["p_tri"],
            rep=range(1, config["hmm_sims"]["mixed_sims"]["reps"] + 1),
            m=config["hmm_sims"]["mixed_sims"]["m"],
            pi0=config["hmm_sims"]["mixed_sims"]["pi0"],
            sigma=config["hmm_sims"]["mixed_sims"]["std_dev"],
            skew=config["hmm_sims"]["mixed_sims"]["skew"],
            n=config["hmm_sims"]["mixed_sims"]["ncells"],
        ),
    output:
        tot_mixed_tsv="results/mixed_ploidy_sims.mosaic_est.tsv.gz",
    run:
        dfs = []
        for p in input.mixed_ploidy_tsvs_mono:
            df = pd.read_csv(p, sep="\t")
            dfs.append(df)
        for p in input.mixed_ploidy_tsvs_tri:
            df = pd.read_csv(p, sep="\t")
            dfs.append(df)
        tot_df = pd.concat(dfs)
        tot_df.to_csv(output.tot_mixed_tsv, sep="\t", index=None)
