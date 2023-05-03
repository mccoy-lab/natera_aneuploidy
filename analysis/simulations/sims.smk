#!python3
# The main entry point of your workflow.

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
    TARGETS.append("results/total_hmm_ploidy.tsv")
if config["hmm_sims"]["mixed_ploidy"]:
    TARGETS.append("results/mixed_ploidy_sims.tsv")
if config["hmm_sims"]["fpr_sims"]:
    TARGETS.append("results/fpr_sims_hmm_ploidy.tsv")


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
        baf=temp(
            "results/hmm_simulations/ploidy{k}/sim{rep}_m{m}_pi{pi0}_sigma{sigma}_skew{skew}.a{a}.npz"
        ),
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
        seed=lambda wildcards: int(wildcards.rep),
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
        hmm_out=temp(
            "results/hmm_ploidy_comp/ploidy{k}/sim{rep}_m{m}_pi{pi0}_sigma{sigma}_skew{skew}.a{a}.phase_error{p}.lrr{lrr}.model_comp.npz"
        ),
    wildcard_constraints:
        lrr="0|1",
        p="0|1",
    resources:
        time="0:30:00",
        mem_mb="2G",
    params:
        model_comp=True,
        eps=-4,
        unphased=False,
        lrr=lambda wildcards: wildcards.lrr == "1",
        phase_error=lambda wildcards: wildcards.p == "1",
        mother_id=lambda wildcards: f"p{wildcards.p}_m{wildcards.rep}",
        father_id=lambda wildcards: f"p{wildcards.p}_m{wildcards.rep}",
        child_id=lambda wildcards: f"p{wildcards.p}_m{wildcards.rep}",
    script:
        "scripts/baf_hmm.py"


rule hmm_combine_baf_lrr:
    """Local rule that collapses all ploidy assignments to a single estimand for mixed-ploidy biopsies."""
    input:
        hmm_model="results/hmm_ploidy_comp/ploidy{k}/sim{rep}_m{m}_pi{pi0}_sigma{sigma}_skew{skew}.a{a}.phase_error{p}.lrr{lrr}.model_comp.npz",
    output:
        ploidy=temp(
            "results/hmm_ploidy_comp/ploidy{k}/sim{rep}_m{m}_pi{pi0}_sigma{sigma}_skew{skew}.a{a}.phase_error{p}.lrr{lrr}.model_comp.tsv"
        ),
    run:
        sigma = int(wildcards.sigma) / 100
        pi0 = int(wildcards.pi0) / 100
        a = int(wildcards.a) / 100
        phase_err = int(wildcards.p) == 1
        lrr = int(wildcards.lrr) == 1
        with open(output.ploidy, "w") as out:
            data = np.load(input.hmm_model)
            if lrr:
                out.write(
                    "mother\tfather\tchild\taploid\trep\tm\tlrr\ta\tsigma\tpi0\tsigma_est\tpi0_est\t0\t1m\t1p\t2\t2m\t2p\t3m\t3p\n"
                )
                out.write(
                    f"{data['mother_id']}\t{data['father_id']}\t{data['child_id']}\t{data['aploid']}\t{wildcards.rep}\t{wildcards.m}\t{lrr}\t{a}\t{sigma}\t{pi0}\t{data['sigma_est']}\t{data['pi0_est']}\t{data['0']}\t{data['1m']}\t{data['1p']}\t{data['2']}\t{data['2m']}\t{data['2p']}\t{data['3m']}\t{data['3p']}\n"
                )
            else:
                out.write(
                    "mother\tfather\tchild\taploid\trep\tm\tlrr\ta\tsigma\tpi0\tsigma_est\tpi0_est\t0\t1m\t1p\t2\t3m\t3p\n"
                )
                out.write(
                    f"{data['mother_id']}\t{data['father_id']}\t{data['child_id']}\t{data['aploid']}\t{wildcards.rep}\t{wildcards.m}\t{lrr}\t{a}\t{sigma}\t{pi0}\t{data['sigma_est']}\t{data['pi0_est']}\t{data['0']}\t{data['1m']}\t{data['1p']}\t{data['2']}\t{data['3m']}\t{data['3p']}\n"
                )


rule collect_hmm_model_baf_lrr:
    input:
        hmm_tsvs=expand(
            "results/hmm_ploidy_comp/ploidy{k}/sim{rep}_m{m}_pi{pi0}_sigma{sigma}_skew{skew}.a{a}.phase_error{p}.lrr{lrr}.model_comp.tsv",
            k=config["hmm_sims"]["simple_sims"]["ploidy"],
            rep=range(1, config["hmm_sims"]["simple_sims"]["reps"] + 1),
            m=config["hmm_sims"]["simple_sims"]["m"],
            pi0=config["hmm_sims"]["simple_sims"]["pi0"],
            sigma=config["hmm_sims"]["simple_sims"]["std_dev"],
            skew=config["hmm_sims"]["simple_sims"]["skew"],
            a=[10, 30, 50],
            p=[1],
            lrr=[0, 1],
        ),
    output:
        tot_hmm_tsv="results/total_hmm_ploidy.tsv",
    run:
        dfs = []
        for p in input.hmm_tsvs:
            df = pd.read_csv(p, sep="\t")
            dfs.append(df)
        tot_df = pd.concat(dfs)
        tot_df.to_csv(output.tot_hmm_tsv, sep="\t", index=None)


rule collect_fpr_baf_model_data:
    """Collect FPR-data simulations for evaluating posterior cutoffs."""
    input:
        hmm_tsvs=expand(
            "results/hmm_ploidy_comp/ploidy{k}/sim{rep}_m{m}_pi{pi0}_sigma{sigma}_skew{skew}.a{a}.phase_error{p}.lrr{lrr}.model_comp.tsv",
            k=config["hmm_sims"]["fpr_sims"]["ploidy"],
            rep=range(1, config["hmm_sims"]["fpr_sims"]["reps"] + 1),
            m=config["hmm_sims"]["fpr_sims"]["m"],
            pi0=config["hmm_sims"]["fpr_sims"]["pi0"],
            sigma=config["hmm_sims"]["fpr_sims"]["std_dev"],
            skew=config["hmm_sims"]["fpr_sims"]["skew"],
            a=100,
            p=1,
            lrr=0,
        ),
    output:
        tot_hmm_tsv="results/fpr_sims_hmm_ploidy.tsv",
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
        baf=temp(
            "results/hmm_simulations/mix_ploidy_{p}/sim{rep}_m{m}_n{n}_pi{pi0}_sigma{sigma}_skew{skew}.npz"
        ),
    resources:
        time="0:10:00",
        mem_mb="1G",
    params:
        sfs=config["afs"],
        n=lambda wildcards: int(wildcards.n),
        p=lambda wildcards: int(wildcards.p) / 100,
        m=lambda wildcards: int(wildcards.m),
        sigma=lambda wildcards: float(wildcards.sigma) / 100,
        pi0=lambda wildcards: float(wildcards.pi0) / 100,
        seed=lambda wildcards: int(wildcards.rep),
        mat_skew=lambda wildcards: float(wildcards.skew) / 100,
        mother_id=lambda wildcards: f"p{wildcards.p}_m{wildcards.rep}",
        father_id=lambda wildcards: f"p{wildcards.p}_m{wildcards.rep}",
        child_id=lambda wildcards: f"p{wildcards.p}_m{wildcards.rep}",
        mixed_ploidy=True,
        alpha=1.0,
    script:
        "scripts/sim_data.py"


rule hmm_model_comparison_mixed:
    input:
        baf=rules.sim_mixed_ploidy.output.baf,
    output:
        hmm_out=temp(
            "results/hmm_ploidy_comp/mix_ploidy_{p}/sim{rep}_m{m}_n{n}_pi{pi0}_sigma{sigma}_skew{skew}.hmm.model_comp.npz"
        ),
    resources:
        time="0:30:00",
        mem_mb="2G",
    params:
        model_comp=True,
        eps=-4,
        unphased=False,
        lrr=False,
        phase_error=True,
        mother_id=lambda wildcards: f"p{wildcards.p}_m{wildcards.rep}",
        father_id=lambda wildcards: f"p{wildcards.p}_m{wildcards.rep}",
        child_id=lambda wildcards: f"p{wildcards.p}_m{wildcards.rep}",
    script:
        "scripts/baf_hmm.py"


rule hmm_model_chromosomes_mixed:
    """Local rule that collapses all ploidy assignments to a single estimand for mixed-ploidy biopsies."""
    input:
        hmm_model="results/hmm_ploidy_comp/mix_ploidy_{p}/sim{rep}_m{m}_n{n}_pi{pi0}_sigma{sigma}_skew{skew}.hmm.model_comp.npz",
        sims="results/hmm_simulations/mix_ploidy_{p}/sim{rep}_m{m}_n{n}_pi{pi0}_sigma{sigma}_skew{skew}.npz",
    output:
        ploidy=temp(
            "results/hmm_ploidy_comp/mix_ploidy_{p}/sim{rep}_m{m}_n{n}_pi{pi0}_sigma{sigma}_skew{skew}.hmm.ploidy.tsv"
        ),
    run:
        sigma = int(wildcards.sigma) / 100
        pi0 = int(wildcards.pi0) / 100
        prop = int(wildcards.p) / 100
        n = int(wildcards.n)
        with open(output.ploidy, "w") as out:
            out.write(
                "mother\tfather\tchild\trep\tm\tploidies\tn\tp_aneuploid\tsigma\tpi0\tsigma_est\tpi0_est\t1m\t1p\t2\t3m\t3p\n"
            )
            data = np.load(input.hmm_model)
            sim_data = np.load(input.sims)
            ploid_str = ",".join([str(x) for x in sim_data["ploidies"]])
            out.write(
                f"{data['mother_id']}\t{data['father_id']}\t{data['child_id']}\t{wildcards.rep}\t{wildcards.m}\t{ploid_str}\t{n}\t{prop}\t{sigma}\t{pi0}\t{data['sigma_est']}\t{data['pi0_est']}\t{data['1m']}\t{data['1p']}\t{data['2']}\t{data['3m']}\t{data['3p']}\n"
            )


rule collect_mixed_ploidy:
    input:
        mixed_ploidy_tsvs=expand(
            "results/hmm_ploidy_comp/mix_ploidy_{p}/sim{rep}_m{m}_n{n}_pi{pi0}_sigma{sigma}_skew{skew}.hmm.ploidy.tsv",
            p=config["hmm_sims"]["mixed_sims"]["p_aneuploid"],
            rep=range(1, config["hmm_sims"]["mixed_sims"]["reps"] + 1),
            m=config["hmm_sims"]["mixed_sims"]["m"],
            pi0=config["hmm_sims"]["mixed_sims"]["pi0"],
            sigma=config["hmm_sims"]["mixed_sims"]["std_dev"],
            skew=config["hmm_sims"]["mixed_sims"]["skew"],
            n=config["hmm_sims"]["mixed_sims"]["ncells"],
        ),
    output:
        tot_mixed_tsv="results/mixed_ploidy_sims.tsv",
    run:
        dfs = []
        for p in input.mixed_ploidy_tsvs:
            df = pd.read_csv(p, sep="\t")
            dfs.append(df)
        tot_df = pd.concat(dfs)
        tot_df.to_csv(output.tot_mixed_tsv, sep="\t", index=None)
