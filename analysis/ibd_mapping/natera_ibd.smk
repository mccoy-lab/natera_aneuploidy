#!python3

# ------- Workflow for looking at IBD-based mapping in the Natera dataset ------- #

import numpy as np
import pandas as pd

import pickle, gzip
from pathlib import Path
from io import StringIO


# ---- Parameters for inference in Natera Data ---- #
metadata_file = "../../data/spectrum_metadata_merged.csv"

# Create the VCF data dictionary for each chromosome ...
vcf_dict = {}
chroms = [f"chr{i}" for i in range(1, 23)]
for i, c in enumerate(range(1, 23)):
    vcf_dict[chroms[i]] = (
        f"/data/rmccoy22/natera_spectrum/genotypes/opticall_parents_100423/genotypes/eagle_phased_hg38/natera_parents.b38.chr{c}.vcf.gz"
    )


# Read in the aggregate metadata file
# meta_df = pd.read_csv(metadata_file)


rule all:
    input:
        expand(
            "results/natera_parents.{chrom}.refined_endpoints.ibd.gz",
            chrom=[f"chr{i}" for i in range(1, 23)],
        ),


rule download_hapIBD:
    output:
        "bin/hap-ibd.jar",
        "bin/ibd-ends.jar",
    shell:
        """
        mkdir -p bin/
        wget https://faculty.washington.edu/browning/hap-ibd.jar -O bin/hap-ibd.jar
        wget https://faculty.washington.edu/browning/ibd-ends.jar -O bin/ibd-ends.jar
        """


rule download_genmaps:
    output:
        expand(
            "results/genmaps/genmap.{chrom}.GRCh38.corrected.map",
            chrom=[f"chr{x}" for x in range(1, 23)],
        ),
    shell:
        """
        mkdir -p results/genmaps/ ; cd results/genmaps/
        wget https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip
        unzip plink.GRCh38.map.zip
        for i in {{1..22}}; do awk \'{{$1="chr"$1; OFS=\"\t\"; print $1,$2,$3,$4}}\' plink.chr${{i}}.GRCh38.map > genmap.chr${{i}}.GRCh38.corrected.map; done
        """


rule call_ibd:
    """Call IBD segments across all pairwise samples."""
    input:
        vcf=lambda wildcards: vcf_dict[wildcards.chrom],
        genmap="results/genmaps/genmap.{chrom}.GRCh38.corrected.map",
    output:
        ibd="results/natera_parents.{chrom}.ibd.gz",
        hbd="results/natera_parents.{chrom}.hbd.gz",
    params:
        outfix=lambda wildcards: f"results/natera_parents.{wildcards.chrom}",
    threads: 8
    resources:
        mem_mb="16G",
    shell:
        "java -Xmx16g -jar bin/hap-ibd.jar gt={input.vcf} map={input.genmap} out={params.outfix} nthreads={threads}"


rule refine_ends:
    """Refine the ending placements of IBD-segments."""
    input:
        vcf=lambda wildcards: vcf_dict[wildcards.chrom],
        genmap="results/genmaps/genmap.{chrom}.GRCh38.corrected.map",
        ibd="results/natera_parents.{chrom}.ibd.gz",
    output:
        ibd="results/natera_parents.{chrom}.refined_endpoints.ibd.gz",
    params:
        seed=42,
        outfix=lambda wildcards: f"results/natera_parents.{wildcards.chrom}.refined_endpoints",
    threads: 8
    shell:
        "java -jar bin/ibd-ends.jar gt={input.vcf} ibd={input.ibd} map={input.genmap} out={params.outfix} seed={params.seed} nthreads={threads}"
