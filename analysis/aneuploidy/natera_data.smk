#!python3

import numpy as np
import pandas as pd
from tqdm import tqdm
from pathlib import Path
from io import StringIO

# ---- Parameters for inference in Natera Data ---- #
metadata_file = "../../data/spectrum_metadata_merged.csv"
alleles_file = "/data/rmccoy22/natera_spectrum/data/illumina_files/humancytosnp-12v2-1_h.update_alleles.txt"
cluster_file = "/scratch16/rmccoy22/abiddan1/natera_spectrum/Rexpected/HumanCytoSNP-12v2-1_NS550.cluster.tsv.gz"
meanr_file = "/scratch16/rmccoy22/abiddan1/natera_spectrum/Rexpected/meanR_child.txt"
strand_file = "/data/rmccoy22/natera_spectrum/data/illumina_files/humancytosnp-12v2-1_h-b37.strand"
strand_refalt = "/data/rmccoy22/natera_spectrum/data/illumina_files/humancytosnp-12v2-1_h-b37.strand.RefAlt"
cytosnp_map_v12 = (
    "/data/rmccoy22/natera_spectrum/data/illumina_files/snp_map_cyto12b_f004.txt"
)
lrrs = ["none"]

# Create the VCF data dictionary for each chromosome ...
vcf_dict = {}
chroms = [f"chr{i}" for i in range(1, 23)]
for i, c in enumerate(range(1, 23)):
    vcf_dict[
        chroms[i]
    ] = f"/data/rmccoy22/natera_spectrum/genotypes/opticall_parents_031423/genotypes/eagle_phased_hg38/natera_parents.b38.chr{c}.vcf.gz"


# Read in the aggregate metadata file
meta_df = pd.read_csv(metadata_file)


def find_child_data(
    child_id, meta_dfi=meta_df, raw_data_path="/data/rmccoy22/natera_spectrum/data/"
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
    valid_trios = []
    unique_mothers = np.unique(
        meta_df[meta_df.family_position == "mother"].array.values
    )
    for m in tqdm(unique_mothers):
        cases = np.unique(meta_df[meta_df.array == m].casefile_id.values)
        cur_df = meta_df[np.isin(meta_df.casefile_id, cases)]
        fathers = np.unique(cur_df[cur_df.family_position == "father"].array.values)
        if fathers.size > 1:
            print(f"More than one partner for {m}")
            for fat in fathers:
                cur_cases = np.unique(meta_df[meta_df.array == fat].casefile_id.values)
                cur_df = meta_df[np.isin(meta_df.casefile_id, cur_cases)]
                for c in cur_df[cur_df.family_position == "child"].array.values:
                    valid_trios.append((m, fat, c))
        elif fathers.size == 1:
            for c in cur_df[cur_df.family_position == "child"].array.values:
                valid_trios.append((m, fathers[0], c))
    parents = [line.rstrip() for line in open(sample_file, "r")]
    # Applies a set of filters here
    valid_filt_trios = []
    for m, f, c in tqdm(valid_trios):
        if (
            (m in parents)
            and (f in parents)
            and find_child_data(c, meta_df, raw_data_path)[1]
        ):
            valid_filt_trios.append((m, f, c))

    return valid_filt_trios


total_data = []
if Path("results/natera_inference/valid_trios.txt").is_file():
    with open("results/natera_inference/valid_trios.txt", "r") as fp:
        for i, line in enumerate(fp):
            [m, f, c] = line.rstrip().split()
            for l in lrrs:
                total_data.append(
                    f"results/natera_inference/{m}+{f}/{c}.{l}.total.posterior.tsv.gz"
                )

# ------- Rules Section ------- #


localrules:
    all,
    hmm_model_chromosomes,


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
        baf_npz=expand(
            "results/natera_inference/{{mother_id}}+{{father_id}}/{{child_id}}.{chrom}.bafs.npz",
            chrom=chroms,
        ),
    resources:
        time="2:00:00",
        mem_mb="5G",
    wildcard_constraints:
        chrom="|".join(chroms),
    run:
        shell(
            "mkdir -p results/natera_inference/{wildcards.mother_id}+{wildcards.father_id}/"
        )
        for v, o in zip(input.vcf_file, output.baf_npz):
            shell(
                "python3 scripts/preprocess_natera.py --child_csv {input.child_data} --cytosnp_map {input.cytosnp_map} --alleles_file {input.alleles_file} --cytosnp_cluster {input.egt_cluster}  --mother_id {wildcards.mother_id} --father_id {wildcards.father_id} --vcf_file {v} --meanr {input.meanr_file} --outfile {o}"
            )


rule hmm_model_comparison:
    """Apply the ploidy HMM to the pre-processed BAF data for this embryo."""
    input:
        baf=expand(
            "results/natera_inference/{{mother_id}}+{{father_id}}/{{child_id}}.{chrom}.bafs.npz",
            chrom=chroms,
        ),
    output:
        hmm_out=temp(
            expand(
                "results/natera_inference/{{mother_id}}+{{father_id}}/{{child_id}}.{chrom}.{{lrr}}.hmm_model.npz",
                chrom=chroms,
            )
        ),
    wildcard_constraints:
        lrr="(none|raw|norm)",
    resources:
        time="1:00:00",
        mem_mb="4G",
    params:
        unphased=False,
        eps=-6,
        lrr=lambda wildcards: f"{wildcards.lrr}",
        mother_id=lambda wildcards: f"{wildcards.mother_id}",
        father_id=lambda wildcards: f"{wildcards.father_id}",
        child_id=lambda wildcards: f"{wildcards.child_id}",
    script:
        "scripts/baf_hmm_bulk.py"


rule hmm_model_chromosomes:
    """Local rule that collapses all ploidy assignments into a single table."""
    input:
        hmm_models=expand(
            "results/natera_inference/{{mother_id}}+{{father_id}}/{{child_id}}.{chrom}.{{lrr}}.hmm_model.npz",
            chrom=chroms,
        ),
    output:
        ploidy="results/natera_inference/{mother_id}+{father_id}/{child_id}.{lrr}.total.ploidy.tsv",
    resources:
        time="0:10:00",
        mem_mb="1G",
    params:
        lrr=lambda wildcards: wildcards.lrr != "none",
    run:
        with open(output.ploidy, "w") as out:
            if not params["lrr"]:
                out.write(
                    "mother\tfather\tchild\tchrom\tsigma_baf\tpi0_baf\tpi0_lrr\tlrr_mu\tlrr_sd\t0\t1m\t1p\t2\t3m\t3p\n"
                )
                for c, x in zip(chroms, input.hmm_models):
                    data = np.load(x)
                    out.write(
                        f"{data['mother_id']}\t{data['father_id']}\t{data['child_id']}\t{c}\t{data['sigma_baf']}\t{data['pi0_baf']}\t{data['pi0_lrr']}\t{data['lrr_mu']}\t{data['lrr_sd']}\t{data['0']}\t{data['1m']}\t{data['1p']}\t{data['2']}\t{data['3m']}\t{data['3p']}\n"
                    )
            else:
                out.write(
                    "mother\tfather\tchild\tchrom\tsigma_baf\tpi0_baf\tpi0_lrr\tlrr_mu\tlrr_sd\t0\t1m\t1p\t2m\t2p\t2\t3m\t3p\n"
                )
                for c, x in zip(chroms, input.hmm_models):
                    data = np.load(x)
                    out.write(
                        f"{data['mother_id']}\t{data['father_id']}\t{data['child_id']}\t{c}\t{data['sigma_baf']}\t{data['pi0_baf']}\t{data['pi0_lrr']}\t{data['lrr_mu']}\t{data['lrr_sd']}\t{data['0']}\t{data['1m']}\t{data['1p']}\t{data['2m']}\t{data['2p']}\t{data['2']}\t{data['3m']}\t{data['3p']}\n"
                    )


rule generate_posterior_table:
    """Generates a full TSV with posterior probabilities for each embryo across ploidy states.
    
    The columns of the TSV contain the specific karyotype states like 0, 1m, 1p, 2m, 2p, 2, 3m, 3p.
    
    Each row corresponds to a specific SNP position.
    """
    input:
        baf_data=expand(
            "results/natera_inference/{{mother_id}}+{{father_id}}/{{child_id}}.{chrom}.bafs.npz",
            chrom=chroms,
        ),
        hmm_models=expand(
            "results/natera_inference/{{mother_id}}+{{father_id}}/{{child_id}}.{chrom}.{{lrr}}.hmm_model.npz",
            chrom=chroms,
        ),
        ploidy="results/natera_inference/{mother_id}+{father_id}/{child_id}.{lrr}.total.ploidy.tsv",
    output:
        posterior="results/natera_inference/{mother_id}+{father_id}/{child_id}.{lrr}.total.posterior.tsv.gz",
    resources:
        time="0:30:00",
        mem_mb="1G",
    wildcard_constraints:
        lrr="(none|raw|norm)",
    run:
        tot_dfs = []
        for c, fp, baf in tqdm(zip(chroms, input.hmm_models, input.baf_data)):
            data = np.load(fp)
            baf_data = np.load(baf)
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
