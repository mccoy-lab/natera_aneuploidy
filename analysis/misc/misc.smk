#!python3

"""Snakefile for miscellaneous analyses in the Natera dataset """


import numpy as np
import pandas as pd
from tqdm import tqdm
import pickle
import gzip

# Create the VCF data dictionary for each chromosome ...
vcf_dict_natera_parents = {}
vcf_dict_1kg_phase3 = {}
chroms = [f"chr{i}" for i in range(1, 23)]
for c in chroms:
    vcf_dict_natera_parents[
        c
    ] = f"/data/rmccoy22/natera_spectrum/genotypes/opticall_parents_100423/genotypes/eagle_phased_hg38/natera_parents.b38.{c}.vcf.gz"
    vcf_dict_1kg_phase3[
        c
    ] = f"/scratch4/rmccoy22/sharedData/populationDatasets/1KGP_phase3/GRCh38_phased_vcfs/ALL.{c}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"


localrules:
    all,


rule all:
    input:
        "results/pca/kg_phase3.natera_merged.grch38.autosomes.phased.snvs.eigenvec",
        "results/ivan_validation_data/ivan_family_embryos.tsv",
        "results/ivan_validation_data/karyohmm_validation_ivan.tsv",


# ----------- 1. PCA  of Parental Genotypes along with 1KG samples ---------- #
rule merge_1kg_phase3_natera:
    input:
        vcfgz_natera=lambda wildcards: vcf_dict_natera_parents[wildcards.chrom],
        vcfgz_1kg_phase3=lambda wildcards: vcf_dict_1kg_phase3[wildcards.chrom],
        chr_rename="data/chr_rename.txt",
    output:
        vcf_rename=temp("results/pca/kg_phase3.grch38.{chrom}.phased.snvs.vcf.gz"),
        vcf_rename_tbi=temp(
            "results/pca/kg_phase3.grch38.{chrom}.phased.snvs.vcf.gz.tbi"
        ),
        vcf_merged="results/pca/kg_phase3.natera_merged.grch38.{chrom}.phased.snvs.vcf.gz",
        vcf_merged_tbi="results/pca/kg_phase3.natera_merged.grch38.{chrom}.phased.snvs.vcf.gz.tbi",
    threads: 8
    wildcard_constraints:
        chrom="|".join(chroms),
    resources:
        time="1:00:00",
        mem_mb="4G",
    shell:
        """
        bcftools annotate --rename-chrs {input.chr_rename} {input.vcfgz_1kg_phase3} | bcftools view -v snps -c 5 -m 2 -M 2 --threads {threads} | bgzip -@{threads} > {output.vcf_rename}   
        tabix -f {output.vcf_rename}
        bcftools merge {output.vcf_rename} {input.vcfgz_natera} --threads {threads} | bcftools view -v snps -m 2 -M 2 -i 'F_MISSING < 0.01' --threads {threads} | bgzip -@{threads} > {output.vcf_merged}
        tabix -f {output.vcf_merged}
        """


rule concat_autosomes:
    input:
        merged_vcfs=expand(
            "results/pca/kg_phase3.natera_merged.grch38.{chrom}.phased.snvs.vcf.gz",
            chrom=chroms,
        ),
    output:
        concat_vcf="results/pca/kg_phase3.natera_merged.grch38.autosomes.phased.snvs.vcf.gz",
        concat_vcf_tbi="results/pca/kg_phase3.natera_merged.grch38.autosomes.phased.snvs.vcf.gz.tbi",
    threads: 12
    shell:
        "bcftools concat --threads {threads} {input.merged_vcfs} | bgzip -@{threads} > {output.concat_vcf}; tabix -f {output.concat_vcf}"


rule run_plink_pca:
    input:
        concat_vcf="results/pca/kg_phase3.natera_merged.grch38.autosomes.phased.snvs.vcf.gz",
        concat_vcf_tbi="results/pca/kg_phase3.natera_merged.grch38.autosomes.phased.snvs.vcf.gz.tbi",
    output:
        eigenvec="results/pca/kg_phase3.natera_merged.grch38.autosomes.phased.snvs.eigenvec",
        eigenval="results/pca/kg_phase3.natera_merged.grch38.autosomes.phased.snvs.eigenval",
        log="results/pca/kg_phase3.natera_merged.grch38.autosomes.phased.snvs.log",
    resources:
        time="1:00:00",
        mem_mb="4G",
    params:
        outfix="results/pca/kg_phase3.natera_merged.grch38.autosomes.phased.snvs",
        pcs=20,
    threads: 24
    shell:
        "plink2 --vcf {input.concat_vcf} --pca {params.pcs} approx --threads {threads} --out {params.outfix}"


# -------- 2. Validation experiment with data from Ivan Vogel & Eva Hoffmann -------- #


rule collect_ivan_data:
    input:
        ivan_csvs=expand(
            "/scratch4/rmccoy22/ivogel/PGD_all_families_csv/PGD_family_{n}.csv",
            n=range(1, 245),
        ),
    output:
        ivan_tot_tsv="results/ivan_validation_data/ivan_family_embryos.tsv",
    run:
        with open(output["ivan_tot_tsv"], "w+") as out:
            for i, t in tqdm(enumerate(input["ivan_csvs"])):
                df = pd.read_csv(
                    t,
                    nrows=500000,
                    dtype={
                        "Chr": str,
                        "mother_gtype": str,
                        "father_gtype": str,
                        "anonymized_sample_id": int,
                        "anonymized_snp_id": int,
                        "b_allele_freq": float,
                        "category_num": int,
                    },
                )
                ind = df["anonymized_sample_id"].unique()
                for x in ind:
                    out.write(f"{i+1}\t{x}\n")


rule run_karyohmm_ivan:
    input:
        ivan_tot_csv="results/ivan_validation_data/ivan_family_embryos.tsv",
        ivan_csv="/scratch4/rmccoy22/ivogel/PGD_all_families_csv/PGD_family_{n}.csv",
    output:
        hmm_pkl="results/ivan_validation_data/PGD_family_{n}.ind{i}.hmm_results.pkl.gz",
    resources:
        time="1:00:00",
        mem_mb="5G",
    script:
        "scripts/karyohmm_ivan_pgd_data.py"


rule process_hmm_ivan:
    input:
        ivan_tot_csv="results/ivan_validation_data/ivan_family_embryos.tsv",
        ivan_hmm_result="results/ivan_validation_data/PGD_family_{n}.ind{i}.hmm_results.pkl.gz",
    output:
        hmm_tsv="results/ivan_validation_data/PGD_family_{n}.ind{i}.hmm_results.tsv",
    resources:
        time="0:10:00",
        mem_mb="1G",
    run:
        dfs = []
        cats = np.array(["0", "1m", "1p", "2", "3m", "3p"])
        full_hmm_output = pickle.load(gzip.open(input.ivan_hmm_result, "r"))
        for c in full_hmm_output:
            for i in full_hmm_output[c]:
                cur_res = {
                    "chrom": c,
                    "indiv": i,
                    "ivan_cat": full_hmm_output[c][i]["ivan_max_cat"],
                    "pi0_est": full_hmm_output[c][i]["pi0_baf"],
                    "sigma_hat": full_hmm_output[c][i]["sigma_baf"],
                }
                cur_res.update({x: full_hmm_output[c][i][x] for x in cats})
                # Defining the maximal category from Ivan?
                df = pd.DataFrame(cur_res, index=[0])
                dfs.append(df)
        tot_df = pd.concat(dfs)
        tot_df.to_csv(output.hmm_tsv, sep="\t", na_rep="NA", index=None)


def collect_hmm_results():
    """Collect all the HMM results that we expect.

    NOTE: this allows us to parallelize across an entire family...
    """
    total_data = []
    if Path("results/ivan_validation_data/ivan_family_embryos.tsv").is_file():
        with open("results/ivan_validation_data/ivan_family_embryos.tsv", "r") as fp:
            for line in fp:
                [n, i] = line.rstrip().split()
                total_data.append(
                    f"results/ivan_validation_data/PGD_family_{n}.ind{i}.hmm_results.tsv"
                )
    return total_data


rule collect_karyohmm_ivan:
    """We have ~245 families that we are validating our embryo calls on."""
    input:
        ivan_tsv_data="results/ivan_validation_data/ivan_family_embryos.tsv",
        karyohmm_calls=collect_hmm_results(),
    output:
        res_tsv="results/ivan_validation_data/karyohmm_validation_ivan.tsv",
    resources:
        time="0:30:00",
        mem_mb="4G",
    run:
        dfs = []
        cats = np.array(["0", "1m", "1p", "2", "3m", "3p"])
        for fp in input["karyohmm_calls"]:
            df = pd.read_csv(fp, sep="\t")
            dfs.append(df)
        tot_df = pd.concat(dfs)
        tot_df.to_csv(output["res_tsv"], sep="\t", index=None)
