#!python3

# Usage: conda activate natera-aneuploidy-gwas \ ml snakemake

# =================
# author: Sara A. Carioscia, Biology Dept., Johns Hopkins University
# email: scarios1@jhu.edu
# last update: July 16, 2025
# aim: run GWAS for single-chromosome aneuploidy with each chromosome (test for cis-effects)
# =================

# -------- Parameters ------- #
configfile: "../gwas/config.yaml"


# Dictionary of number of files to chunk each vcf into in `split`
chunks_dict = {
    "chr1": 600,
    "chr2": 600,
    "chr3": 600,
    "chr4": 600,
    "chr5": 600,
    "chr6": 600,
    "chr7": 480,
    "chr8": 480,
    "chr9": 240,
    "chr10": 480,
    "chr11": 480,
    "chr12": 480,
    "chr13": 240,
    "chr14": 240,
    "chr15": 200,
    "chr16": 160,
    "chr17": 160,
    "chr18": 160,
    "chr19": 120,
    "chr20": 120,
    "chr21": 40,
    "chr22": 80,
    "chr23": 240,
}

# Phenotype for each chromosome 
# single_chr_phenotypes = [
#     "chr1_aneuploidy",
#     "chr2_aneuploidy",
#     "chr3_aneuploidy",
#     "chr4_aneuploidy",
#     "chr5_aneuploidy",
#     "chr6_aneuploidy",
#     "chr7_aneuploidy",
#     "chr8_aneuploidy",
#     "chr9_aneuploidy",
#     "chr10_aneuploidy",
#     "chr11_aneuploidy",
#     "chr12_aneuploidy",
#     "chr13_aneuploidy",
#     "chr14_aneuploidy",
#     "chr15_aneuploidy",
#     "chr16_aneuploidy",
#     "chr17_aneuploidy",
#     "chr18_aneuploidy",
#     "chr19_aneuploidy",
#     "chr20_aneuploidy",
#     "chr21_aneuploidy",
#     "chr22_aneuploidy",
#     "chr23_aneuploidy",
# ]

phenotype = [f"chr{i}_aneuploidy" for i in range(1, 24)]
phenotype_chrom_map = {f"chr{i}_aneuploidy": str(i) for i in range(1, 24)}


parents = "mother"
#dataset_type = ["discovery", "test"]
dataset_type = "discovery"
chroms = range(1,24)
population = config["population"]

# functions to determine chromosome of interest, bed and bim file to use for GWAS
def chr_number_from_phenotype(phenotype):
    return phenotype.split("_")[0].replace("chr", "")

def bed_input(wildcards):
    chrom = chr_number_from_phenotype(wildcards.phenotype)
    return f"../gwas/results/gwas/subsets/spectrum_imputed_chr{chrom}_rehead_filter_cpra_{wildcards.chunk}.bed"

def bim_input(wildcards):
    chrom = chr_number_from_phenotype(wildcards.phenotype)
    return f"../gwas/results/gwas/subsets/spectrum_imputed_chr{chrom}_rehead_filter_cpra_{wildcards.chunk}.bim"



# -------- 0. Define rule all -------- #
# rule all_chromosomes_gwas_complete:
#     input:
#         expand(
#             "results/gwas/summary_stats/{population}/lmm_gwas_{phenotype}_by_{parent}_{dataset_type}_{population}_total.tsv",
#             population=config["population"],
#             parent="mother",
#             dataset_type=dataset_type,
#             phenotype=phenotype,
#             chrom=[phenotype_chrom_map[p] for p in phenotype],
#         )


# -------- 0. Define rule all -------- #
rule all_chromosomes_gwas_complete:
    input:
        expand(
            "results/gwas/summary_stats/{population}/lmm_gwas_{phenotype}_by_{parent}_{dataset_type}_{population}_total.tsv",
            population=config["population"],
            parent="mother",
            dataset_type=dataset_type,
            phenotype=phenotype,
            chrom=[phenotype_chrom_map[p] for p in phenotype],
        )

# -------- 1. Create phenotypes for each single chromosome aneuploidy -------- #
rule generate_phenotypes:
    """Make file for each phenotype"""
    input:
        rscript="../gwas/scripts/phenotypes/generate_phenotype_files.R",
        ploidy_calls=config["ploidy_calls"],
        segmental_calls=config["segmental_calls"],
        metadata=config["metadata"],
    output:
        phenotype_file="results/phenotypes/{phenotype}_by_{parent}.csv",
    wildcard_constraints:
        phenotype="|".join([f"chr{i}_aneuploidy" for i in range(1, 24)]),
    resources:
        time="0:10:00",
        mem_mb="10G",
    params:
        filter_day_5="TRUE",
        bayes_factor_cutoff=2,
        nullisomy_threshold=5,
        min_prob=0.9,
        max_meiotic=5,
        min_ploidy=15,
        filter_mosaics="TRUE",
    shell:
        """
        ml gcc r/4.3.0
        Rscript --vanilla {input.rscript} {input.ploidy_calls} {input.segmental_calls} {wildcards.parent} {input.metadata} {wildcards.phenotype} {params.filter_day_5} {params.bayes_factor_cutoff} {params.nullisomy_threshold} {params.min_prob} {params.max_meiotic} {params.min_ploidy} {params.filter_mosaics} {output.phenotype_file}
        """


# -------- 2. Run chr-wide association study for each single chromosome aneuploidy -------- #
rule run_gwas_lmm_autosome_subset:
    """Run GWAS LMM for each phenotype on each chunk of the correct chromosome"""
    input:
        gwas_rscript="../gwas/scripts/gwas/gwas_lmm_chunks.R",
        metadata=config["metadata"],
        discovery_test="../gwas/results/gwas/intermediate_files/discover_validate_split_{parent}.txt",
        parental_pcs="../gwas/results/gwas/intermediate_files/merged_imputed.eigenvec",
        phenotype_file="results/phenotypes/{phenotype}_by_{parent}.csv",
        assigned_pops=config["assigned_pops"],
        bed=bed_input,
        bim=bim_input,
    output:
        gwas_output=temp(
            "results/gwas/summary_stats/{population}/lmm_subset_gwas_{phenotype}_by_{parent}_{dataset_type}_{population}_{chunk}.tsv"
        ),
    threads: 16
    resources:
        time="0:50:00",
        mem_mb="60G",
    params:
        population=config["population"],
    wildcard_constraints:
        dataset_type="discovery|test",
        parent="mother",
        phenotype="|".join([f"chr{i}_aneuploidy" for i in range(1, 24)]),
    shell:
        """
        ml gcc r/4.3.0
        Rscript --vanilla {input.gwas_rscript} {input.metadata} {input.bed} {input.discovery_test} {input.parental_pcs} {input.phenotype_file} {input.bim} {wildcards.dataset_type} {input.assigned_pops} {params.population} {wildcards.phenotype} {wildcards.parent} {threads} {output.gwas_output}
        """


# -------- 3. Merge each chromosome's GWAS across chunks -------- #
rule merge_lmm_subsets:
    """Create single file for GWAS LMM for each chromosome, merging all subsets"""
    input:
        lambda wildcards: expand(
            "results/gwas/summary_stats/{{population}}/lmm_subset_gwas_{{phenotype}}_by_{{parent}}_{{dataset_type}}_{{population}}_{chunk}.tsv",
            population=config["population"],
            phenotype=wildcards.phenotype,
            parent=wildcards.parent,
            dataset_type=wildcards.dataset_type,
            chunk=range(chunks_dict.get(f"chr{chr_number_from_phenotype(wildcards.phenotype)}", 0)),
        ),
    output:
        gwas_output="results/gwas/summary_stats/{population}/lmm_gwas_{phenotype}_by_{parent}_{dataset_type}_{population}_total.tsv",
    wildcard_constraints:
        chrom="|".join(map(str, range(1, 24))),
    shell:
        "cat {input} > {output.gwas_output}"