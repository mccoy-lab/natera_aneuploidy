#!python3

# Usage: conda activate natera-aneuploidy-gwas \ ml snakemake
# Usage on dev node: nohup snakemake -p --cores 48 -j 12 --snakefile gwas.smk > nohup_date.out 2>&1 &
# Usage on rockfish: nohup snakemake -p --snakefile gwas.smk -j 200 --profile ~/code/rockfish_smk_profile/ &
# Optional: add -n for a dry run
# Executed from /scratch16/rmccoy22/scarios1/natera_aneuploidy/analysis/gwas/

# -------- Parameters ------- #
configfile: "config.yaml"

# Dictionary of number of files to chunk each vcf into in `split`
chunks_dict = {
    "chr1": 60,
    "chr2": 60,
    "chr3": 60,
    "chr4": 60,
    "chr5": 60,
    "chr6": 60,
    "chr7": 60,
    "chr8": 60,
    "chr9": 60,
    "chr10": 60,
    "chr11": 60,
    "chr12": 60,
    "chr13": 60,
    "chr14": 60,
    "chr15": 50,
    "chr16": 40,
    "chr17": 40,
    "chr18": 40,
    "chr19": 30,
    "chr20": 30,
    "chr21": 20,
    "chr22": 20,
    "chr23": 60,
}

# Parameters pipeline will run on
phenotypes = [
    "embryo_count",
    "maternal_age",
    "maternal_meiotic_aneuploidy",
    "haploidy",
    "triploidy",
]
parents = ["mother", "father"]
dataset_type = ["discovery", "test"]
chroms = range(1, 24)

# shell.prefix("set -o pipefail; ")


# -------- Rules section -------- #
rule all:
    input:
        expand(
            "results/gwas/summary_stats/gwas_{phenotype}_by_{parent}_{dataset_type}_total.tsv.gz",
            phenotype="embryo_count",
            parent="mother",
            dataset_type=dataset_type,
        ),
# rule all:
#     input:
#         expand(
#             "results/gwas/summary_stats/gwas_{phenotype}_by_{parent}_{dataset_type}_total.tsv.gz",
#             phenotype="embryo_count",
#             parent="mother",
#             dataset_type="discovery",
#         ),


# -------- 0. Preprocess genetic data -------- #
rule vcf2pgen:
    """Convert imputed VCF files to PGEN format for use in KING and PCA"""
    input:
        input_vcf="/data/rmccoy22/natera_spectrum/genotypes/imputed_parents_101823_cpra/spectrum_imputed_chr{chrom}_rehead_filter_cpra.vcf.gz",
    output:
        pgen=temp("results/gwas/intermediate_files/spectrum_imputed_chr{chrom}.pgen"),
        psam=temp("results/gwas/intermediate_files/spectrum_imputed_chr{chrom}.psam"),
        pvar=temp("results/gwas/intermediate_files/spectrum_imputed_chr{chrom}.pvar"),
        log=temp("results/gwas/intermediate_files/spectrum_imputed_chr{chrom}.log"),
    threads: 24
    wildcard_constraints:
        chrom = "|".join(map(str, range(1, 23))),
    resources:
    	mem_mb="100G",
    params:
        outfix=lambda wildcards: f"results/gwas/intermediate_files/spectrum_imputed_chr{wildcards.chrom}",
    shell:
        "plink2 --vcf {input.input_vcf} dosage=DS --double-id --maf 0.005 --threads {threads} --make-pgen --out {params.outfix}"


rule merge_full_pgen:
    """Merge PGEN file from each chromosome into a consolidated PGEN file"""
     input:
        expand("results/gwas/intermediate_files/spectrum_imputed_chr{chrom}.pgen", chrom=range(1, 23)),
        expand("results/gwas/intermediate_files/spectrum_imputed_chr{chrom}.psam", chrom=range(1, 23)),
        expand("results/gwas/intermediate_files/spectrum_imputed_chr{chrom}.pvar", chrom=range(1, 23)),
    output:
        tmp_merge_file=temp("results/gwas/intermediate_files/merged_pgen.txt"),
        pgen="results/gwas/intermediate_files/merged_imputed.pgen",
        psam="results/gwas/intermediate_files/merged_imputed.psam",
        pvar="results/gwas/intermediate_files/merged_imputed.pvar",
    params:
        outfix="results/gwas/intermediate_files/merged_imputed"
    resources:
        time="6:00:00",
        mem_mb="100G",
    threads: 24
    shell:
        """
        for i in $(seq 1 22); do echo \"results/gwas/intermediate_files/spectrum_imputed_chr${{i}}\" ; done > {output.tmp_merge_file}
        plink2 --pmerge-list {output.tmp_merge_file} --maf 0.005 --threads {threads} --make-pgen --out {params.outfix}
        """

# -------- 1. Generate files necessary for GWAS -------- #
rule compute_pcs:
    """Compute PCA using genotype data"""
    input:
        pgen=rules.merge_full_pgen.output.pgen,
        psam=rules.merge_full_pgen.output.psam,
        pvar=rules.merge_full_pgen.output.pvar,
    output:
        keep_variants="results/gwas/intermediate_files/merged_imputed.prune.in",
        remove_variants=temp("results/gwas/intermediate_files/merged_imputed.prune.out"),
        evecs="results/gwas/intermediate_files/merged_imputed.eigenvec",
        evals="results/gwas/intermediate_files/merged_imputed.eigenval",
    resources:
        time="6:00:00",
        mem_mb="50G",
    params:
        npcs=20,
        outfix="results/gwas/intermediate_files/merged_imputed"
    threads: 24
    shell:
        """
        plink2 --pgen {input.pgen} --psam {input.psam} --pvar {input.pvar} --threads {threads} --maf 0.01 --indep-pairwise 200 25 0.4 --out {params.outfix}
        plink2 --pgen {input.pgen} --psam {input.psam} --pvar {input.pvar} --extract {output.keep_variants} --pca {params.npcs} approx --threads {threads} --out {params.outfix}
        """


rule king_related_individuals:
    """Isolate related individuals (up to second degree) to exclude from GWAS"""
    input:
        pgen=rules.merge_full_pgen.output.pgen,
        psam=rules.merge_full_pgen.output.psam,
        pvar=rules.merge_full_pgen.output.pvar,
    output:
        king_include=temp("results/gwas/intermediate_files/king_result.king.cutoff.in.id"),
        king_exclude="results/gwas/intermediate_files/king_result.king.cutoff.out.id",
    resources:
        time="6:00:00",
        mem_mb="50G",
    threads: 24
    params:
        king_thresh=0.125,
        outfix="results/gwas/intermediate_files/king_result"
    shell:
        "plink2 --pgen {input.pgen} --psam {input.psam} --pvar {input.pvar} --threads {threads} --maf 0.01 --king-cutoff {params.king_thresh} --out {params.outfix}"


rule discovery_validate_split:
    """Split families into discovery/validation sets for use in GWAS"""
    input:
        discovery_validate_R="scripts/discovery_test_split.R",
        metadata=config['metadata'],
        fam_file=config['opticall'],
        king_to_remove=rules.king_related_individuals.output.king_exclude,
    output:
        metadata_weighted_ages="results/gwas/intermediate_files/spectrum_metadata_weighted_ages.tsv",
        discovery_validate_maternal="results/gwas/intermediate_files/discover_validate_split_mother.txt",
        discovery_validate_paternal="results/gwas/intermediate_files/discover_validate_split_father.txt",
    shell:
        "Rscript --vanilla {input.discovery_validate_R} {input.metadata} {input.fam_file} {input.king_to_remove} {output.metadata_weighted_ages} {output.discovery_validate_maternal} {output.discovery_validate_paternal}"


# -------- 2. Subset genetic data for each chromosome to decrease computation time -------- #
rule get_chrom_pos:
    input:
        input_vcf="/data/rmccoy22/natera_spectrum/genotypes/imputed_parents_101823_cpra/spectrum_imputed_chr{chrom}_rehead_filter_cpra.vcf.gz",
    output:
        chrom_mapfile="results/gwas/subsets/spectrum_imputed_chr{chrom}_chrom_pos.txt",
    resources:
        mem_mb="2G",
    threads: 1
    wildcard_constraints:
        chrom = "|".join(map(str, range(1, 23))),
    shell:
        "bcftools query -f'%CHROM\t%POS\n' {input.input_vcf} > {output.chrom_mapfile}"


rule make_vcf_regions:
    input:
        chrom_mapfile=rules.get_chrom_pos.output.chrom_mapfile,
        get_regions="scripts/gwas/get_regions.py",
    output:
        regions_file="results/gwas/subsets/spectrum_imputed_chr{chrom}_regions.txt",
    resources:
        mem_mb="2G",
    params:
        nchunks=lambda wildcards: chunks_dict[f"chr{wildcards.chrom}"],
    threads: 1
    wildcard_constraints:
        chrom = "|".join(map(str, range(1, 23))),
    shell:
        "python3 {input.get_regions} {params.nchunks} {input.chrom_mapfile} {output.regions_file}"


rule bed_split_vcf:
    input:
        regions_file=rules.make_vcf_regions.output.regions_file,
        input_vcf="/data/rmccoy22/natera_spectrum/genotypes/imputed_parents_101823_cpra/spectrum_imputed_chr{chrom}_rehead_filter_cpra.vcf.gz",
    output:
        bcf="results/gwas/subsets/spectrum_imputed_chr{chrom}_rehead_filter_cpra_{chunk}.bcf",
        bed="results/gwas/subsets/spectrum_imputed_chr{chrom}_rehead_filter_cpra_{chunk}.bed",
        bim="results/gwas/subsets/spectrum_imputed_chr{chrom}_rehead_filter_cpra_{chunk}.bim",
        fam="results/gwas/subsets/spectrum_imputed_chr{chrom}_rehead_filter_cpra_{chunk}.fam",
        log="results/gwas/subsets/spectrum_imputed_chr{chrom}_rehead_filter_cpra_{chunk}.log",
    resources:
        mem_mb="3G",
    params:
        nchunks=lambda wildcards: chunks_dict[f"chr{wildcards.chrom}"],
        outfix="results/gwas/subsets/spectrum_imputed_chr{chrom}_rehead_filter_cpra_{chunk}",
    threads: 16
    wildcard_constraints:
        chrom = "|".join(map(str, range(1, 23))),
    shell:
        """
        region=$(awk -v n={wildcards.chunk} "NR==n+1 {{print}}" {input.regions_file})
        bcftools view -r $region -Ob {input.input_vcf} > {output.bcf}
        plink2 --memory 9000 --bcf {output.bcf} --double-id --allow-extra-chr --make-bed --out {params.outfix}
        """


# -------- 3. Create aneuploidy phenotypes -------- #
rule generate_aneuploidy_phenotypes:
    """Make file for each aneuploidy phenotype"""
    input:
        rscript="scripts/phenotypes/generate_phenotype_files.R",
        ploidy_calls=config['ploidy_calls'],
        segmental_calls=config['segmental_calls'],
        metadata=config['metadata'],
    output:
        phenotype_file="results/phenotypes/{phenotype}_by_{parent}.csv",
    wildcard_constraints:
        phenotype="embryo_count|maternal_age|maternal_meiotic_aneuploidy|haploidy|triploidy",
        parent="mother|father",
    resources:
        time="0:30:00",
        mem_mb="10G",
    params:
        filter_day_5="TRUE",
        bayes_factor_cutoff=2,
        nullisomy_threshold=5,
        min_prob=0.9,
        max_meiotic=3,
        min_ploidy=15,
    shell:
        """
        ml gcc r/4.0.2
        Rscript --vanilla {input.rscript} {input.ploidy_calls} {input.segmental_calls} {wildcards.parent} {input.metadata} {wildcards.phenotype} {params.filter_day_5} {params.bayes_factor_cutoff} {params.nullisomy_threshold} {params.min_prob} {params.max_meiotic} {params.min_ploidy} {output.phenotype_file}
        """


# -------- 4. Execute GWAS and concatenate files -------- #
rule run_gwas_subset:
    """Run GWAS for each set of parameters, using the subsetted bed files"""
    input:
        gwas_rscript="scripts/gwas/gwas_all.R",
        metadata=config['metadata'],
        bed=rules.bed_split_vcf.output.bed,
        discovery_test="results/gwas/intermediate_files/discover_validate_split_{parent}.txt",
        parental_pcs=rules.compute_pcs.output.evecs,
        phenotype_file=rules.generate_aneuploidy_phenotypes.output.phenotype_file,
        bim=rules.bed_split_vcf.output.bim,
    output:
        gwas_output=temp("results/gwas/summary_stats/subset_gwas_{phenotype}_by_{parent}_{dataset_type}_{chrom}_{chunk}.tsv"),
    threads: 16
    resources:
        time="0:30:00",
        mem_mb="10G",
    wildcard_constraints:
        dataset_type="discovery|test",
        phenotype="embryo_count|maternal_age|maternal_meiotic_aneuploidy|haploidy|triploidy",
        parent="mother|father",
        chrom = "|".join(map(str, range(1, 23))),
    shell:
        """
        ml gcc r/4.0.2
        Rscript --vanilla {input.gwas_rscript} {input.metadata} {input.bed} {input.discovery_test} {input.parental_pcs} {input.phenotype_file} {input.bim} {wildcards.dataset_type} {wildcards.phenotype} {wildcards.parent} {threads} {output.gwas_output}
        """


rule merge_subsets:
    """Create single file for GWAS for each chromosome, merging all subsets"""
    input:
        lambda wildcards: expand(
            "results/gwas/summary_stats/subset_gwas_{{phenotype}}_by_{{parent}}_{{dataset_type}}_{{chrom}}_{chunk}.tsv",
            phenotype=wildcards.phenotype,
            parent=wildcards.parent,
            dataset_type=wildcards.dataset_type,
            chrom=wildcards.chrom,
            chunk=range(chunks_dict.get(f"chr{wildcards.chrom}", 0)),
        ),
    output:
        gwas_output=temp("results/gwas/summary_stats/gwas_{phenotype}_by_{parent}_{dataset_type}_{chrom}.tsv"),
    wildcard_constraints:
        chrom = "|".join(map(str, range(1, 23))),
    shell:
        "cat {input} > {output.gwas_output}"


rule gwas_x_chrom: 
    """Compute GWAS for the whole X chromosome."""
    input:
        gwas_rscript="scripts/gwas/gwas_all.R",
        metadata=config['metadata'],
        bed="/data/rmccoy22/natera_spectrum/genotypes/imputed_parents_101823_cpra/spectrum_imputed_chr23_rehead_filter_plink_cpra.bed",
        discovery_test="results/gwas/intermediate_files/discover_validate_split_{parent}.txt",
        parental_pcs=rules.compute_pcs.output.evecs,
        phenotype_file=rules.generate_aneuploidy_phenotypes.output.phenotype_file,
        bim="/data/rmccoy22/natera_spectrum/genotypes/imputed_parents_101823_cpra/spectrum_imputed_chr23_rehead_filter_plink_cpra.bim",
    output:
        gwas_output=temp("results/gwas/summary_stats/gwas_{phenotype}_by_{parent}_{dataset_type}_23.tsv"),
    threads: 16
    resources:
        time="6:00:00",
        mem_mb="100G",
    wildcard_constraints:
        dataset_type="discovery|test",
        phenotype="embryo_count|maternal_age|maternal_meiotic_aneuploidy|haploidy|triploidy",
        parent="mother|father",
    shell:
        """
        ml gcc r/4.0.2
        Rscript --vanilla {input.gwas_rscript} {input.metadata} {input.bed} {input.discovery_test} {input.parental_pcs} {input.phenotype_file} {input.bim} {wildcards.dataset_type} {wildcards.phenotype} {wildcards.parent} {threads} {output.gwas_output}
        """


rule merge_chroms:
    """Create single file for each phenotype/parent/dataset, merging all chromosomes"""
    input:
        expand(
            "results/gwas/summary_stats/gwas_{{phenotype}}_by_{{parent}}_{{dataset_type}}_{chrom}.tsv",
            chrom=range(1, 24),
        ),
    output:
        merged_file="results/gwas/summary_stats/gwas_{phenotype}_by_{parent}_{dataset_type}_total.tsv.gz",
    shell:
        "cat {input} | gzip > {output.merged_file}"
