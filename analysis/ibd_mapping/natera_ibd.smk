#!python3

# ------- Workflow for looking at IBD-based mapping in the Natera dataset ------- #
import pickle, gzip
from pathlib import Path
from io import StringIO


configfile: "config.yaml"


rule call_ibd_total:
    """Call and refine IBD segments across parents."""
    input:
        expand(
            "results/natera_parents.{chrom}.refined_endpoints.ibd.gz",
            chrom=config["vcf"].keys(),
        ),


rule call_ibd_clusters_total:
    """Estimate IBD clusters within data to use for downstream association testing."""
    input:
        "results/natera_parents.autosomes.ibd_clusters.bed",


rule download_hapIBD:
    output:
        "bin/hap-ibd.jar",
        "bin/ibd-ends.jar",
        "bin/ibd-cluster.jar",
    shell:
        """
        mkdir -p bin/
        wget https://faculty.washington.edu/browning/hap-ibd.jar -O bin/hap-ibd.jar
        wget https://faculty.washington.edu/browning/ibd-ends.jar -O bin/ibd-ends.jar
        wget https://faculty.washington.edu/browning/ibd-cluster.jar -O bin/ibd-cluster.jar
        """


rule download_genmaps:
    output:
        expand(
            "results/genmaps/genmap.{chrom}.GRCh38.corrected.map",
            chrom=config["vcf"].keys(),
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
        vcf=lambda wildcards: config["vcf"][wildcards.chrom],
        genmap="results/genmaps/genmap.{chrom}.GRCh38.corrected.map",
        hap_ibd="bin/hap-ibd.jar",
    output:
        ibd="results/natera_parents.{chrom}.ibd.gz",
        hbd="results/natera_parents.{chrom}.hbd.gz",
    params:
        outfix=lambda wildcards: f"results/natera_parents.{wildcards.chrom}",
    threads: 8
    resources:
        mem_mb="16G",
    shell:
        "java -Xmx16g -jar {input.hap_ibd} gt={input.vcf} map={input.genmap} out={params.outfix} nthreads={threads}"


rule refine_ends:
    """Refine the ending placements of IBD-segments."""
    input:
        vcf=lambda wildcards: config["vcf"][wildcards.chrom],
        genmap="results/genmaps/genmap.{chrom}.GRCh38.corrected.map",
        ibd="results/natera_parents.{chrom}.ibd.gz",
        ibd_ends="bin/ibd-ends.jar",
    output:
        ibd="results/natera_parents.{chrom}.refined_endpoints.ibd.gz",
    params:
        seed=42,
        outfix=lambda wildcards: f"results/natera_parents.{wildcards.chrom}.refined_endpoints",
    threads: 8
    resources:
        mem_mb="16G",
    shell:
        "java -Xmx16g -jar {input.ibd_ends} gt={input.vcf} ibd={input.ibd} map={input.genmap} out={params.outfix} seed={params.seed} nthreads={threads}"


rule ibd_cluster:
    """Run IBD-clustering using the IBD-cluster method of Browning & Browning 2023."""
    input:
        vcf=lambda wildcards: config["vcf"][wildcards.chrom],
        genmap="results/genmaps/genmap.{chrom}.GRCh38.corrected.map",
        ibd_cluster="bin/ibd-cluster.jar",
    output:
        ibdclust="results/natera_parents.{chrom}.ibd_cluster.ibdclust.gz",
    threads: 8
    resources:
        mem_mb="16G",
    params:
        min_maf=0.01,
        length=3,
        outfix=lambda wildcards: f"results/natera_parents.{wildcards.chrom}.ibd_cluster",
    shell:
        "java -Xmx16g -jar {input.ibd_cluster} gt={input.vcf} map={input.genmap} min-maf={params.min_maf} length={params.length} nthreads={threads} out={params.outfix}"


rule ibdclust_to_tped:
    """Convert IBD clusters to TPED format."""
    input:
        ibd_clust=rules.ibd_cluster.output.ibdclust,
    output:
        tped=temp("results/natera_parents.{chrom}.ibd_cluster.tped"),
        tfam=temp("results/natera_parents.{chrom}.ibd_cluster.tfam"),
    params:
        clique_size=5,
        max_clique_size=100,
        max_n_cliques=50,
    script:
        "scripts/ibdclust_to_tped.py"


rule tped_to_bed:
    """Convert tped to bed files."""
    input:
        tped="results/natera_parents.{chrom}.ibd_cluster.tped",
        tfam="results/natera_parents.{chrom}.ibd_cluster.tfam",
    output:
        temp("results/natera_parents.{chrom}.ibd_cluster.bed"),
        temp("results/natera_parents.{chrom}.ibd_cluster.bim"),
        temp("results/natera_parents.{chrom}.ibd_cluster.fam"),
    params:
        outfix=lambda wildcards: f"results/natera_parents.{wildcards.chrom}.ibd_cluster",
    resources:
        mem_mb="8G",
    threads: 8
    shell:
        "plink --tped {input.tped} --tfam {input.tfam} --memory 8000 --threads {threads} --make-bed --out {params.outfix}"


rule create_mergelist:
    """Merge all autosomal BED files."""
    input:
        bedfiles=expand(
            "results/natera_parents.{chrom}.ibd_cluster.bed",
            chrom=config["vcf"].keys(),
        ),
        bimfiles=expand(
            "results/natera_parents.{chrom}.ibd_cluster.bim",
            chrom=config["vcf"].keys(),
        ),
        famfiles=expand(
            "results/natera_parents.{chrom}.ibd_cluster.fam",
            chrom=config["vcf"].keys(),
        ),
    output:
        bfile_merge=temp("results/natera_parents.merge.autosomes.files.txt"),
    run:
        with open(output["bfile_merge"], "w+") as out:
            for chrom in config["vcf"].keys():
                out.write(
                    f"results/natera_parents.{chrom}.ibd_cluster.bed results/natera_parents.{chrom}.ibd_cluster.bim results/natera_parents.{chrom}.ibd_cluster.fam\n"
                )


rule merge_bfile:
    input:
        bedfiles=expand(
            "results/natera_parents.{chrom}.ibd_cluster.bed",
            chrom=config["vcf"].keys(),
        ),
        bimfiles=expand(
            "results/natera_parents.{chrom}.ibd_cluster.bim",
            chrom=config["vcf"].keys(),
        ),
        famfiles=expand(
            "results/natera_parents.{chrom}.ibd_cluster.fam",
            chrom=config["vcf"].keys(),
        ),
        bfile_merge="results/natera_parents.merge.autosomes.files.txt",
    output:
        bed="results/natera_parents.autosomes.ibd_clusters.bed",
        bim="results/natera_parents.autosomes.ibd_clusters.bim",
        fam="results/natera_parents.autosomes.ibd_clusters.fam",
    params:
        outfix="results/natera_parents.autosomes.ibd_clusters",
    threads: 12
    resources:
        mem_mb="8G",
    shell:
        "plink --bed {input.bedfiles[0]} --bim {input.bimfiles[0]} --fam {input.famfiles[0]} --threads {threads} --memory 8000 --merge-list {input.bfile_merge} --make-bed --out {params.outfix}"
