#!python3


configfile: "config.yaml"


rule all:
    input:
        expand(
            "results/singer_trees/{rsid}/{rsid}.{subsample}.singer_0.trees",
            rsid=config["focal_variants"].keys(),
            subsample=config["subsamples"].keys(),
        ),
        expand(
            "results/betascan/{rsid}/{population}.bscan.tsv",
            rsid=config["focal_variants"].keys(),
            population=config["betascan_stats"].keys(),
        ),


rule create_locus_specific_vcf:
    input:
        vcf=lambda wildcards: config["vcfs"][
            config["focal_variants"][wildcards.rsid]["chrom"]
        ],
        subsamples=lambda wildcards: config["subsamples"][wildcards.subsample],
    output:
        region_vcf="results/singer_trees/{rsid}/{rsid}.{subsample}.regional.vcf",
    params:
        start=lambda wildcards: int(config["focal_variants"][wildcards.rsid]["pos"])
        - 2000000,
        end=lambda wildcards: int(config["focal_variants"][wildcards.rsid]["pos"])
        + 2000000,
        chrom=lambda wildcards: config["focal_variants"][wildcards.rsid]["chrom"],
    shell:
        """
        bcftools view -r {params.chrom}:{params.start}-{params.end} -S {input.subsamples} --force-samples {input.vcf} | bcftools view -c 1 > {output.region_vcf}
        """


rule run_singer:
    """Run SINGER to generate posterior samplings of regional ARGs.
    
    NOTE: should estimate Ne from pi in the region as well.
    """
    input:
        convert2tskit=config["programs"]["convert2tskit"],
        singer=config["programs"]["singer"],
        vcf=rules.create_locus_specific_vcf.output.region_vcf,
    output:
        expand(
            "results/singer_trees/{{rsid}}/{{rsid}}.{{subsample}}.singer_{x}.trees",
            x=range(20),
        ),
    params:
        start=lambda wildcards: int(config["focal_variants"][wildcards.rsid]["pos"])
        - 2000000,
        end=lambda wildcards: int(config["focal_variants"][wildcards.rsid]["pos"])
        + 2000000,
        chrom=lambda wildcards: config["focal_variants"][wildcards.rsid]["chrom"],
        mu=1.25e-8,
        thin=20,
        endx=20 - 1,
        n=100,
        vcf_fix=lambda wildcards: f"results/singer_trees/{wildcards.rsid}/{wildcards.rsid}.{wildcards.subsample}.regional",
        outfix=lambda wildcards: f"results/singer_trees/{wildcards.rsid}/{wildcards.rsid}.{wildcards.subsample}.singer",
    shell:
        """
        {input.singer} -N 2e4 -m {params.mu} -vcf {params.vcf_fix} -output {params.outfix} -n {params.n} -thin {params.thin} -seed 42 -start {params.start} -end {params.end}
        {input.convert2tskit} -input {params.outfix} -output {params.outfix} -start 0 -end {params.endx}
        """


rule extract_betascan:
    """Extract the Beta-scan statistics centered around variants."""
    input:
        betascan_stats=lambda wildcards: config["betascan_stats"][wildcards.population]
        + f"/{config['focal_variants'][wildcards.rsid]['chrom']}_B2std.out",
    output:
        regional_betascan="results/betascan/{rsid}/{population}.bscan.tsv",
    params:
        start=lambda wildcards: int(config["focal_variants"][wildcards.rsid]["pos"])
        - 1e6,
        end=lambda wildcards: int(config["focal_variants"][wildcards.rsid]["pos"]) + 1e6,
    run:
        import polars as pl

        betascan_region_df = (
            pl.read_csv(input.betascan_stats, separator="\t")
            .filter(pl.col("Position") >= params.start)
            .filter(pl.col("Position") <= params.end)
        )
        betascan_region_df.write_csv(output.regional_betascan, separator="\t")
