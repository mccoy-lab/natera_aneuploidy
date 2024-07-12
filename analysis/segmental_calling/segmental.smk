#!python3

import numpy as np
import pandas as pd
import polars as pl
import pickle, gzip
from tqdm import tqdm
from pathlib import Path
from io import StringIO


configfile: "config.yaml"


# Create the VCF data dictionary for each chromosome ...
vcf_dict = {}
chroms = [f"chr{i}" for i in range(1, 23)]
for i, c in enumerate(range(1, 23)):
    vcf_dict[chroms[i]] = (
        f"/data/rmccoy22/natera_spectrum/genotypes/opticall_parents_100423/genotypes/eagle_phased_hg38/natera_parents.b38.chr{c}.vcf.gz"
    )

# Read in the aggregate metadata file
meta_df = pd.read_csv(config["metadata"])


def find_child_data(
    child_id, meta_df=meta_df, raw_data_path="/data/rmccoy22/natera_spectrum/data/"
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
    assert "family_position" in meta_df.columns
    assert "casefile_id" in meta_df.columns
    assert "array" in meta_df.columns
    grouped_df = (
        meta_df.groupby(["casefile_id", "family_position"])["array"]
        .agg(lambda x: list(x))
        .reset_index()
    )
    valid_trios = []
    for case in tqdm(np.unique(grouped_df.casefile_id)):
        cur_df = grouped_df[grouped_df.casefile_id == case]
        for m in cur_df[cur_df.family_position == "mother"].array.values[0]:
            for f in cur_df[cur_df.family_position == "father"].array.values[0]:
                for c in cur_df[cur_df.family_position == "child"].array.values[0]:
                    valid_trios.append((case, m, f, c))

    valid_df = pd.DataFrame(
        valid_trios, columns=["casefile_id", "mother", "father", "child"]
    )
    parents = [line.rstrip() for line in open(sample_file, "r")]
    valid_df["parents_in"] = valid_df.mother.isin(parents) & valid_df.father.isin(
        parents
    )
    valid_df["child_found"] = [
        find_child_data(c)[1] for c in tqdm(valid_df.child.values)
    ]
    valid_filt_df = valid_df[
        valid_df.parents_in & valid_df.child_found
    ].drop_duplicates()[["mother", "father", "child"]]
    valid_filt_trios = [
        (m, f, c)
        for (m, f, c) in zip(
            valid_filt_df.mother.values,
            valid_filt_df.father.values,
            valid_filt_df.child.values,
        )
    ]
    return valid_filt_trios


# Run inference for all valid triplets (restricting to euploid chromosomes)
total_data = []
if Path("results/segmental_inference/valid_trios.txt").is_file():
    with open("results/segmental_inference/valid_trios.txt", "r") as fp:
        for i, line in enumerate(fp):
            [m, f, c] = line.rstrip().split()
            if m != "mother":
                total_data.append(
                    f"results/segmental_inference/{m}+{f}+{c}.segmental.tsv"
                )
        total_data = np.unique(total_data).tolist()


# ------- Rules Section ------- #
localrules:
    all,
    generate_parent_sample_list,
    obtain_valid_trios,


rule all:
    input:
        "results/segmental_inference/valid_trios.txt",
        total_data,
        "results/tables/segmental_calls_raw.tsv.gz",
        "results/tables/segmental_calls_filt.tsv.gz",
        "results/tables/segmental_calls_postqc.tsv.gz",
        "results/tables/segmental_calls_postqc_refined.tsv.gz",


rule generate_parent_sample_list:
    """Generate the parental samples list."""
    input:
        vcf=[vcf_dict[c] for c in chroms],
    output:
        "results/segmental_inference/parent_samples.txt",
    resources:
        time="0:30:00",
        mem_mb="1G",
    shell:
        "bcftools query -l {input.vcf} | sort | uniq > {output}"


rule obtain_valid_trios:
    """Obtain the valid trios here."""
    input:
        metadata_tsv=config["metadata"],
        parent_samples="results/segmental_inference/parent_samples.txt",
    output:
        valid_trios="results/segmental_inference/valid_trios.txt",
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
            out.write("mother\tfather\tchild\n")
            for m, f, c in valid_trios:
                out.write(f"{m}\t{f}\t{c}\n")


def define_triplets(
    mother_id,
    father_id,
    child_id,
    trio_file="results/segmental_inference/valid_trios.txt",
    base_path="/home/abiddan1/scratch16/natera_aneuploidy/analysis/aneuploidy/results/natera_inference",
):
    trio_df = pd.read_csv(trio_file, sep="\t")
    filt_df = trio_df[
        (trio_df.mother == mother_id)
        & (trio_df.father == father_id)
        & (trio_df.child == child_id)
    ]
    assert filt_df.shape[0] > 0
    return f"{base_path}/{mother_id}+{father_id}/{child_id}.hmm_model.pkl.gz"


rule est_segmental_changepoints:
    """Estimate segmental aneuploidy using changepoint detection on the HMM trace."""
    input:
        triplets="results/segmental_inference/valid_trios.txt",
        hmm_trace=lambda wildcards: define_triplets(
            mother_id=wildcards.mother,
            father_id=wildcards.father,
            child_id=wildcards.child,
        ),
    output:
        changept_calls="results/segmental_inference/{mother}+{father}+{child}.segmental.tsv",
    resources:
        time="1:00:00",
        mem_mb="5G",
    params:
        chroms=chroms,
        min_size=20,
        penalty=10,
    script:
        "scripts/changept_segmental.py"


rule collapse_segmental_calls:
    """Collapse all of the segmental calls into their postQC components."""
    input:
        total_data,
    output:
        "results/raw_segmental_calls/segmental_calls.tsv.gz",
    shell:
        "find results/segmental_inference/ -name \"*.tsv\"  | while read line; do cat $line; done | awk '!visited[$0]++' | gzip > {output}"


rule gzip_table:
    """Gzipping of tables - since polars does not by default."""
    input:
        segmental_calls="results/tables/{tag}.tsv",
    output:
        segmental_calls="results/tables/{tag}.tsv.gz",
    shell:
        "gzip -c {input} > {output}"


rule merge_aneuploidy_w_segmental:
    """Merging full aneuploidy calls and segmental findings."""
    input:
        segmental_calls="results/raw_segmental_calls/segmental_calls.tsv.gz",
        aneuploidy_calls=config["aneuploidy_calls"],
    output:
        raw_segmental_calls=temp("results/tables/segmental_calls_raw.tsv"),
        postqc_segmental_calls=temp("results/tables/segmental_calls_filt.tsv"),
    params:
        ppThresh=0.80,
        lengthThresh=5e6,
        nsnps=100,
    run:
        seg_df = pl.read_csv(input.segmental_calls, null_values=["NA"], separator="\t")
        aneu_df = pl.read_csv(
            input.aneuploidy_calls, null_values=["NA"], separator="\t"
        )
        seg_df = seg_df.with_columns(
            (pl.col("end") - pl.col("start")).alias("segment_size")
        )
        merged_seg_df = seg_df.join(aneu_df, on=["mother", "father", "child", "chrom"])
        merged_seg_df.write_csv(
            output.raw_segmental_calls, separator="\t", null_value="NA"
        )
        filt_seg_df = merged_seg_df.filter(
            (pl.col("mean_posterior") >= params.ppThresh)
            & (pl.col("segment_size") >= params.lengthThresh)
            & (pl.col("nsnps") > params.nsnps)
            & (pl.col("karyotype") != pl.col("bf_max_cat"))
        )
        filt_seg_df.write_csv(
            output.postqc_segmental_calls, separator="\t", null_value="NA"
        )


rule split_filt_segmental_data:
    """Split the filtered segmental data for analyses."""
    input:
        segmental_calls=rules.merge_aneuploidy_w_segmental.output.postqc_segmental_calls,
    output:
        segmental_calls_split=temp(
            expand(
                "results/tables_split/segmental_calls_filt.{x}.tsv",
                x=range(nchunks + 1),
            )
        ),
    run:
        df = pl.read_csv(input.segmental_calls, null_values=["NA"], separator="\t")
        n = df.shape[0]
        chunk_size = int(n / nchunks)
        for ix, frame in enumerate(df.iter_slices(n_rows=chunk_size)):
            frame.write_csv(
                f"results/tables_split/segmental_calls_filt.{ix}.tsv", separator="\t"
            )


rule segment_endpoint_refinement:
    """Use the karyohmm path-traces to map the quantiles of the segment endpoints."""
    input:
        segmental_calls="results/tables_split/segmental_calls_filt.{x}.tsv",
    output:
        segmental_endpoints="results/tables_split/segmental_postqc_probendpoints.{x}.tsv",
    resources:
        time="0:30:00",
    params:
        base_path="/scratch16/rmccoy22/abiddan1/natera_aneuploidy/analysis/aneuploidy/results/natera_inference",
    script:
        "scripts/endpt_estimation.py"


rule bph_sph_segmental:
    """Estimate the BPH and SPH statuses of each trisomy."""
    input:
        segmental_calls="results/tables_split/segmental_calls_filt.{x}.tsv",
    output:
        segmental_bph_sph="results/tables_split/segmental_postqc_bph_sph.{x}.tsv",
    resources:
        time="0:30:00",
    params:
        base_path="/scratch16/rmccoy22/abiddan1/natera_aneuploidy/analysis/aneuploidy/results/natera_inference",
    script:
        "scripts/bph_sph_segmental.py"


rule collect_endpoint_refinement:
    input:
        expand(
            "results/tables_split/segmental_postqc_probendpoints.{x}.tsv",
            x=range(nchunks + 1),
        ),
    output:
        "results/tables/segmental_postqc_probendpoints.tsv",
    shell:
        "awk '!a[$0]++' {input} > {output}"


rule collect_bph_sph_segmental:
    input:
        expand(
            "results/tables_split/segmental_postqc_bph_sph.{x}.tsv",
            x=range(nchunks + 1),
        ),
    output:
        "results/tables/segmental_postqc_bph_sph.tsv",
    shell:
        "awk '!a[$0]++' {input} > {output}"


rule merge_qc_tables:
    """Merge all of the post-QC analyses to a much larger table."""
    input:
        segmental_calls=rules.merge_aneuploidy_w_segmental.output.postqc_segmental_calls,
        segmental_bph_sph="results/tables/segmental_postqc_bph_sph.tsv",
        segmental_endpoints="results/tables/segmental_postqc_probendpoints.tsv",
    output:
        post_qc_segmental_calls="results/tables/segmental_calls_postqc.tsv",
    run:
        seg_df = pl.read_csv(
            input.segmental_calls,
            separator="\t",
            infer_schema_length=10000,
            null_values=["NA"],
        )
        bph_sph_df = pl.read_csv(
            input.segmental_bph_sph,
            separator="\t",
            infer_schema_length=10000,
            ignore_errors=True,
            null_values=["NA"],
        )
        endpoints_df = pl.read_csv(
            input.segmental_endpoints,
            separator="\t",
            infer_schema_length=10000,
            ignore_errors=True,
            null_values=["NA"],
        )
        merge1_df = seg_df.join(
            endpoints_df,
            how="left",
            on=["mother", "father", "child", "chrom", "start", "end", "karyotype"],
        )
        merge2_df = merge1_df.join(
            bph_sph_df,
            how="left",
            on=["mother", "father", "child", "chrom", "start", "end"],
        )
        merge2_df.write_csv(
            output.post_qc_segmental_calls,
            separator="\t",
            null_value="NA",
        )


rule final_postqc_filter:
    """Applying one final filtering step for segmental calls."""
    input:
        post_qc_segmental_calls="results/tables/segmental_calls_postqc.tsv",
    output:
        post_qc_segmental_calls_refined="results/tables/segmental_calls_postqc_refined.tsv",
    params:
        ppThresh=0.90,
        lengthThresh=5e6,
        nsnps=100,
    run:
        seg_df = pl.read_csv(
            input.post_qc_segmental_calls,
            separator="\t",
            infer_schema_length=10000,
            null_values=["NA"],
        )
        filt_seg_df = seg_df.with_columns(
            (pl.col("end_50") - pl.col("start_50")).alias("segment_size")
        ).filter(
            (pl.col("segment_size") >= params.lengthThresh)
            & (pl.col("nsnps") > params.nsnps)
            & (pl.col("mean_posterior") > params.ppThresh)
            & (pl.col("karyotype") != pl.col("bf_max_cat"))
        )
        filt_seg_df = filt_seg_df.with_columns(
            pl.concat_str(
                pl.col("mother"), pl.col("father"), pl.col("child"), separator="+"
            ).alias("uid")
        )
        filt_seg_df.write_csv(
            output.post_qc_segmental_calls_refined, separator="\t", null_value="NA"
        )
