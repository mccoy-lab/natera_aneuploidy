import numpy as np
import pandas as pd
from cyvcf2 import VCF
from tqdm import tqdm
import gzip as gz
import pickle
import sys


def obtain_parental_genotypes(vcf_file, mother_id, father_id, af=0.001, threads=4):
    """Obtain the parental genotypes and check that these are in the dataset.

    Returns:
     - rsids
     - pos
     - refs
     - alts
     - mat_haps
     - pat_haps
    """
    rsids = []
    pos = []
    ref = []
    alt = []
    mat_haps = []
    pat_haps = []
    for variant in tqdm(
        VCF(vcf_file, gts012=True, samples=[mother_id], threads=threads)
    ):
        if (variant.aaf > af) | (1 - variant.aaf > af):
            rsids.append(variant.ID)
            pos.append(variant.POS)
            ref.append(variant.REF)
            alt.append(variant.ALT[0])
            mat_haps.append(variant.genotypes[0][:2])
    for variant in tqdm(
        VCF(vcf_file, gts012=True, samples=[father_id], threads=threads)
    ):
        if (variant.aaf > af) | (1 - variant.aaf > af):
            pat_haps.append(variant.genotypes[0][:2])

    # Convert to numpy objects for easier downstream processing
    rsids = np.array(rsids, dtype=str)
    pos = np.array(pos, dtype=np.uint64)
    ref = np.array(ref, dtype=str)
    alt = np.array(alt, dtype=str)
    mat_haps = np.array(mat_haps).T
    pat_haps = np.array(pat_haps).T
    return rsids, pos, ref, alt, mat_haps, pat_haps


def obtain_child_data(child_csv_file, cytosnp_map, allele_file):
    """Read in the child CSV file which contains the x,y,b. fields.

    This assumes that x,y,b are fields in a CSV or .csv.gz file.
    NOTE: we only end up using the B-allele frequency as the major output here.
    NOTE: this does not use the raw xy, b files supplied to us from Natera
    """
    cytosnp_map_df = pd.read_csv(cytosnp_map, sep="\t")
    allele_cytosnp_df = pd.read_csv(
        allele_file, sep="\t|\s", engine="python", header=None
    )
    allele_cytosnp_df.columns = ["rsid", "X1", "X2", "A", "B"]
    child_df = pd.read_csv(child_csv_file)
    child_anno_df = child_df.join(cytosnp_map_df)[
        ["Name", "ChrPosition", "rsid", "b", "x", "y"]
    ].merge(allele_cytosnp_df, how="inner")
    return child_anno_df


def valid_allele(allele):
    """Validate that the allele is not some other character."""
    return allele in ["A", "C", "G", "T"]


def complement(allele):
    """Take the complement of the allele."""
    if allele == "A":
        return "T"
    elif allele == "T":
        return "A"
    elif allele == "C":
        return "G"
    elif allele == "G":
        return "C"
    else:
        raise ValueError(f"Not a correct allele ")


def filter_parent_child_data(
    child_df, mat_haps, pat_haps, rsids, pos, ref, alt, chrY=False
):
    """Filter the resultant parent-child data for this chromosome."""
    assert rsids.size == ref.size
    assert ref.size == alt.size
    assert pos.size == alt.size
    assert (mat_haps.ndim == 2) and (pat_haps.ndim == 2)
    assert mat_haps.shape[1] == rsids.size
    assert pat_haps.shape[1] == rsids.size
    bafs = np.zeros(len(rsids))
    print(child_df[child_df.rsid.isin(rsids)])
    rsid_dict = {}
    for r, baf, B in tqdm(
        zip(child_df.rsid.values, child_df.b.values, child_df.B.values)
    ):
        rsid_dict[r] = (baf, B)
    for i, (r, rx, ax) in tqdm(enumerate(zip(rsids, ref, alt))):
        (cur_baf, b_allele) = rsid_dict[r]
        if chrY:
            cond = valid_allele(b_allele) and (np.sum(pat_haps[:, i]) in [0, 1, 2])
        else:
            cond = (
                valid_allele(b_allele)
                and (np.sum(mat_haps[:, i]) in [0, 1, 2])
                and (np.sum(pat_haps[:, i]) in [0, 1, 2])
            )
        if cond:
            if (b_allele == ax) | (b_allele == complement(ax)):
                bafs[i] = cur_baf
            elif (b_allele == rx) | (b_allele == complement(rx)):
                bafs[i] = 1.0 - cur_baf
            else:
                bafs[i] = np.nan
        else:
            bafs[i] = np.nan
    idx = ~np.isnan(bafs)
    bafs = bafs[idx]
    mat_haps = mat_haps[:, idx]
    pat_haps = pat_haps[:, idx]
    pos = pos[idx]
    ref = ref[idx]
    alt = alt[idx]
    rsids = rsids[idx]
    return bafs, mat_haps, pat_haps, rsids, pos, ref, alt


def main(
    child_csv,
    cytosnp_map,
    alleles_file,
    cytosnp_cluster,
    vcf_file,
    mother_id,
    father_id,
    chrY=False,
):
    """
    Main function for BAF processing.

    Arguments:
        - child_csv: Embryo CSV file.
        - cytosnp_map: CytoSNP v12 Mapping file.
        - alleles_file: Alleles for specific cytosnp probes.
        - cytosnp_cluster: Cytosnp clusters from the EGT file.
        - vcf_file: VCF File containing parental genotypes.
        - mother_id: Mother ID.
        - father_id: Father ID
        - outfile: Output File containing parental haplotypes and embryo BAF
    """
    # Read in the child embryo data
    child_df = obtain_child_data(
        child_csv_file=child_csv,
        cytosnp_map=cytosnp_map,
        allele_file=alleles_file,
    )

    # Obtain the parental genotypes
    print("Obtaining parental genotypes...", file=sys.stderr)
    rsids, pos, ref, alt, mat_haps, pat_haps = obtain_parental_genotypes(
        vcf_file,
        mother_id=mother_id,
        father_id=father_id,
    )
    print("Finished obtaining parental genotypes!", file=sys.stderr)
    print("Filtering parent & embryo data ... ", file=sys.stderr)
    # Filter the parent child data
    baf, mat_haps, pat_haps, rsids, pos, refs, alts = filter_parent_child_data(
        child_df, mat_haps, pat_haps, rsids, pos, ref, alt, chrY=chrY
    )
    print("Finished parent & embryo filtering!", file=sys.stderr)
    # Save the data to an npz file or table
    res_dict = {
        "baf_embryo": baf,
        "mat_haps": mat_haps,
        "pat_haps": pat_haps,
        "rsids": rsids,
        "pos": pos,
        "refs": refs,
        "alts": alts,
        "aploid": "real_data",
    }
    return res_dict


if __name__ == "__main__":
    meta_dict = {}
    for v, c in zip(snakemake.input["vcf_file"], snakemake.params["chroms"]):
        print(f"Processing {c}...", file=sys.stderr)
        chrom_dict = main(
            child_csv=snakemake.input["child_data"],
            cytosnp_map=snakemake.input["cytosnp_map"],
            alleles_file=snakemake.input["alleles_file"],
            cytosnp_cluster=snakemake.input["egt_cluster"],
            vcf_file=v,
            mother_id=snakemake.wildcards["mother_id"],
            father_id=snakemake.wildcards["father_id"],
            chrY=(c == "chrY"),
        )
        chrom_dict["mother_id"] = snakemake.wildcards["mother_id"]
        chrom_dict["father_id"] = snakemake.wildcards["father_id"]
        chrom_dict["child_id"] = snakemake.wildcards["child_id"]
        meta_dict[c] = chrom_dict
    pickle.dump(meta_dict, gz.open(snakemake.output["baf_pkl"], "wb"))
