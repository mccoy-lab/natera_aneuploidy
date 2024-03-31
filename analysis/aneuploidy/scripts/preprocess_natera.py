import click
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from cyvcf2 import VCF
from tqdm import tqdm
import gzip as gz
import pickle
import sys


def obtain_parental_genotypes(vcf_file, mother_id, father_id, af=0.01, threads=4):
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


def obtain_R_expected_raw(child_anno_df, mean_R_fp=None):
    """Obtain R-expected on the raw-intensity scale to pre-process.
    NOTE: this uses the pre-created files from Daniel now in the `R_expected` directory
    """
    if mean_R_fp is None:
        child_anno_df["LRR_raw"] = np.nan
        return child_anno_df
    else:
        meanR_df = pd.read_csv(mean_R_fp, sep="\t")
        meanR_df.columns = ["Name", "ChrPosition", "R_baf13", "R_baf_mid", "R_baf_23"]
        child_anno_df = child_anno_df.merge(meanR_df, how="inner")
        r_raw_expected = np.empty(child_anno_df.shape[0])
        r_raw_expected[:] = np.nan
        i = 0
        for b, r13, rmid, r23 in tqdm(
            zip(
                child_anno_df.b.values,
                child_anno_df.R_baf13.values,
                child_anno_df.R_baf_mid.values,
                child_anno_df.R_baf_23.values,
            )
        ):
            f_r = interp1d([0, 0.3, 0.5, 0.6, 1.0], [r13, r13, rmid, r23, r23])
            r_raw_expected[i] = f_r(b)
            i += 1
        child_anno_df["R_exp_raw"] = r_raw_expected
        child_anno_df["R_raw"] = child_anno_df.x.values + child_anno_df.y.values
        child_anno_df["LRR_raw"] = np.log2(
            child_anno_df["R_raw"].values / child_anno_df["R_exp_raw"].values
        )
        return child_anno_df


def baf2thetas(baf, thetas=np.array([0.0, 0.5, 1.0])):
    """Transform BAF back into the theta-space approximately."""
    assert thetas.size == 3
    if baf == 0.0:
        return thetas[0]
    elif baf == 1.0:
        return thetas[2]
    else:
        theta = (baf * (thetas[1] - thetas[0])) / 0.5 + thetas[0]
        if (theta >= thetas[0]) and (theta <= thetas[1]):
            return theta
        theta = ((baf - 0.5) * (thetas[2] - thetas[1])) / 0.5 + thetas[1]
        if (theta >= thetas[1]) and (theta <= thetas[2]):
            return theta
        else:
            return baf


def calculate_R_embryo(child_anno_df, cluster_df, norm_xy=None, raw_xy=None):
    """Calculate the underlying R-value for an embryo using the normalized space."""
    if (norm_xy is None) or (raw_xy is None):
        return child_anno_df
    else:
        rsdict = {}
        for r, aa_r, ab_r, bb_r in tqdm(
            zip(
                cluster_df.rsid.values,
                cluster_df.AA_R_mean.values,
                cluster_df.AB_R_mean.values,
                cluster_df.BB_R_mean.values,
            )
        ):
            rsdict[r] = np.array([aa_r, ab_r, bb_r])
        norm_x_dict = {}
        norm_y_dict = {}
        raw_x_dict = {}
        raw_y_dict = {}
        # NOTE: these gzip operations could potentially be made using pigz
        with gz.open(norm_xy, "rt") as norm_intensities:
            norm_intensities.readline()
            for line in tqdm(norm_intensities):
                lnsplt = line.split()
                rsid_norm = lnsplt[0]
                A_intensities_norm = np.array(lnsplt[3:][::2], dtype=np.float32)
                B_intensities_norm = np.array(lnsplt[3:][1::2], dtype=np.float32)
                if rsid_norm in rsdict:
                    r_real = A_intensities_norm + B_intensities_norm
                    scaling = np.median(r_real) / np.median(rsdict[rsid_norm])
                    A_intensities_scaled = A_intensities_norm / scaling
                    B_intensities_scaled = B_intensities_norm / scaling
                    norm_x_dict[rsid_norm] = A_intensities_scaled
                    norm_y_dict[rsid_norm] = B_intensities_scaled
        with gz.open(raw_xy, "rt") as raw_intensities:
            raw_intensities.readline()
            for line in tqdm(raw_intensities):
                lnsplt = line.split()
                rsid_raw = lnsplt[0]
                if rsid_raw in norm_x_dict:
                    A_intensities_raw = np.array(lnsplt[3:][::2], dtype=np.float32)
                    B_intensities_raw = np.array(lnsplt[3:][1::2], dtype=np.float32)
                    raw_x_dict[rsid_raw] = A_intensities_raw
                    raw_y_dict[rsid_raw] = B_intensities_raw
        r_embryo = np.empty(child_anno_df.shape[0])
        x_norm = np.empty(child_anno_df.shape[0])
        y_norm = np.empty(child_anno_df.shape[0])
        r_embryo[:] = np.nan
        x_norm[:] = np.nan
        y_norm[:] = np.nan
        i = 0
        for r, b, x, y in tqdm(
            zip(
                child_anno_df.rsid.values,
                child_anno_df.b.values,
                child_anno_df.x.values,
                child_anno_df.y.values,
            )
        ):
            if (r in norm_x_dict) and (r in raw_x_dict):
                fx = interp1d(raw_x_dict[r], norm_x_dict[r], fill_value="extrapolate")
                fy = interp1d(raw_y_dict[r], norm_y_dict[r], fill_value="extrapolate")
                x_n, y_n = fx(x), fy(y)
                x_norm[i] = x_n
                y_norm[i] = y_n
                r_embryo[i] = x_n + y_n
            i += 1
        # Set these terms in the dataframe moving forward
        child_anno_df["x_norm"] = x_norm
        child_anno_df["y_norm"] = y_norm
        child_anno_df["R_norm"] = r_embryo
        return child_anno_df


def interpolate_r_expected(theta, thetas, rs):
    """Interpolate the R-expected values."""
    assert thetas.size > 0
    assert thetas.size == rs.size
    f = interp1d(thetas, rs, fill_value="extrapolate")
    return f(theta)


def obtain_r_expected(child_anno_df, cluster_df, norm_xy=None):
    """Calculate the R-expected value for the embryo data."""
    if norm_xy is None:
        child_anno_df["LRR_norm"] = np.nan
        return child_anno_df
    else:
        # Iterate through the variants
        for x in ["rsid", "b", "x", "y"]:
            assert x in child_anno_df.columns
        for x in [
            "rsid",
            "AA_theta_mean",
            "AA_R_mean",
            "AB_theta_mean",
            "AB_R_mean",
            "BB_theta_mean",
            "BB_R_mean",
        ]:
            assert x in cluster_df.columns
        r_expected = []
        thetas_dict = {}
        rs_dict = {}
        for r, aa_r, ab_r, bb_r, aa_theta, ab_theta, bb_theta in tqdm(
            zip(
                cluster_df.rsid.values,
                cluster_df.AA_R_mean.values,
                cluster_df.AB_R_mean.values,
                cluster_df.BB_R_mean.values,
                cluster_df.AA_theta_mean.values,
                cluster_df.AB_theta_mean.values,
                cluster_df.BB_theta_mean.values,
            )
        ):
            rs_dict[r] = np.array([aa_r, ab_r, bb_r])
            thetas_dict[r] = np.array([aa_theta, ab_theta, bb_theta])
        for r, b in tqdm(zip(child_anno_df.rsid.values, child_anno_df.b.values)):
            if r in thetas_dict:
                # Do the interpolation in the theta-space
                theta = baf2thetas(b, thetas_dict[r])
                cur_r_exp = interpolate_r_expected(theta, thetas_dict[r], rs_dict[r])
                r_expected.append(cur_r_exp)
            else:
                r_expected.append(np.nan)
        r_expected = np.array(r_expected)
        child_anno_df["R_expected_norm"] = r_expected
        if "R_norm" in child_anno_df.columns:
            child_anno_df["LRR_norm"] = np.log2(
                child_anno_df.R_norm.values / child_anno_df.R_expected_norm.values
            )
        else:
            child_anno_df["LRR_norm"] = np.nan
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


def filter_parent_child_data(child_df, mat_haps, pat_haps, rsids, pos, ref, alt):
    """Filter the resultant parent-child data for this chromosome."""
    assert rsids.size == ref.size
    assert ref.size == alt.size
    assert pos.size == alt.size
    assert (mat_haps.ndim == 2) and (pat_haps.ndim == 2)
    assert mat_haps.shape[1] == rsids.size
    assert pat_haps.shape[1] == rsids.size
    bafs = np.zeros(len(rsids))
    lrrs_raw = np.empty(len(rsids))
    lrrs_norm = np.empty(len(rsids))
    rsid_dict = {}
    for r, lrr_raw, lrr_norm, baf, B in tqdm(
        zip(
            child_df.rsid.values,
            child_df.LRR_raw.values,
            child_df.LRR_norm.values,
            child_df.b.values,
            child_df.B.values,
        )
    ):
        rsid_dict[r] = (lrr_raw, lrr_norm, baf, B)
    for i, (r, rx, ax) in tqdm(enumerate(zip(rsids, ref, alt))):
        (lrr_raw, lrr_norm, cur_baf, b_allele) = rsid_dict[r]
        lrrs_raw[i] = lrr_raw
        lrrs_norm[i] = lrr_norm
        if (
            valid_allele(b_allele)
            and (np.sum(mat_haps[:, i]) in [0, 1, 2])
            and (np.sum(pat_haps[:, i]) in [0, 1, 2])
        ):
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
    lrrs_raw = lrrs_raw[idx]
    lrrs_norm = lrrs_norm[idx]
    mat_haps = mat_haps[:, idx]
    pat_haps = pat_haps[:, idx]
    pos = pos[idx]
    ref = ref[idx]
    alt = alt[idx]
    rsids = rsids[idx]
    return bafs, lrrs_raw, lrrs_norm, mat_haps, pat_haps, rsids, pos, ref, alt


# @click.command()
# @click.option('--child_csv', required=True, type=str, help='Embryo CSV file.')
# @click.option('--cytosnp_map', required=True, type=str, help='CytoSNP v12 Mapping file.')
# @click.option('--alleles_file', required=True, type=str, help='Alleles for specific cytosnp probes.')
# @click.option('--cytosnp_cluster', required=True, type=str, help='Cytosnp clusters from the EGT file.')
# @click.option('--norm_xy', required=False, default=None, type=str, help='Normalized XY-intensities.')
# @click.option('--raw_xy', required=False, default=None, type=str, help='Raw XY-intensities.')
# @click.option('--meanr', required=False, default=None, type=str, help='mean-R tables based on raw intensity.')
# @click.option('--vcf_file', required=True, type=str, help='VCF File containing parental genotypes.')
# @click.option('--mother_id', required=True, type=str, help='Mother ID.')
# @click.option('--father_id', required=True, type=str, help='Father ID.')
# @click.option('--outfile', required=True, type=str, help='Output File containing SNP values.')
def main(
    child_csv,
    cytosnp_map,
    alleles_file,
    cytosnp_cluster,
    norm_xy,
    raw_xy,
    meanr,
    vcf_file,
    mother_id,
    father_id,
):
    # Read in the child embryo data
    child_df = obtain_child_data(
        child_csv_file=child_csv,
        cytosnp_map=cytosnp_map,
        allele_file=alleles_file,
    )
    print("Calculating R-values in raw-space ...", file=sys.stderr)
    child_df = obtain_R_expected_raw(child_df, mean_R_fp=meanr)
    print("Finished calculating LRR in raw-space!", file=sys.stderr)
    # Process the LRR for the embryo
    cluster_df = pd.read_csv(cytosnp_cluster, sep="\t")
    print("Calculating R-values in normalized space ...", file=sys.stderr)
    child_df = calculate_R_embryo(child_df, cluster_df, raw_xy=raw_xy, norm_xy=norm_xy)
    print("Finished calculating R-values in normalized space!")
    print("Calculating LRR normalized values ... ", file=sys.stderr)
    child_df = obtain_r_expected(child_df, cluster_df, norm_xy=norm_xy)
    print("Finished calculating LRR normalized values!", file=sys.stderr)

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
    (
        baf,
        lrrs_raw,
        lrrs_norm,
        mat_haps,
        pat_haps,
        rsids,
        pos,
        refs,
        alts,
    ) = filter_parent_child_data(child_df, mat_haps, pat_haps, rsids, pos, ref, alt)
    print("Finished parent & embryo filtering!", file=sys.stderr)
    # Save the data to an npz file or table
    res_dict = {
        "baf_embryo": baf,
        "lrr_embryo_raw": lrrs_raw,
        "lrr_embryo_norm": lrrs_norm,
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
            norm_xy=None,
            raw_xy=None,
            meanr=None,
            vcf_file=v,
            mother_id=snakemake.wildcards["mother_id"],
            father_id=snakemake.wildcards["father_id"],
        )
        meta_dict[c] = chrom_dict
    meta_dict["mother"] = snakemake.wildcards["mother_id"]
    meta_dict["father"] = snakemake.wildcards["father_id"]
    meta_dict["child"] = snakemake.wildcards["child_id"]
    pickle.dump(meta_dict, gz.open(snakemake.output["baf_pkl"], "wb"))
