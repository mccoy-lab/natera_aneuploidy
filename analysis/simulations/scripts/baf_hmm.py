import sys

import numpy as np
from karyohmm import MetaHMM

if __name__ == "__main__":
    # Read in the input data and params ...
    baf_data = np.load(snakemake.input["baf"])
    eps = 10 ** snakemake.params["eps"]
    logr = snakemake.params["lrr"]
    hmm = MetaHMM(logr=logr)
    print("Running meta HMM parameter estimation ...", file=sys.stderr)
    n01 = np.nansum((baf_data["baf_embryo"] == 1) | (baf_data["baf_embryo"] == 0))
    m = baf_data["baf_embryo"].size
    if n01 == m:
        print("Warning: all BAF values were either [0,1].", file=sys.stderr)
        res_dict = {
            "0": np.nan,
            "1m": np.nan,
            "1p": np.nan,
            "2": np.nan,
            "2p": np.nan,
            "2m": np.nan,
            "3m": np.nan,
            "3p": np.nan,
            "sigma_est": np.nan,
            "pi0_est": pi0_est,
            "aploid": baf_data["aploid"],
        }
    else:
        if snakemake.params["phase_error"]:
            mat_haps = baf_data["mat_haps_prime"]
            pat_haps = baf_data["pat_haps_prime"]
        else:
            mat_haps = baf_data["mat_haps"]
            pat_haps = baf_data["pat_haps"]
        pi0_est, sigma_est = hmm.est_sigma_pi0(
            bafs=baf_data["baf_embryo"],
            lrrs=baf_data["lrr_embryo"],
            mat_haps=mat_haps,
            pat_haps=pat_haps,
            eps=eps,
            unphased=snakemake.params["unphased"],
            logr=False,
        )
        print("Finished meta HMM parameter estimation!", file=sys.stderr)
        print("Starting meta HMM forward-backward algorithm.", file=sys.stderr)
        if snakemake.params["lrr"]:
            pi0_lrr, lrrs_mu, lrrs_sd, _ = hmm.est_lrr_sd(lrrs=baf_data["lrr_embryo"])
            gammas, states, karyotypes = hmm.forward_backward(
                bafs=baf_data["baf_embryo"],
                lrrs=baf_data["lrr_embryo"],
                lrr_mu=lrrs_mu,
                lrr_sd=lrrs_sd,
                mat_haps=mat_haps,
                pat_haps=pat_haps,
                pi0=pi0_est,
                std_dev=sigma_est,
                pi0_lrr=pi0_lrr,
                eps=eps,
                unphased=snakemake.params["unphased"],
                logr=True,
            )
        else:
            gammas, states, karyotypes = hmm.forward_backward(
                bafs=baf_data["baf_embryo"],
                lrrs=baf_data["lrr_embryo"],
                mat_haps=mat_haps,
                pat_haps=pat_haps,
                pi0=pi0_est,
                std_dev=sigma_est,
                eps=eps,
                unphased=snakemake.params["unphased"],
                logr=False,
            )
        print(
            "Finished running meta HMM forward-backward algorithm ...",
            file=sys.stderr,
        )
        res_dict = hmm.posterior_karyotypes(gammas, karyotypes)
        res_dict["sigma_est"] = sigma_est
        res_dict["pi0_est"] = pi0_est
        res_dict["aploid"] = baf_data["aploid"]
        res_dict["karyotypes"] = karyotypes
        res_dict["gammas"] = gammas.astype(np.float16)
    try:
        res_dict["mother_id"] = snakemake.params["mother_id"]
        res_dict["father_id"] = snakemake.params["father_id"]
        res_dict["child_id"] = snakemake.params["child_id"]
    except KeyError:
        pass
    np.savez_compressed(snakemake.output["hmm_out"], **res_dict)
