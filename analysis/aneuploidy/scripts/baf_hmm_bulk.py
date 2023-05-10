import sys
import numpy as np
from karyohmm import  MetaHMM

if __name__ == "__main__":
    # Read in the input data and params ...
    eps = 10 ** snakemake.params["eps"]
    lrr = snakemake.params["lrr"]
    for i, b in enumerate(snakemake.input["baf"]):
        baf_data = np.load(b)
        print("Running meta HMM parameter estimation ...", file=sys.stderr)
        n01 = np.nansum((baf_data["baf_embryo"] == 1) | (baf_data["baf_embryo"] == 0))
        m = baf_data["baf_embryo"].size
        if (n01 == m):
            print("Warning: all BAF values were either [0,1].", file=sys.stderr)
            res_dict = {
                "0": np.nan,
                "1m": np.nan,
                "1p": np.nan,
                "2m": np.nan,
                "2p": np.nan,
                "2": np.nan,
                "3m": np.nan,
                "3p": np.nan,
                "sigma_baf": np.nan,
                "pi0_baf": np.nan,
                "pi0_lrr": np.nan,
                "lrr_mu": np.nan,
                "lrr_sd": np.nan,
                "aploid": baf_data["aploid"]
            }
        else:
            logr=False
            lrrs = np.ones(baf_data['baf_embryo'].size)
            if lrr == 'raw':
                lrrs = np.nan_to_num(baf_data['lrr_embryo_raw'])
                logr = True
            elif lrr == 'norm':
                lrrs = np.nan_to_num(baf_data['lrr_embryo_norm'])
                logr = True
            hmm = MetaHMM(logr=logr)
            # NOTE: this naively just takes every other snp to reduce runtimes ...
            pi0_est, sigma_est = hmm.est_sigma_pi0(
                bafs=baf_data["baf_embryo"][::2],
                lrrs=lrrs[::2],
                mat_haps=baf_data["mat_haps"][:, ::2],
                pat_haps=baf_data["pat_haps"][:, ::2],
                eps=eps,
                unphased=snakemake.params["unphased"],
                logr=False
            )
            pi0_lrr = np.nan
            lrr_mu = None
            lrr_sd = None
            if logr:
                pi0_lrr, lrr_mu, lrr_sd, _ = hmm.est_lrr_sd(lrrs, niter=50)
            print("Finished meta HMM parameter estimation ...", file=sys.stderr)
            print("Starting meta HMM forward-backward algorithm", file=sys.stderr)
            gammas, states, karyotypes = hmm.forward_backward(
                bafs=baf_data["baf_embryo"],
                lrrs=lrrs,
                mat_haps=baf_data["mat_haps"],
                pat_haps=baf_data["pat_haps"],
                pi0=pi0_est,
                std_dev=sigma_est,
                pi0_lrr=pi0_lrr,
                lrr_mu=lrr_mu,
                lrr_sd=lrr_sd,
                eps=eps,
                unphased=snakemake.params["unphased"],
                logr=logr
            )
            print(
                "Finished running meta HMM forward-backward algorithm ...",
                file=sys.stderr,
            )
            res_dict = hmm.posterior_karyotypes(gammas, karyotypes)
            res_dict["sigma_baf"] = sigma_est
            res_dict["pi0_baf"] = pi0_est
            res_dict["pi0_lrr"] = pi0_lrr
            res_dict["lrr_mu"] = np.nan if lrr_mu is None else ','.join([str(m) for m in lrr_mu])
            res_dict["lrr_sd"] = np.nan if lrr_sd is None else ','.join([str(s) for s in lrr_sd])
            res_dict["aploid"] = baf_data["aploid"]
            res_dict["karyotypes"] = karyotypes
            res_dict["gammas"] = gammas.astype(np.float16)
            res_dict["states"] = np.array(
                [hmm.get_state_str(s) for s in states], dtype=str
            )
        try:
            res_dict["mother_id"] = snakemake.params["mother_id"]
            res_dict["father_id"] = snakemake.params["father_id"]
            res_dict["child_id"] = snakemake.params["child_id"]
        except KeyError:
            pass
        np.savez_compressed(snakemake.output["hmm_out"][i], **res_dict)
