import sys
import numpy as np
import pickle
import gzip as gz
from karyohmm import MetaHMM

if __name__ == "__main__":
    # Read in the input data and params ...
    data = pickle.load(gz.open(snakemake.input["baf_pkl"], "r"))
    full_chrom_hmm_dict = {}
    for c in snakemake.params["chroms"]:
        baf_data = data[c]
        print("Running meta HMM parameter estimation ...", file=sys.stderr)
        n01 = np.nansum((baf_data["baf_embryo"] == 1) | (baf_data["baf_embryo"] == 0))
        m = baf_data["baf_embryo"].size
        if n01 == m:
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
                "aploid": baf_data["aploid"],
            }
        else:
            hmm = MetaHMM()
            pi0_est, sigma_est = hmm.est_sigma_pi0(
                bafs=baf_data["baf_embryo"][::5],
                mat_haps=baf_data["mat_haps"][:, ::5],
                pat_haps=baf_data["pat_haps"][:, ::5],
                pos=baf_data["pos"][::5],
                unphased=snakemake.params["unphased"],
            )
            print("Finished meta HMM parameter estimation ...", file=sys.stderr)
            print("Starting meta HMM forward-backward algorithm", file=sys.stderr)
            gammas, states, karyotypes = hmm.forward_backward(
                bafs=baf_data["baf_embryo"],
                mat_haps=baf_data["mat_haps"],
                pat_haps=baf_data["pat_haps"],
                pos=baf_data["pos"],
                pi0=pi0_est,
                std_dev=sigma_est,
                unphased=snakemake.params["unphased"],
            )
            print(
                "Finished running meta HMM forward-backward algorithm ...",
                file=sys.stderr,
            )
            res_dict = hmm.posterior_karyotypes(gammas, karyotypes)
            res_dict["sigma_baf"] = sigma_est
            res_dict["pi0_baf"] = pi0_est
            res_dict["aploid"] = baf_data["aploid"]
            res_dict["karyotypes"] = karyotypes
            res_dict["gammas"] = gammas.astype(np.float16)
            res_dict["states"] = np.array(
                [hmm.get_state_str(s) for s in states], dtype=str
            )
            res_dict["pos"] = baf_data["pos"]
        try:
            res_dict["mother_id"] = snakemake.params["mother_id"]
            res_dict["father_id"] = snakemake.params["father_id"]
            res_dict["child_id"] = snakemake.params["child_id"]
        except KeyError:
            pass
        full_chrom_hmm_dict[c] = res_dict
    pickle.dump(full_chrom_hmm_dict, gz.open(snakemake.output["hmm_pkl"], "wb"))
