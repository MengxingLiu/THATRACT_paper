import os
os.environ['FOR_DISABLE_CONSOLE_CTRL_HANDLER'] = '1'
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import rcParams
import random, matplotlib, itertools, glob, platform, getpass
import scipy.stats
import numpy as np
import matplotlib.gridspec as gridspec
from pathlib import Path

if getpass.getuser() == "mengxing":
    git_dir = Path("/home/mengxing/GIT/THATRACT_paper")
elif getpass.getuser() == "lmengxing":
    if platform.system() == "Linux":
        git_dir = Path("/bcbl/home/home_g-m/lmengxing/TESTDATA/GIT/THATRACT_paper")
    elif platform.system() == "Windows":
        git_dir = Path("F:\TESTDATA\GIT\THATRACT_paper")
raw_csv_dir = Path(f"{git_dir}/raw_csv")
# read tractparams to get the target label dictionary
tractparams = pd.read_csv(git_dir / "tractparams_thatract_paper.csv")
tract_dic = dict(zip(tractparams["label"],tractparams["slabel"] ))


## load pairwise_TRT
pairwise = pd.read_csv(raw_csv_dir / "pairwise.csv")
# select only TRT
pairwise = pairwise[pairwise["btw"] == "T01vsT02"]


pairwise = pairwise[~(pairwise["bundle_adjacency_voxels"]=="bundle_adjacency_voxels")]
pairwise["bundle_adjacency_voxels"] = pairwise["bundle_adjacency_voxels"].astype(float)
pairwise["dice_voxels"] = pairwise["dice_voxels"].astype(float)
pairwise["density_correlation"] = pairwise["density_correlation"].astype(float)

## organise noise file
noise = pd.read_csv(raw_csv_dir / "noise.csv")
noise = noise.drop(noise[noise["SUBID"]=="SUBID"].index )
noise["noise"] = noise["noise"].apply(lambda x: x.lstrip("b'").rstrip(r" \n'"))

noise = noise.pivot(index=["SUBID", "TCK"], columns = "ses", values = "noise")
noise = noise.reset_index()
noise = noise.rename(columns = {"T01": "noise_T01", "T02": "noise_T02"})
noise.to_csv(raw_csv_dir / "noise_clean.csv", index=False)

# tckstats organize
tckstats = pd.read_csv(raw_csv_dir / 
                "tckstats_classic.csv")

tckstats["stats"] = tckstats.stats.map(lambda x:x.lstrip("b'").rstrip(r"\\n'"))
tckstats[["tck_mean", "tck_std", "tck_min", "tck_max", "tck_count"]] = (
            tckstats["stats"].str.split(' ',4,  expand=True)
)

tckstats.to_csv(raw_csv_dir / 
            "tckstats_clean.csv", index=False)

## concatenate pairwise agreement of test-retest and only AL06vsAL07 to compaire
## with noise and other stuff

pairwise = pairwise[["bundle_adjacency_voxels", "dice_voxels", 
                'density_correlation', 'TCK', 'SUBID', 'btw']]
# pairwise = pairwise[(pairwise["btw"]=="T01vsT02") | 
#                     (pairwise["btw"]=="comAL_06vscomAL_07")]

pairwise = pairwise.pivot(index=["SUBID", "TCK"], columns = "btw")
pairwise.columns = ["_".join(i) for i in pairwise.columns.to_flat_index()]
pairwise = pairwise.reset_index()
    

# for density correlation, do Fisher z transform
i = "density_correlation"
pairwise[i+"_T01vsT02_Fisher_z"] = np.arctanh(pairwise[i+"_T01vsT02"])

i = "dice_voxels"
pairwise[i+"_T01vsT02_ln"] = np.log10(
        pairwise[i+"_T01vsT02"]/(1-pairwise[i+"_T01vsT02"]))

# concatenate pairwise and noise
noise = pd.read_csv(raw_csv_dir / "noise_clean.csv")
pairwise_noise = pd.merge(pairwise, noise,  how = "outer", on = ["SUBID", "TCK"])

# load head motion, the column 6 is FWD
head_motion = pd.read_csv(raw_csv_dir / "head_motion_all.csv", header=None)
head_motion = head_motion[[6,7,8]]
head_motion.columns = ["fwd", "SUBID", "ses"]
# average head motion across volumes
a = head_motion.groupby(["SUBID", "ses"]).mean().reset_index()
head_motion = a.pivot(index="SUBID", columns = "ses", values = "fwd").reset_index()
head_motion.columns = ["SUBID", "motion_T01", "motion_T02"]

# concatenate pairwise, noise and head_motion
pw_ns_mt = pd.merge(pairwise_noise, head_motion, how = "outer", on = ["SUBID"])


# load fa correlation
fa_corr = pd.read_csv(raw_csv_dir / "correlation_fa.csv")
fa_corr["TCK"]=fa_corr.TCK.str[1:]
fa_corr = fa_corr.replace({"TCK":tract_dic})
fa_corr = fa_corr.rename(columns={"subID":"SUBID"})
#fa_corr = fa_corr.pivot(index=["SUBID", "TCK"], columns= "btw", values = "fa_y")
fa_corr = fa_corr.reset_index()
fa_corr = fa_corr.reindex(sorted(fa_corr.columns), axis=1)

# Fisher z transform
fa_corr["test-retest_AL_07_fa_Fisher_Z"] = np.arctanh(fa_corr["corr"])

# concatenate pairwise, noise, head_motion and fa correlation
pw_fa_ns_mt = pd.merge(pw_ns_mt, fa_corr, how = "outer", on = ["SUBID", "TCK"])

# load tck stats
tckstats = pd.read_csv(raw_csv_dir / "tckstats_clean.csv")
tckstats = tckstats.pivot(index = ["SUBID", "TCK"], columns = "ses")
tckstats.columns = ["_".join(i) for i in tckstats.columns.to_flat_index()]
tckstats = tckstats.reset_index()

# concatenate pairwise, noise, head_motion, fa correlation and tckstats
pw_fa_ns_mt_st = pd.merge(pw_fa_ns_mt, tckstats, 
                        how="outer", on = ["SUBID", "TCK"])
# calculate noise difference between T01 and T02
pw_fa_ns_mt_st["noise_diff_T01-T02"] = pw_fa_ns_mt_st["noise_T01"] - pw_fa_ns_mt_st["noise_T02"]
pw_fa_ns_mt_st.to_csv(raw_csv_dir / "pw_fa_ns_mt_st.csv", index=False)



pairwise = pd.read_csv(raw_csv_dir / "pairwise_agreement_all_final.csv")
# calculate description for pairwise agreement TRT
inds = ["bundle_adjacency_voxels", "dice_voxels", 'density_correlation']
pairwise[pairwise["btw"]=="T01vsT02"].groupby(["TCK"]
            )[inds].describe().to_csv(raw_csv_dir / "pairwise_description_TRT.csv")
# calculate description for pairwise agreement computational
pairwise[~(pairwise["btw"]=="T01vsT02")].groupby(["TCK"]
            )[inds].describe(
            ).to_csv(raw_csv_dir / "pairwise_description_computation.csv")

# calculate description for FA correlation
fa_corr = pd.read_csv(raw_csv_dir / "correlation_all_index.csv")
fa_corr[fa_corr["btw"]=="test-retest_AL_07"].groupby(["TCK"]
            ).describe(
            )["fa_y"].to_csv(raw_csv_dir / "correlation_description_TRT.csv")
fa_corr[~(fa_corr["btw"]=="test-retest_AL_07")].groupby(["TCK"]
            ).describe(
            )["fa_y"].to_csv(raw_csv_dir / "correlation_description_computation.csv")

# organize noise file
noise = pd.read_csv(raw_csv_dir / "noise.csv")
noise = noise.drop(noise[noise["SUBID"]=="SUBID"].index )
noise["noise"] = noise["noise"].apply(lambda x: x.lstrip("b'").rstrip(r" \n'"))
noise = noise.replace({"TCK":tract_dic})
noise.to_csv(git_dir / "noise_clean.csv", index=False)