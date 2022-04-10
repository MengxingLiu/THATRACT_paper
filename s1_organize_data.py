import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import rcParams
import scipy,random, matplotlib, itertools, glob, os, platform, getpass
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
raw_csv = Path(f"{git_dir}/raw_csv")
# read tractparams to get the target label dictionary
tractparams = pd.read_csv(raw_csv / "tractparams_THATRACT.csv")
tract_dic = dict(zip(tractparams["slabel"], tractparams["roi2"]))
tract_dic = {k:v for k, v in tract_dic.items() if "KN" in k}
## load pairwise_TRT
pairwise_TRT = pd.read_csv(raw_csv / "pairwise_agreement.csv")
pairwise_TRT["btw"] = "T01vsT02"
pairwise_TRT["analysis"] = "01"

pairwise_0607 = pd.read_csv(raw_csv / "pairwise_agreement_compute.csv")
pairwise = pd.concat([pairwise_TRT, pairwise_0607])
#pairwise = pairwise.replace({"tract":tract_dic})
pairwise = pairwise.rename(columns={"tract":"TCK"})
pairwise = pairwise.drop(pairwise[pairwise_0607["bundle_adjacency_voxels"].str.contains("b", na=True)].index)
pairwise["bundle_adjacency_voxels"] = pairwise["bundle_adjacency_voxels"].astype(float)
pairwise["dice_voxels"] = pairwise["dice_voxels"].astype(float)
pairwise["density_correlation"] = pairwise["density_correlation"].astype(float)
pairwise.to_csv(git_dir / "pairwise.csv", index=False)

pairwise = pd.read_csv(git_dir / "pairwise.csv")
pairwise = pairwise[["bundle_adjacency_voxels", "dice_voxels", 
                'density_correlation', 'TCK', 'SUBID', 'btw']]

pairwise = pairwise.pivot(index=["SUBID", "TCK"], columns = "btw")
pairwise.columns = ["_".join(i) for i in pairwise.columns.to_flat_index()]
pairwise = pairwise.reset_index()

## load noise
noise = pd.read_csv(raw_csv / "noise.csv")
noise = noise.drop(noise[noise["SUBID"]=="SUBID"].index )
noise["noise"] = noise["noise"].apply(lambda x: x.lstrip("b'").rstrip(r" \n'"))
#noise = noise.replace({"TCK":tract_dic})
noise = noise.pivot(index=["SUBID", "TCK"], columns = "ses", values = "noise")
noise = noise.reset_index()
noise = noise.rename(columns = {"T01": "noise_T01", "T02": "noise_T02"})
noise["noise_T01"] = noise["noise_T01"].astype(float)
noise["noise_T02"] = noise["noise_T02"].astype(float)
noise.to_csv(git_dir / "noise_clean.csv", index=False)

pairwise_noise = pd.merge(pairwise, noise,  how = "outer", on = ["SUBID", "TCK"])

plot_repro_factors(pairwise_noise,"dice_voxels_T01vsT02", "noise_T01")

def plot_repro_factors(data, repro, factor):
    data = data[((data[repro].notna()) & (data[factor].notna()))]
    a = data.groupby("TCK").mean()
    a = a.sort_values(by=repro).reset_index()
    r = a[repro].corr(a[factor]) 
    r_10 = a[10:-10][repro].corr(a[10:-10][factor])
    b = a[~a["TCK"].str.contains("Mamm")]
    r_no_mam = b[repro].corr(b[factor])
    print(f"{r}, {r_10}, if remove 10 highest and lowest ")
    print(f"{r_no_mam} if remove mammilothalamic tract")
    fig, axes = plt.subplots(figsize=(38,20))
    axes2 = axes.twiny()
    # fa correlation with noise
    sns.stripplot(y="TCK", x=repro, order = a["TCK"], 
                    data=data, alpha = 0.5, ax = axes)
    sns.pointplot(y="TCK", x=repro, order = a["TCK"], 
                    data=a, alpha = 0.55, ax = axes, 
                    join=False, color="Blue")
    sns.pointplot(y="TCK", x=factor, order = a["TCK"], 
                    data=a, alpha = 0.55, ax = axes2, 
                    join=False, color="Green")
    axes.set_title(f"correlation between {repro} and {factor}/ r = {r}")
    plt.xticks(rotation = -90)
    #plt.show() #plt.show()
