# -*- coding: utf-8 -*-
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import rcParams
import scipy,random, matplotlib, itertools, glob, os, platform, getpass
import numpy as np
import matplotlib.gridspec as gridspec
from pathlib import Path
from statsmodels.stats import weightstats

if getpass.getuser() == "mengxing":
    git_dir = Path("/home/mengxing/GIT/THATRACT_paper")
elif getpass.getuser() == "lmengxing":
    if platform.system() == "Linux":
        git_dir = Path("/bcbl/home/home_g-m/lmengxing/TESTDATA/GIT/THATRACT_paper")
    elif platform.system() == "Windows":
        git_dir = Path("F:\TESTDATA\GIT\THATRACT_paper")
raw_csv = Path(f"{git_dir}/raw_csv")
fig_dir = Path(f"{git_dir}/fig_dir")

correlation = pd.read_csv( raw_csv / "correlation_fa_TRT_compute.csv")
pairwise = pd.read_csv( raw_csv / "pairwise_all.csv")

# between computational 
pairs = {"L_OR_05":"R_OR_05", "L_AR_A1-4":"R_AR_A1-4", 
            "L_MR_M1":"R_MR_M1", "L_MR_S1":"R_MR_S1",
            "L_DT-4":"R_DT-4"}
effect = ["corr", "density_correlation"]
n = 0
for TCK, btw in itertools.product(pairs, ["T01", "com"]):
    a = correlation[(correlation.TCK==TCK) & 
            (correlation.btw.str.contains(btw))]
    b = correlation[(correlation.TCK==pairs[TCK]) & 
            (correlation.btw.str.contains(btw))]
    x = pd.merge(left=a, right=b, 
                left_on=["SUBID", "btw"], right_on=["SUBID", "btw"])
    a_z = np.arctanh(x["corr_x"])
    b_z = np.arctanh(x["corr_y"])
    t_val, p_val = scipy.stats.ttest_rel(a_z, b_z)
    if p_val<0.05/20: print(TCK, btw, t_val, p_val, x["corr_x"].mean(), x["corr_y"].mean())
    n += 1
n = 0    
for TCK, btw in itertools.product(pairs, ["T01", "com"]):
    a = pairwise[(pairwise.TCK==TCK) & 
            (pairwise.btw.str.contains(btw))]
    b = pairwise[(pairwise.TCK==pairs[TCK]) & 
            (pairwise.btw.str.contains(btw))]
    x = pd.merge(left=a, right=b, 
                left_on=["SUBID", "btw"], right_on=["SUBID", "btw"])
    a_z = np.arctanh(x["density_correlation_x"])
    b_z = np.arctanh(x["density_correlation_y"])
    t_val, p_val = scipy.stats.ttest_rel(a_z, b_z)
    if p_val<0.05/20: print(TCK, btw, np.mean(a_z), np.mean(b_z))
    n += 1
    
n=0


for TCK, btw, ind in itertools.product(pairs, 
                            ["T01", "com"], ["dice_voxels", "bundle_adjacency_voxels"]):
    a = pairwise[(pairwise.TCK==TCK) & 
            (pairwise.btw.str.contains(btw))]
    b = pairwise[(pairwise.TCK==pairs[TCK]) & 
            (pairwise.btw.str.contains(btw))]
    x = pd.merge(left=a, right=b, 
                left_on=["SUBID", "btw"], right_on=["SUBID", "btw"])
    
    
    t_val, p_val = scipy.stats.ttest_rel(x[f"{ind}_x"], x[f"{ind}_y"])
    if p_val<0.05/20: print(TCK, btw,ind, x[f"{ind}_x"].mean(), x[f"{ind}_y"].mean())
    n += 1

    
