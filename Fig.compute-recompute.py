# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 14:05:49 2021

@author: lmengxing
"""
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
fig_dir = Path(f"{git_dir}/fig_dir")

# plot Fig. 2 
profile = pd.read_csv( raw_csv / "RTP_Profile_compute_newlabel.csv")
correlation = pd.read_csv( raw_csv / "correlation_fa_TRT_compute.csv")
pairwise = pd.read_csv( raw_csv / "pairwise_all.csv")

tck_to_plot = ["L_OR_05", "R_OR_05", "L_AR_A1-4", "R_AR_A1-4", 
               "L_MR_M1", "R_MR_M1", "L_MR_S1", "R_MR_S1","L_DT-4", "R_DT-4"]

def plot_fig(btw, profile=profile, correlation=correlation, 
                pairwise=pairwise, tck_to_plot=tck_to_plot):
    
    profile = profile[profile["TCK"].isin(tck_to_plot)]
    correlation = correlation[correlation["TCK"].isin(tck_to_plot)]
    pairwise = pairwise[pairwise["TCK"].isin(tck_to_plot)]
    size=3
    if btw == "compute":
        tmp = profile[(profile["analysis"].isin([1,2])) & (profile["TCK"]=='L_MR_M1') &
              (profile["ses"]=="T01")]
        tmp["analysis"] = tmp["analysis"].map(lambda x : f"compute_0{str(x)}")
        hue_pro = "analysis"
        correlation = correlation[  ~(correlation["btw"]=="T01vsT02")]
        pairwise = pairwise[  ~(pairwise["btw"]=="T01vsT02")]

    elif btw == "test-retest":
        tmp = profile[(profile["TCK"]=='L_MR_M1')]
        tmp = tmp[tmp["SUBID"].isin(tmp[tmp["ses"]=="T02"].SUBID.unique())]
        tmp = tmp[tmp["analysis"]==1]
        hue_pro = "ses"
        correlation = correlation[  (correlation["btw"]=="T01vsT02")]
        pairwise = pairwise[  (pairwise["btw"]=="T01vsT02")]
    sns.set_style("darkgrid")
    fig2 = plt.figure(constrained_layout=True)
    gs = fig2.add_gridspec(6, 2)
    f2_ax1 = fig2.add_subplot(gs[0:3, 0])
    palette = sns.color_palette("Paired")
    order = tck_to_plot
    # FA profile 


    f2_ax1 = sns.lineplot(data=tmp, x="ind", y="fa", hue=hue_pro,
                    ci="sd", style = hue_pro, palette = ["Grey", "Green"] )
    f2_ax1.set_title('gs[0, :]')
    f2_ax2 = fig2.add_subplot(gs[3:6, 0])

    # FA correlation dist
    f2_ax2 = sns.stripplot(x="TCK", y = "corr", data = correlation,
                        order = order, palette = palette, rasterized=True, size=size)
    f2_ax2.set(xlabel="fiber group label")

    f2_ax3 = fig2.add_subplot(gs[0:2, 1])
    f2_ax3 = sns.stripplot(x="TCK", y = "bundle_adjacency_voxels", data = pairwise,
                        order= order, palette = palette,
                        rasterized=True, size=size)
    f2_ax3.set(xticklabels=[])
    f2_ax3.set(xlabel=None)

    f2_ax4 = fig2.add_subplot(gs[2:4, 1])
    f2_ax4 = sns.stripplot(x="TCK", y = "dice_voxels", data = pairwise,
                        order = order, palette = palette, rasterized=True, size=size)
    f2_ax4.set(xticklabels=[])
    f2_ax4.set(xlabel=None)

    f2_ax5 = fig2.add_subplot(gs[4:6, 1])
    f2_ax5 = sns.stripplot(x="TCK", y = "density_correlation", data = pairwise,
                        order = order, palette = palette, rasterized=True, size=size)
    f2_ax5.set(xlabel="fiber group label")
    return fig2 

fig2 = plot_fig("compute", profile, correlation, pairwise, tck_to_plot)

fig2.set_size_inches(9.88,5.93)
fig2.savefig( fig_dir / "Fig2_computational_new.svg", dpi=300, bbox_inches='tight')


fig3 = plot_fig("test-retest", profile, correlation, pairwise, tck_to_plot)
fig3.set_size_inches(9.88,5.93)
fig3.savefig( fig_dir / "Fig3_test-retest_new.svg", dpi=300, bbox_inches='tight')

"""

 # plot Fig.3 
csv_dir = r"F:\TESTDATA\GIT\THATRACT_paper\csv"
profile = pd.read_csv(f"{csv_dir}\RTP_Profile.csv")
correlation = pd.read_csv(f"{csv_dir}\correlation_fa.csv")
pairwise = pd.read_csv(f"{csv_dir}\pairwise_agreement.csv")

# change tck labels
profile["TCK"] = profile.TCK.map(lambda x : "L"+x)
profile = profile.replace({"TCK":tractDic} )
pairwise = pairwise.replace({"tract":tractDic} )

tck_to_plot = ["L_OR_05", "R_OR_05", "L_AR_A1-4", "R_AR_A1-4", 
               "L_MR_M1", "R_MR_M1", "L_DT-4", "R_DT-4"]
profile = profile[profile["TCK"].isin(tck_to_plot)]


profile = profile[profile["subID"].isin(profile[profile["ses"]=="T02"].subID.unique())]
correlation = correlation.replace({"TCK":{"L_MT_M1":"L_MR_M1","R_MT_M1":"R_MR_M1"}})
correlation = correlation[correlation["TCK"].isin(tck_to_plot)]
correlation["TCK"] = correlation["TCK"].map(lambda x : x[0:4])
                                                 
pairwise = pairwise[pairwise["tract"].isin(tck_to_plot)]
pairwise["tract"] = pairwise["tract"].map(lambda x : x[0:4])

pairwise["bundle_adjacency_voxels"] = pairwise["bundle_adjacency_voxels"].astype(float)
pairwise["dice_voxels"] = pairwise["dice_voxels"].astype(float)
pairwise["density_correlation"] = pairwise["density_correlation"].astype(float)

fig3 = plt.figure(constrained_layout=True)
gs = fig3.add_gridspec(6, 2)
f3_ax1 = fig3.add_subplot(gs[0:3, 0])

# FA profile 
#tmp = profile[(profile["TCK"]=='L_MR_M1') & (profile["subID"]=="S038")]
tmp = profile[(profile["TCK"]=='L_MR_M1')]
tmp = tmp[tmp["subID"].isin(tmp[tmp["ses"]=="T02"].subID.unique())]
f3_ax1 = sns.lineplot(data=tmp, x="ind", y="fa", hue="ses",
                   style="ses", ci="sd", palette = ["Grey", "Green"] )
# FA correlation dist
f3_ax2 = fig3.add_subplot(gs[3:6, 0])
f3_ax2 = sns.stripplot(x="TCK", y = "corr", data = correlation,
                       order = order, palette = palette, rasterized=True)
f3_ax2.set(xlabel="fiber group label")

f3_ax3 = fig3.add_subplot(gs[0:2, 1])

f3_ax3 = sns.stripplot(x="tract", y = "bundle_adjacency_voxels",
                       data = pairwise, order = order, palette = palette,
                       rasterized=True)
f3_ax3.set(xticklabels=[])
f3_ax3.set(xlabel=None)

f3_ax4 = fig3.add_subplot(gs[2:4, 1])
f3_ax4 = sns.stripplot(x="tract", y = "dice_voxels", data = pairwise, 
                       order = order, palette = palette, rasterized=True)
f3_ax4.set(xticklabels=[])
f3_ax4.set(xlabel=None)

f3_ax5 = fig3.add_subplot(gs[4:6, 1])
f3_ax5 = sns.stripplot(x="tract", y = "density_correlation", data = pairwise,
                       order = order, palette = palette, rasterized=True)
f3_ax5.set(xlabel="fiber group label")

fig3.set_size_inches(9.88,4.93)
fig3.savefig(f"{fig_dir}\Fig3_test-retest_new.svg", dpi=300, format = "svg", bbox_inches='tight')
"""