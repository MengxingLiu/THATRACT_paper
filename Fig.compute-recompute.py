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
import matplotlib.ticker as ticker

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
               "L_MR_M1", "R_MR_M1", "L_DT", "R_DT", "L_SR", "R_SR"]

def plot_fig(btw, profile=profile, correlation=correlation, 
                pairwise=pairwise, tck_to_plot=tck_to_plot):
    
    profile = profile[profile["TCK"].isin(tck_to_plot)]
    correlation = correlation[correlation["TCK"].isin(tck_to_plot)]
    pairwise = pairwise[pairwise["TCK"].isin(tck_to_plot)]
    # add two columns to separate fiber group and hemi
    profile[["hemi","group"]] = profile["TCK"].str.split("_",1, expand=True)
    correlation[["hemi","group"]] = correlation["TCK"].str.split("_",1, expand=True)
    pairwise[["hemi","group"]] = pairwise["TCK"].str.split("_",1, expand=True)
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
    new_palette = palette.copy()
    order = tck_to_plot.copy()
    i = 2
    
    while i < len(order):
        # add space between bundle groups
        order.insert(i, "")
        order.insert(i+1, "")
        new_palette.insert(i, palette[0]);new_palette.insert(i, palette[0])
        i += (3+1)

    # FA profile 
    print(order)
    dodge = 0.2
    f2_ax1 = sns.lineplot(data=tmp, x="ind", y="fa", hue=hue_pro,
                    ci="sd", style = hue_pro, palette = ["Grey", "Green"] )
    f2_ax1.set_title('gs[0, :]')
    f2_ax2 = fig2.add_subplot(gs[3:6, 0])

    # FA correlation dist
    f2_ax2 = sns.stripplot(x="TCK", y = "corr", data = correlation,
                        order = order, palette = new_palette,
                        rasterized=True, size=size)
    f2_ax2.set(xlabel="fiber group label")

    f2_ax3 = fig2.add_subplot(gs[0:2, 1])
    f2_ax3 = sns.stripplot(x="TCK", y = "bundle_adjacency_voxels", data = pairwise,
                        order= order, palette = new_palette,
                        rasterized=True, size=size)
    f2_ax3.set(xticklabels=[])
    f2_ax3.set(xlabel=None)
    if btw=="compute": f2_ax3.set_yticks([0,0.5,1.0,1.5])

    f2_ax4 = fig2.add_subplot(gs[2:4, 1])
    f2_ax4 = sns.stripplot(x="TCK", y = "dice_voxels", data = pairwise,
                        order = order, palette = new_palette, 
                        rasterized=True, size=size)
    f2_ax4.set(xticklabels=[])
    f2_ax4.set(xlabel=None)

    f2_ax5 = fig2.add_subplot(gs[4:6, 1])
    f2_ax5 = sns.stripplot(x="TCK", y = "density_correlation", data = pairwise,
                        order = order,  palette = new_palette, 
                        rasterized=True, size=size)
    f2_ax5.set(xlabel="fiber group label")
    return fig2 

fig2 = plot_fig("compute", profile, correlation, pairwise, tck_to_plot)

fig2.set_size_inches(9.88,5.93)
fig2.savefig( fig_dir / "Fig2_computational_new.svg", dpi=300, bbox_inches='tight')


fig3 = plot_fig("test-retest", profile, correlation, pairwise, tck_to_plot)
fig3.set_size_inches(9.88,5.93)
fig3.savefig( fig_dir / "Fig3_test-retest_new.svg", dpi=300, bbox_inches='tight')

