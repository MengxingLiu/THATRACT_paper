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

""" calculate results for compute vs recompute
tractDic = {"KN27":"L_OR_05", "KN28":"R_OR_05", "KN29":"L_OR_1", "KN30":"R_OR_1",
            "KN31":"L_AR_belt-3", "KN32":"L_AR_belt-4",
            "KN33":"L_AR_A1-3", "KN34": "L_AR_A1-4",
            "KN35":"R_AR_belt-3", "KN36":"R_AR_belt-4",
            "KN37":"R_AR_A1-3", "KN38":"R_AR_A1-4",
            "KN39":"L_DT-3", "KN40":"L_DT-4",
            "KN41":"R_DT-3", "KN42":"R_DT-4",
            "KN43":"L_MR_M1", "KN44":"R_MR_M1",
            "KN45":"L_MR_dlPreM", "KN46":"R_MR_dlPreM", 
            "KN47":"L_MR_S1", "KN48":"R_MR_S1"      }


list1 = [f"{i:02d}" for i in range(1,11)]
list1 = range(1,11)
a = itertools.combinations(list1, 2)
csv_dir = "F:\TESTDATA\GIT\THATRACT_paper\csv"
df = pd.read_csv(f"{csv_dir}\RTP_Profile_compute.csv")

con_fa_ana = pd.DataFrame(columns = ["subID", "TCK", "corr", "btw", "ind"])
inds = ['ad', 'cl', 'md', 'volume', 'curvature',
       'rd', 'fa', 'torsion']
for i in a:

    df_ana_x = df[ (df["ses"] == "T01") & (df["analysis"] == i[0])]
    df_ana_y = df[ (df["ses"] == "T01") & (df["analysis"] == i[1])]
    TCKS = df_ana_y.TCK.unique()
    ses = df_ana_y.ses.unique()
    SUBS = df_ana_y.subID.unique()
    
    # correlation of fa
    for tck, sub in itertools.product(TCKS, SUBS):
        a = df_ana_x.loc[(df_ana_x["TCK"]==tck) & (df_ana_x["ses"]=="T01") 
                    & (df_ana_x["subID"]==sub), inds]
        a = a.reset_index(drop=True)
        if len(a)==0 or a.isnull().values.any():
            continue
        b = df_ana_y.loc[(df_ana_y["TCK"]==tck) & (df_ana_y["ses"]=="T01") 
                    & (df_ana_y["subID"]==sub), inds]
        if len(b)==0 or b.isnull().values.any():
            continue
        b = b.reset_index(drop=True)
        c = a.corrwith(b)
        d = a.corrwith(b.reindex(index=b.index[::-1]).reset_index(drop=True))
        result = c if c.sum()> d.sum() else d
        print(tck, sub, f"{i[0]} and {i[1]}")
        result["subID"]=sub;result["TCK"]=tck;
        result["btw"]=f"{i[0]}vs{i[1]}";
        result = pd.DataFrame(result).transpose()
        result = pd.melt(result, id_vars=["subID","TCK","btw"], 
                         value_vars = inds, value_name = "corr", var_name="ind" )
        
        con_fa_ana = con_fa_ana.append(result, ignore_index=True)
        
con_fa_ana["TCK"] = con_fa_ana.TCK.map(lambda x : "L"+x)
con_fa_ana = con_fa_ana.replace({"TCK":tractDic} )
con_fa_ana = con_fa_ana.replace({"TCK":newlabel})
con_fa_ana.to_csv(f"{csv_dir}\correlation_all_indices_compute.csv", index=False)

" calculate descriptions"

correlation = pd.read_csv(f"{csv_dir}\correlation_all_indices_compute.csv")

tck_to_plot = ["L_OR_05", "R_OR_05", "L_AR_A1-4", "R_AR_A1-4", 
               "L_MR_M1", "R_MR_M1", "L_DT-4", "R_DT-4"]

correlation = correlation[correlation["TCK"].isin(tck_to_plot)]
correlation["TCK"] = correlation["TCK"].map(lambda x : x[0:4])
                                                 

correlation.groupby("TCK").apply(np.mean)
correlation.groupby("TCK").describe().to_csv("correlation_compute_description.csv")
pairwise.groupby("tract").describe().to_csv("pairwise_agreement_compute_description.csv")
         
"""

""" calculate results for test vs retest
csv_dir = "F:\TESTDATA\THATRACT_paper\csv"

df = pd.read_csv(f"{csv_dir}\RTP_Profile.csv")

con_fa_ana = pd.DataFrame(columns = ["subID", "TCK", "corr", "btw"])

df_ana_x = df[ (df["ses"] == "T01") & (df["analysis"] == 1)]
df_ana_y = df[ (df["ses"] == "T02") & (df["analysis"] == 1)]
TCKS = df_ana_y.TCK.unique()
ses = df_ana_y.ses.unique()
SUBS = df_ana_y.subID.unique()
inds = ['ad', 'cl', 'md', 'volume', 'curvature',
       'rd', 'fa', 'torsion']
# correlation of fa
for tck, sub in itertools.product(TCKS, SUBS):
    a = df_ana_x.loc[(df_ana_x["TCK"]==tck) & (df_ana_x["ses"]=="T01") 
                & (df_ana_x["subID"]==sub), inds]
    if len(a)==0 or a.isnull().values.any():
        continue
    b = df_ana_y.loc[(df_ana_y["TCK"]==tck) & (df_ana_y["ses"]=="T02") 
                & (df_ana_y["subID"]==sub), inds]
    if len(b)==0 or b.isnull().values.any():
        continue
    a = a.reset_index(drop=True)
    b = b.reset_index(drop=True)

    c = a.corrwith(b)
    d = a.corrwith(b.reindex(index=b.index[::-1]).reset_index(drop=True))
    result = c if c.sum()> d.sum() else d
    print(tck, sub)
    result["subID"]=sub;result["TCK"]=tck;
    result["btw"]=f"T01vsT02";
    result = pd.DataFrame(result).transpose()
    result = pd.melt(result, id_vars=["subID","TCK","btw"], 
                     value_vars = inds, value_name = "corr", var_name="ind" )
    con_fa_ana = con_fa_ana.append(result, ignore_index=True)
    

con_fa_ana = con_fa_ana.replace({"TCK":tractDic} )
con_fa_ana = con_fa_ana.replace({"TCK":newlabel})
con_fa_ana.to_csv(f"{csv_dir}\correlation_all_indices_test-retest.csv", index=False)  

" calculate descriptions "

correlation = pd.read_csv("correlation_fa.csv")
pairwise = pd.read_csv("pairwise_agreement.csv")
pairwise = pairwise.replace({"tract":tractDic} )

tck_to_plot = ["L_OR_05", "R_OR_05", "L_AR_A1-4", "R_AR_A1-4", 
               "L_MT_M1", "R_MT_M1", "L_DT-4", "R_DT-4"]
profile = profile[profile["TCK"].isin(tck_to_plot)]
correlation = correlation[correlation["TCK"].isin(tck_to_plot)]
correlation["TCK"] = correlation["TCK"].map(lambda x : x[0:4])
                                                 
pairwise = pairwise[pairwise["tract"].isin(tck_to_plot)]
pairwise["tract"] = pairwise["tract"].map(lambda x : x[0:4])
correlation.groupby("TCK").apply(np.mean)

correlation.groupby("TCK").describe().to_csv("correlation_description.csv")
pairwise.groupby("tract").describe().to_csv("pairwise_agreement_description.csv")
           
"""

tck_to_plot = ["L_OR_05", "R_OR_05", "L_AR_A1-4", "R_AR_A1-4", 
               "L_MR_M1", "R_MR_M1", "L_DT", "R_DT", 
                "L_MR_S1", "R_MR_S1"]

 # plot Fig. 2 

correlation = pd.read_csv(raw_csv / "correlation_all_indices_compute.csv")
correlation_DT = pd.read_csv(raw_csv / "correlation_DT_final.csv")
correlation_DT = correlation_DT.melt(id_vars=["SUBID","TCK","btw"],
                        var_name="ind", value_name="corr")
correlatoin = correlation.rename(columns={"subID":"SUBID"}, inplace=True)
correlation = pd.concat( [correlation,correlation_DT])
correlation["ind"] = correlation["ind"].str[:2]
inds = ['ad', 'md', 'rd']
correlation = correlation[correlation["TCK"].isin(tck_to_plot)]
correlation = correlation[correlation["ind"].isin(inds)]
correlation = correlation[~(correlation.btw=="T01vsT02")]
fig = plt.figure(constrained_layout=True)
ax = sns.stripplot(x="TCK", y = "corr", hue = "ind", dodge = True, 
              data = correlation, order = tck_to_plot, rasterized=True)
ax.set(ylabel="r"); ax.set(xlabel=False)
ax.xaxis.grid(True)
fig.set_size_inches(13,4.93)

fig.savefig(fig_dir / "SM_Fig1_correlation_computation.svg", dpi=300,
            format = "svg", bbox_inches='tight')


### plot Fig.3 
correlation = pd.read_csv(raw_csv / "correlation_all_indices_test-retest.csv")
correlatoin = correlation.rename(columns={"subID":"SUBID"}, inplace=True)
correlation = pd.concat( [correlation,correlation_DT])
correlation["ind"] = correlation["ind"].str[:2]
correlation = correlation[correlation["TCK"].isin(tck_to_plot)]
correlation = correlation[correlation["ind"].isin(inds)]
correlation = correlation[correlation.btw=="T01vsT02"]                                                 

fig = plt.figure(constrained_layout=True)
ax = sns.stripplot(x="TCK", y = "corr", hue = "ind", dodge = 0.05, 
              data = correlation, order = tck_to_plot, rasterized=True)
ax.set(ylabel="r"); ax.set(xlabel=False)
ax.xaxis.grid(True)
fig.set_size_inches(13,4.93)

fig.savefig(fig_dir / "SM_Fig2_correlation_test-retest.svg", dpi=300,
            format = "svg", bbox_inches='tight')

