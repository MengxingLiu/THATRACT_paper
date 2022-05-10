from curses import raw
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

sns.set_theme()


## load pairwise_TRT
pairwise_TRT = pd.read_csv(raw_csv / "pairwise_TRT_DT_2.csv")
pairwise_TRT["btw"] = "T01vsT02"
pairwise_TRT = pairwise_TRT.rename(columns={"tract":"TCK"})
# pairwise_TRT["analysis"] = "01"
pairwise_TRT.groupby(["TCK"]).describe().to_csv(raw_csv / 
                                "pairwise_TRT_DT_description.csv")
def plt_DT(pairwise_TRT, index):
    mean = pairwise_TRT.groupby("TCK").mean()
    mean = mean.sort_values(by = index)
    mean = mean.reset_index()
    order = mean.TCK.unique()
    fig, axes = plt.subplots()
    sns.stripplot(y="TCK", x = index, data = pairwise_TRT, 
                    order=order )
    sns.pointplot(y="TCK", x = index, data = mean, 
                    order = order, join = False)
    plt.xticks(rotation = -90)
    #plt.tight_layout()
    plt.title(index)
    plt.show()

plt_DT(pairwise_TRT, "density_correlation")
plt_DT(pairwise_TRT, "dice_voxels")
plt_DT(pairwise_TRT, "bundle_adjacency_voxels")

pairwise_TRT.groupby("TCK").apply(np.mean).sort_values(by="density_correlation")
pairwise_TRT.groupby("TCK").apply(np.mean).sort_values(by="dice_voxels")
# TRT profile correlation
df = pd.read_csv(raw_csv / "RTP_Profile_DT_final.csv", dtype={"analysis":str})
analysis = df.analysis.unique()
analysis = itertools.combinations(analysis,2)

## calculate RTP profile correlation

between = ["computation", "test-retest"]

inds = ['ad', 'cl', 'md', 'volume', 'curvature',
       'rd', 'fa', 'torsion']
con_fa_ana = pd.DataFrame()

for ana in analysis:
    for BTW in between:
        print(ana, " start")
        if BTW=="computation":
            
            df_ana_x = df[ (df["ses"] == "T01") & (df["analysis"] == ana[0])]
            df_ana_y = df[ (df["ses"] == "T01") & (df["analysis"] == ana[1])]
            btw = f"{ana[0]}vs{ana[1]}"
            df_ana = df_ana_x.merge(df_ana_y, how="inner",
                    on=["subID", "ind", "TCK", "ses"])
        elif ((BTW=="test-retest" ) & (ana == ("31", "32"))):

            df_ana_x = df[ (df["ses"] == "T01") & (df["analysis"] == ana[0])]
            df_ana_y = df[ (df["ses"] == "T02") & (df["analysis"] == ana[0])] 
            df_ana = df_ana_x.merge(df_ana_y, how="inner",
                    on=["subID", "ind", "TCK", "analysis"])
            btw = "T01vsT02"
        else: 
            continue       
        tmp_df = pd.DataFrame()
        for i in inds:
            c = df_ana.groupby(["subID", "TCK"])[i+"_x", i+"_y"].corr().iloc[0::2,-1]
            c = c.reset_index()[["subID", "TCK", i+"_y"]]
            if i == inds[0]:
                tmp_df = c.copy()
            else: tmp_df[i+"_y"] = c[i+"_y"]
        tmp_df["btw"] = btw
        print(ana, " finish")
        con_fa_ana = con_fa_ana.append(tmp_df, ignore_index=True)

con_fa_ana_bk = con_fa_ana.copy()
con_fa_ana_bk.to_csv(raw_csv / "con_fa_ana_DT_raw.csv", index=False)
con_fa_ana = pd.read_csv(raw_csv / "con_fa_ana_DT_raw.csv")
con_filter = con_fa_ana[con_fa_ana["fa_y"]<0.7]
for row in con_filter.itertuples(index=True, name='Pandas'):
    sub = row.subID
    tck = row.TCK
    btw = row.btw
    fa = row.fa_y
    
    if "T01vsT02" in btw:
        df_x = df[((df.TCK==tck) & (df.subID==sub) & (df.ses=="T01") &
                    (df.analysis=="31"))]
        df_y = df[((df.TCK==tck) & (df.subID==sub) & (df.ses=="T02") &
                    (df.analysis=="31"))]
    else:
        df_x = df[((df.TCK==tck) & (df.subID==sub) & (df.ses=="T01") &
                    (df.analysis==btw.split("vs")[0]))]
        df_y = df[((df.TCK==tck) & (df.subID==sub) & (df.ses=="T01") &
                    (df.analysis==btw.split("vs")[1]))]
    b = scipy.stats.pearsonr(df_x["fa"], df_y["fa"][::-1])[0]
    if b < fa:
        
        continue
    else:
        print(sub, tck, btw, fa, b, "changing")
        for i in inds:
            b = scipy.stats.pearsonr(df_x[i], df_y[i][::-1])[0]
            con_fa_ana.loc[((con_fa_ana["TCK"]==tck) & 
                            (con_fa_ana.subID ==sub) &
                            (con_fa_ana.btw ==btw)), i+"_y"] = b
# remove subject 08CAMINO4391 due to failure of reconstruction
con_fa_ana = con_fa_ana[~(con_fa_ana["SUBID"]=="08CAMINO4391")]
plt_DT(con_fa_ana, "fa_y")
con_fa_ana = con_fa_ana.replace({"TCK":{"KN65":"L_DT", "KN67":"R_DT"}} )

con_fa_ana = con_fa_ana.rename(columns={"subID":"SUBID"})
con_fa_ana.to_csv(raw_csv / "correlation_DT_final.csv", index=False)

con_fa_ana.groupby("TCK").describe().to_csv(raw_csv / 
                "correlation_description_TRT_DT_2.csv")

## load pairwise_TRT
pairwise_TRT = pd.read_csv(raw_csv / "pairwise_agreement_SR_TRT.csv")
pairwise_TRT["btw"] = "T01vsT02"
pairwise_TRT = pairwise_TRT.rename(columns={"tract":"TCK"})
# pairwise_TRT["analysis"] = "01"

plt_DT(pairwise_TRT, "density_correlation")
plt_DT(pairwise_TRT, "dice_voxels")
plt_DT(pairwise_TRT, "bundle_adjacency_voxels")




df = pd.read_csv(raw_csv / "RTP_Profile_SR_final.csv", dtype={"analysis":str})
analysis = df.analysis.unique()
analysis = itertools.combinations(analysis,2)

## calculate RTP profile correlation

between = ["computation", "test-retest"]

inds = ['ad', 'cl', 'md', 'volume', 'curvature',
       'rd', 'fa', 'torsion']
con_fa_ana = pd.DataFrame()

for ana in analysis:
    for BTW in between:
        print(ana, " start")
        if BTW=="computation":
            
            df_ana_x = df[ (df["ses"] == "T01") & (df["analysis"] == ana[0])]
            df_ana_y = df[ (df["ses"] == "T01") & (df["analysis"] == ana[1])]
            btw = f"{ana[0]}vs{ana[1]}"
            df_ana = df_ana_x.merge(df_ana_y, how="inner",
                    on=["subID", "ind", "TCK", "ses"])
        elif ((BTW=="test-retest" ) & (ana == ("51", "52"))):

            df_ana_x = df[ (df["ses"] == "T01") & (df["analysis"] == ana[0])]
            df_ana_y = df[ (df["ses"] == "T02") & (df["analysis"] == ana[0])] 
            df_ana = df_ana_x.merge(df_ana_y, how="inner",
                    on=["subID", "ind", "TCK", "analysis"])
            btw = "T01vsT02"
        else: 
            continue       
        tmp_df = pd.DataFrame()
        for i in inds:
            c = df_ana.groupby(["subID", "TCK"])[i+"_x", i+"_y"].corr().iloc[0::2,-1]
            c = c.reset_index()[["subID", "TCK", i+"_y"]]
            if i == inds[0]:
                tmp_df = c.copy()
            else: tmp_df[i+"_y"] = c[i+"_y"]
        tmp_df["btw"] = btw
        print(ana, " finish")
        con_fa_ana = con_fa_ana.append(tmp_df, ignore_index=True)

con_fa_ana_bk = con_fa_ana.copy()
con_fa_ana_bk.to_csv(raw_csv / "con_fa_ana_SR_raw.csv", index=False)
con_fa_ana = pd.read_csv(raw_csv / "con_fa_ana_SR_raw.csv")
con_filter = con_fa_ana[con_fa_ana["fa_y"]<0.7]

plt_DT(con_fa_ana, "fa_y")
con_fa_ana = con_fa_ana.replace({"TCK":{"KN75":"L_SR", "KN76":"R_SR"}} )

con_fa_ana = con_fa_ana.rename(columns={"subID":"SUBID"})
con_fa_ana.to_csv(raw_csv / "correlation_SR_final.csv", index=False)



