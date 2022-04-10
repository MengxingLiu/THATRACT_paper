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
# pairwise_TRT["analysis"] = "01"

pairwise_compute = pd.read_csv(raw_csv / "pairwise_agreement_THATRACT.csv")
pairwise_compute = pairwise_compute.drop(pairwise_compute[pairwise_compute["bundle_adjacency_voxels"].str.contains("bundle", na=True)].index)
pairwise = pd.concat([pairwise_TRT, pairwise_compute])

pairwise = pairwise.rename(columns={"tract":"TCK"})
pairwise["bundle_adjacency_voxels"] = pairwise["bundle_adjacency_voxels"].astype(float)
pairwise["dice_voxels"] = pairwise["dice_voxels"].astype(float)
pairwise["density_correlation"] = pairwise["density_correlation"].astype(float)
pairwise = pairwise.replace({"TCK":tractDic})
pairwise.to_csv(git_dir / "pairwise_all.csv", index=False)

pairwise = pd.read_csv(git_dir / "pairwise.csv")
pairwise = pairwise[["bundle_adjacency_voxels", "dice_voxels", 
                'density_correlation', 'TCK', 'SUBID', 'btw']]

""" calculate results for compute vs recompute
tractDic = {"LKN27":"L_OR_05", "LKN28":"R_OR_05", "LKN29":"L_OR_1", "LKN30":"R_OR_1",
            "LKN31":"L_AR_belt-3", "LKN32":"L_AR_belt-4",
            "LKN33":"L_AR_A1-3", "LKN34": "L_AR_A1-4",
            "LKN35":"R_AR_belt-3", "LKN36":"R_AR_belt-4",
            "LKN37":"R_AR_A1-3", "LKN38":"R_AR_A1-4",
            "LKN39":"L_DT-3", "LKN40":"L_DT-4",
            "LKN41":"R_DT-3", "LKN42":"R_DT-4",
            "LKN43":"L_MR_M1", "LKN44":"R_MR_M1",
            "LKN45":"L_MR_dlPreM", "LKN46":"R_MR_dlPreM", 
            "LKN47":"L_MR_S1", "LKN48":"R_MR_S1"      }


list1 = [f"{i:02d}" for i in range(1,11)]
list1 = range(1,11)
a = itertools.combinations(list1, 2)
df = pd.read_csv("RTP_Profile_compute.csv")

con_fa_ana = pd.DataFrame(columns = ["subID", "TCK", "corr", "btw"])

for i in a:
    
    df_ana_x = df[ (df["ses"] == "T01") & (df["analysis"] == i[0])]
    df_ana_y = df[ (df["ses"] == "T01") & (df["analysis"] == i[1])]
    TCKS = df_ana_y.TCK.unique()
    ses = df_ana_y.ses.unique()
    SUBS = df_ana_y.subID.unique()
    ind = 'fa'
    # correlation of fa
    for tck, sub in itertools.product(TCKS, SUBS):
        a = df_ana_x.loc[(df_ana_x["TCK"]==tck) & (df_ana_x["ses"]=="T01") 
                    & (df_ana_x["subID"]==sub), "fa"]
        if len(a)==0 or a.isnull().values.any():
            continue
        b = df_ana_y.loc[(df_ana_y["TCK"]==tck) & (df_ana_y["ses"]=="T01") 
                    & (df_ana_y["subID"]==sub), "fa"]
        if len(b)==0 or b.isnull().values.any():
            continue
        
        c, _ = scipy.stats.pearsonr(a,b)
        d, _ = scipy.stats.pearsonr(a[::-1],b)
        c = max(c,d)    
        print(tck, sub, f"{i[0]} and {i[1]}")
        con_fa_ana = con_fa_ana.append({"subID":sub, "TCK":tck, 
                                              "corr":c, 
                                              "btw":f"{i[0]}vs{i[1]}"}, 
                                             ignore_index=True)
        
con_fa_ana["TCK"] = con_fa_ana.TCK.map(lambda x : "L"+x)
con_fa_ana = con_fa_ana.replace({"TCK":tractDic} )
con_fa_ana = con_fa_ana.replace({"TCK":newlabel})
con_fa_ana.to_csv("correlation_fa_compute.csv", index=False)

" calculate descriptions"

correlation = pd.read_csv("correlation_fa_compute.csv")
pairwise = pd.read_csv("pairwise_agreement_compute.csv")

pairwise = pairwise.replace({"tract":tractDic} )

tck_to_plot = ["L_OR_05", "R_OR_05", "L_AR_A1-4", "R_AR_A1-4", 
               "L_MR_M1", "R_MR_M1", "L_DT-4", "R_DT-4"]
profile = profile[profile["TCK"].isin(tck_to_plot)]
correlation = correlation[correlation["TCK"].isin(tck_to_plot)]
correlation["TCK"] = correlation["TCK"].map(lambda x : x[0:4])
                                                 
pairwise = pairwise[pairwise["tract"].isin(tck_to_plot)]
pairwise["tract"] = pairwise["tract"].map(lambda x : x[0:4])
pairwise["bundle_adjacency_voxels"] = pairwise["bundle_adjacency_voxels"].astype(float)
pairwise["dice_voxels"] = pairwise["dice_voxels"].astype(float)
pairwise["density_correlation"] = pairwise["density_correlation"].astype(float)

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
ind = 'fa'
# correlation of fa
for tck, sub in itertools.product(TCKS, SUBS):
    a = df_ana_x.loc[(df_ana_x["TCK"]==tck) & (df_ana_x["ses"]=="T01") 
                & (df_ana_x["subID"]==sub), "fa"]
    if len(a)==0 or a.isnull().values.any():
        continue
    b = df_ana_y.loc[(df_ana_y["TCK"]==tck) & (df_ana_y["ses"]=="T02") 
                & (df_ana_y["subID"]==sub), "fa"]
    if len(b)==0 or b.isnull().values.any():
        continue
    
    c, _ = scipy.stats.pearsonr(a,b)
    d, _ = scipy.stats.pearsonr(a[::-1],b)
    c = max(c,d)    
    print(tck, sub, "T01 vs T02")
    con_fa_ana = con_fa_ana.append({"subID":sub, "TCK":tck, 
                                          "corr":c, 
                                          "btw":"T01vsT02"}, 
                                         ignore_index=True)
    
con_fa_ana["TCK"] = con_fa_ana.TCK.map(lambda x : "L"+x)
con_fa_ana = con_fa_ana.replace({"TCK":tractDic} )
con_fa_ana = con_fa_ana.replace({"TCK":newlabel})
con_fa_ana.to_csv("correlation_fa.csv", index=False)  

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


correlation_TRT = pd.read_csv( raw_csv / "correlation_fa.csv")
correlation_compute = pd.read_csv(raw_csv / "correlation_fa_compute.csv")
correlation = pd.concat([correlation_TRT, correlation_compute])
correlation = correlation.rename(columns={"subID":"SUBID"})
correlation = correlation.replace({"TCK":{"L_MT_S1":"L_MR_S1", "L_MT_M1":"L_MR_M1",
                            "R_MT_S1":"R_MR_S1","R_MT_M1":"R_MR_M1"}})
correlation.to_csv(raw_csv / "correlation_fa_TRT_compute.csv", index=False)


profile = pd.read_csv( raw_csv / "RTP_Profile_compute.csv")
profile["TCK"] = profile.TCK.map(lambda x : "L"+x if "KN" in x else x)
profile = profile.replace({"TCK":tractDic})
profile = profile.rename(columns={"subID":"SUBID"})
profile.to_csv(raw_csv / "RTP_Profile_compute_newlabel.csv", index=False)

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
