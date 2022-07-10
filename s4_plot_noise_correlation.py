from telnetlib import AO
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import rcParams
import scipy,random, matplotlib, itertools, glob, os, platform, getpass, re
import numpy as np
import matplotlib.gridspec as gridspec
from pathlib import Path
plt.rcParams['svg.fonttype'] = 'none'
if getpass.getuser() == "mengxing":
    git_dir = Path("/home/mengxing/GIT/THATRACT_paper")
elif getpass.getuser() == "lmengxing":
    if platform.system() == "Linux":
        git_dir = Path("/bcbl/home/home_g-m/lmengxing/TESTDATA/GIT/THATRACT_paper")
    elif platform.system() == "Windows":
        git_dir = Path("F:\TESTDATA\GIT\THATRACT_paper")
raw_csv = Path(f"{git_dir}/raw_csv")



pw_fa_ns_mt_st = pd.read_csv(raw_csv / "pw_fa_ns_mt_st.csv")
# drop those computational stats
for col in pw_fa_ns_mt_st.columns:
    if re.match(".*AL.*AL", col):
        del pw_fa_ns_mt_st[col]

# calculate group mean of each correlation pair
# first plot reproducibility with noise
data = pw_fa_ns_mt_st.copy()
a = data.groupby("TCK").mean()
a = a.reset_index()

groups = pd.read_csv(raw_csv / "groups.csv")
tcks = []
for x in groups.columns:
    tcks = tcks + list(groups[x][groups[x].notna()])
b = a[a["TCK"].isin(tcks[:-2])]
factors = ["noise_T01", "tck_mean_T01", "tck_count_T01"]
repro = ["test-retest_AL_07_fa_Fisher_Z", "dice_voxels_T01vsT02_ln", 
        "density_correlation_T01vsT02_Fisher_z", "bundle_adjacency_voxels_T01vsT02"]
corr_heat = pd.DataFrame()
for key, value in itertools.product(repro, factors):
    r = a[key].corr(a[value])
    corr_heat = corr_heat.append({"repro":key, value:r}, ignore_index=True)
corr_heat = corr_heat.groupby("repro").max()

TRT = ['test-retest_AL_07_fa_Fisher_Z',
        'bundle_adjacency_voxels_T01vsT02',
        'dice_voxels_T01vsT02_ln',
       'density_correlation_T01vsT02_Fisher_z' ]
COM = ['AL_all_fa_Fisher_Z', 'bundle_adjacency_voxels_comAL_all_average',
       'dice_voxels_comAL_all_average_ln',
       'density_correlation_comAL_all_average_Fisher_z']
fig, ax = plt.subplots(1,2)
cmap = "RdBu_r"
cmap = sns.diverging_palette(0, 230, 90, 60, as_cmap=True)
xlabels = ["noise", "length", "number"]
ylabels = ["FA correlation", "adjacency", "density", "Dice"]
sns.heatmap(corr_heat.loc[TRT], annot=True, fmt=".2f", cmap=cmap, 
                    center=0, vmin=-1, vmax=1,
                    xticklabels = xlabels,yticklabels=ylabels,
                    ax = ax[0])

sns.heatmap(corr_heat.loc[COM], annot=True, fmt=".2f", cmap=cmap, 
                    center=0, vmin=-1, vmax=1,
                    xticklabels = xlabels,yticklabels=ylabels,
                    ax = ax[1])
plt.show()

def corrdot(*args, **kwargs):
    corr_r = args[0].corr(args[1], 'pearson')
    corr_text = f"{corr_r:2.2f}".replace("0.", ".")
    ax = plt.gca()
    #ax.set_axis_off()
    marker_size = abs(corr_r) * 10000
    ax.scatter([.5], [.5], marker_size, [corr_r], alpha=0.6, cmap="RdBu_r",
               vmin=-1, vmax=1, transform=ax.transAxes)
    font_size = abs(corr_r) * 40 + 5
    ax.annotate(corr_text, [.5, .5,],  xycoords="axes fraction",
                ha='center', va='center', fontsize=font_size)
sns.set(style='white', font_scale=1.6)
g = sns.PairGrid(b, aspect=1.4, diag_sharey=False, 
                    x_vars=factors, y_vars=TRT)
g.map_offdiag(corrdot)
#g.map_offdiag(sns.regplot)
plt.show()


plt.close()
sns.set(style='white', font_scale=1.6)
fig, axes = plt.subplots(4,3)
for ind, ax in zip(itertools.product(COM, factors), axes.flatten()):
    corr_r = corr_heat.loc[ind[0]][ind[1]]
    # if it's bundle adjacency, reverse symbol
    if "adjacency" in ind[0]:
        corr_r = -corr_r
    corr_text = f"{corr_r:2.2f}".replace("0.", ".")
    if 0.209 < abs(corr_r) < 0.27:
        sig_text = "*"
    elif 0.27 < abs(corr_r) < 0.34:
        sig_text = "**"
    elif 0.34 < abs(corr_r) < 1: 
        sig_text= "***"
    else: sig_text = ""
    marker_size = abs(corr_r)*10000
    ax.set_axis_off()
    ax.set(xticklabels=[], yticklabels=[])
    ax.set(xlabel=ind[1], ylabel=ind[0])
    ax.scatter([.5],[.5], marker_size, [corr_r],
            alpha = 0.6, vmin=-1,vmax=1, cmap = cmap)
    font_size = abs(corr_r) * 40 + 5       
    ax.annotate(corr_text, [.5, .5,],  xycoords="axes fraction",
                ha='center', va='center', fontsize=font_size)
    ax.annotate(sig_text, [.5, .7,],  xycoords="axes fraction",
                ha='center', va='center',
                color='red', fontsize=20)
    
for ax in axes.flat:
    ax.label_outer()

plt.subplots_adjust(hspace=.0, wspace=0)
plt.show()

# print color bar
a = np.array([[1,-1]])
plt.figure(figsize=(9, 1.5))
img = plt.imshow(a, cmap=cmap)
plt.gca().set_visible(False)
cax = plt.axes([0.1, 0.2, 0.8, 0.6])
plt.colorbar(orientation="horizontal", cax=cax)
plt.show()


classic = ['CFMaj', 'CFMin', 'LAF', 'LCC', 'LCH', 'LCST', 'LIFOF', 'LILF',
            'LORV1V2', 'LORV1V2R3', 'LSLF', 'LTR', 'LUF', 'RAF', 'RCC', 
            'RCH', 'RCST', 'RIFO', 'RILF', 'RORV1', 'RORV1V2',
            'RORV1V2R3', 'RSLF', 'RTR', 'RUF']
corr_subwise=pd.DataFrame()
data = data[(~data.noise_T02.isna()) & (~data.noise_T01.isna()) ]
for key, value in itertools.product(repro, factors):
    r = []
    for i in data.SUBID.unique():
        x = data[(data.SUBID==i) & (data.TCK.isin(classic))][key]
        y = data[(data.SUBID==i) & (data.TCK.isin(classic))][value]
        if len(x)<20:continue
        r_tmp,p = scipy.stats.pearsonr(x,y)
        r.append(r_tmp)
    corr_subwise = corr_subwise.append({"repro":key, value:np.mean(r)}, ignore_index=True)
corr_subwise = corr_subwise.groupby("repro").max()

