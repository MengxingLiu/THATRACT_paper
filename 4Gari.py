# -*- coding: utf-8 -*-
"""
Created on Thu May 27 09:27:20 2021

@author: lmengxing
"""
c =  sns.scatterplot(x='x',y='y',data=d)
#c.set_xlim(0.4,0.6)
#c.set_ylim(0.4,0.6)


tmp = profile[(profile["TCK"]=='L_MT_M1')]

tmp = tmp[tmp["subID"].isin(tmp[tmp["ses"]=="T02"].subID.unique())]
ax1 = sns.lineplot(data=tmp, x="ind", y="fa", hue="ses", ci="sd",
                   style="ses")
ax1.set_ylim(0.2789487381975443, 0.607562335478535)
ax1.figure.savefig("test-retest_fa_lineplot.svg")
a = tmp.groupby(by=["subID", "ses"]).mean()

a = a.reset_index()


x = list(a[a["ses"]=="T01"]["fa"])
y = list(a[a["ses"]=="T02"]["fa"])
dic = {"x":x,"y":y}
df_tmp = pd.DataFrame(dic)
plt.gca().set_aspect('equal',adjustable='box')
c = sns.regplot(x="x",y="y",data=df_tmp)
c.set_xlim(0.4,0.53)
c.set_ylim(0.4, 0.53)
# add diagonal line
x0,x1 = c.axes.get_xlim()
y0,y1 = c.axes.get_ylim()
lims = [max(x0,y0),min(x1,y1)]
c.plot(lims,lims)
c.figure.savefig("test-retest_fa_scatterplot.svg")



# for first and second computation
tmp = profile[(profile["analysis"].isin([1,2])) & (profile["TCK"]=='L_MT_M1') &
              (profile["ses"]=="T01")]
tmp["analysis"] = tmp["analysis"].map(lambda x : f"compute_0{str(x)}")
ax1 = sns.lineplot(data=tmp, x="ind", y="fa", hue="analysis", ci="sd",
                   style="analysis")
ax1.figure.savefig("computational_fa_lineplot.svg")

a = tmp.groupby(by=["subID", "analysis"]).mean()

a = a.reset_index()

x = list(a[a["analysis"]=="compute_01"]["fa"])
y = list(a[a["analysis"]=="compute_01"]["fa"])
dic = {"x":x,"y":y}
df_tmp = pd.DataFrame(dic)
plt.gca().set_aspect('equal',adjustable='box')
c = sns.regplot(x="x",y="y",data=df_tmp)


# add diagonal line
x0,x1 = c.axes.get_xlim()
y0,y1 = c.axes.get_ylim()
lims = [max(x0,y0),min(x1,y1)]
c.plot(lims,lims)
c.set_xlim(0.4,0.53)
c.set_ylim(0.4, 0.53)
c.figure.savefig("computational_fa_scatterplot.svg")
