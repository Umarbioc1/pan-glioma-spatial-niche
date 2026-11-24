import pandas as pd, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt, sys
out=sys.argv[1]
D=pd.read_csv(f"{out}/lr_per_section.csv.gz"); M=pd.read_csv(f"{out}/lr_niche_meta.csv")
D=D[D.frac_in>=0.10]
D["grp"]=np.where(D.subtype.astype(str).str.contains("GBM"),"GBM","IDH-mut")
sig=M[M.significant]
# subtype breakdown for significant pairs
T=D.merge(sig[["niche","pair"]]).groupby(["niche","pair","grp"]).agg(n=("auc","size"),median_auc=("auc","median")).unstack("grp")
T.columns=[f"{a}_{b}" for a,b in T.columns]; T=T.reset_index()
T.to_csv(f"{out}/lr_significant_by_subtype.csv",index=False); print(T.to_string())
# dot plot: niches x top pairs
top=sig.sort_values("mean_pct_rank",ascending=False).groupby("niche").head(8)
pairs=list(dict.fromkeys(top.sort_values(["niche","mean_pct_rank"],ascending=[True,False]).pair))
niches=sorted(M.niche.unique(),key=lambda s:int(s.split("_")[1]))
P=M[M.pair.isin(pairs)].set_index(["niche","pair"])
fig,ax=plt.subplots(figsize=(1.0+0.55*len(niches),0.32*len(pairs)+1.5))
for i,n in enumerate(niches):
    for j,p in enumerate(pairs):
        if (n,p) not in P.index: continue
        r=P.loc[(n,p)]
        sc=ax.scatter(i,j,s=max(r["frac_auc_gt0.6"],0.02)*160,c=[r.median_auc],cmap="RdBu_r",vmin=0.35,vmax=0.75,
                      edgecolors="k" if r.significant else "none",linewidths=1.2)
ax.set_xticks(range(len(niches))); ax.set_xticklabels([n.replace("Niche_","N") for n in niches])
ax.set_yticks(range(len(pairs))); ax.set_yticklabels(pairs,fontsize=8); ax.invert_yaxis()
ax.set_xlabel("Niche"); cb=plt.colorbar(sc,ax=ax,shrink=0.4); cb.set_label("median AUC (niche vs rest)")
ax.set_title("Niche-enriched ligand–receptor pairs (CellChatDB)\nsize = fraction of sections with AUC>0.6; black ring = significant",fontsize=9)
plt.tight_layout(); plt.savefig(f"{out}/lr_niche_dotplot.png",dpi=200)
