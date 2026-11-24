import pandas as pd, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
M=pd.read_csv("cgga_out/meta_TCGA_CGGA.csv")
groups=["All gliomas","GBM, IDH-wt","Astrocytoma, IDH-mut","IDH-wt lower-grade","Oligodendroglioma"]
panels=[("HypoxicMES_N6","age (+grade, type)","Hypoxic-MES niche (N6)"),("HypoxicMES_N6","+MP_Hypoxia","N6, beyond hypoxia program"),
        ("Neuronal_N1","age (+grade, type)","Neuronal niche (N1)"),("ANGPTL4","age (+grade, type)","ANGPTL4 (single gene)")]
fig,ax=plt.subplots(1,4,figsize=(17,4.2),sharey=True)
for a,(v,adj,t) in zip(ax,panels):
    x=M[(M.variable==v)&(M.extra_adjustment==adj)].set_index("group").reindex(groups)
    y=np.arange(len(groups))[::-1]
    col=np.where(x.p_meta<0.05,"#B5172F","#5B6C8F")
    a.errorbar(x.HR_meta,y,xerr=[x.HR_meta-x.CI_low,x.CI_high-x.HR_meta],fmt="none",ecolor=col,elinewidth=1.6)
    a.scatter(x.HR_meta,y,c=col,s=40,zorder=3); a.axvline(1,c="grey",ls="--",lw=.8); a.set_xscale("log")
    for yi,(_,r) in zip(y,x.iterrows()): a.text(r.CI_high*1.04 if not np.isnan(r.CI_high) else 1, yi, f"{r.HR_meta:.2f}, p={r.p_meta:.1g}", va="center", fontsize=8)
    a.set_xticks([0.5,1,2,4]); a.set_xticklabels(["0.5","1","2","4"]); a.minorticks_off(); a.set_title(t,fontsize=10); a.set_xlabel("HR per SD (meta, 95% CI)"); a.set_xlim(0.4,5)
ax[0].set_yticks(np.arange(len(groups))[::-1]); ax[0].set_yticklabels(groups)
plt.suptitle("Overall survival, fixed-effect meta-analysis of TCGA + CGGA-693 + CGGA-325 (oligodendroglioma: TCGA + CGGA-693)",fontsize=11)
plt.tight_layout(); plt.savefig("cgga_out/meta_forest.png",dpi=200)
