import pandas as pd, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
A=pd.read_csv("ncdeg_out/ncdeg_meta.csv"); A=A[(A.focal=="Malignant")&(A.niche=="Niche_6")]
B=pd.read_csv("ncdeg_stateadj/ncdeg_meta.csv"); B=B[B.niche=="Niche_6"]
J=A.merge(B[["gene","mean_beta","p","fdr","sign_consistency"]],on="gene",suffixes=("","_adj"))
J.to_csv("ncdeg_out/niche6_malignant_unadjusted_vs_stateadjusted.csv",index=False)
fig,ax=plt.subplots(1,2,figsize=(13,5.6))
c=np.where(~J.NCDEG,"#C8CDD3",np.where(J.mean_beta>0,"#B5172F","#2A6FB0"))
ax[0].scatter(J.mean_beta,-np.log10(J.fdr),s=9,c=c,lw=0)
for _,r in J[J.NCDEG].sort_values("fdr").head(22).iterrows(): ax[0].annotate(r.gene,(r.mean_beta,-np.log10(r.fdr)),fontsize=7.5,xytext=(3,2),textcoords="offset points")
ax[0].axhline(-np.log10(0.05),ls="--",c="grey",lw=.8); ax[0].axvline(.2,ls=":",c="grey",lw=.8); ax[0].axvline(-.2,ls=":",c="grey",lw=.8)
ax[0].set_xlabel("Niche-6 effect in malignant-rich spots (mean β, log1p CP10k)"); ax[0].set_ylabel("−log10 FDR (section-level test)")
ax[0].set_title("A  Malignant NCDEGs, niche 6 vs other niches\n(adjusted for non-malignant admixture)",loc="left",fontsize=10)
S=J[J.NCDEG]
ax[1].scatter(S.mean_beta,S.mean_beta_adj,c=np.where(S.mean_beta>0,"#B5172F","#2A6FB0"),s=22)
for _,r in S.iterrows():
    if r.gene in ["ANGPTL4","NDRG1","VEGFA","CAV1","CHI3L1","MT1X","BCAN","OLIG1","ID3","ADM","PDPN"]: ax[1].annotate(r.gene,(r.mean_beta,r.mean_beta_adj),fontsize=7.5,xytext=(3,2),textcoords="offset points",fontweight="bold" if r.gene=="ANGPTL4" else "normal")
lim=[-.55,.8]; ax[1].plot(lim,lim,ls="--",c="grey",lw=.8); ax[1].axhline(0,c="k",lw=.5); ax[1].axvline(0,c="k",lw=.5)
ax[1].set_xlabel("Effect, admixture-adjusted"); ax[1].set_ylabel("Effect, additionally adjusted for malignant state mix")
ax[1].set_title("B  How much of the niche-6 signal is beyond MES-like state?\n(47 NCDEGs; dashed = no change after state adjustment)",loc="left",fontsize=10)
plt.tight_layout(); plt.savefig("ncdeg_out/ncdeg_niche6_malignant.png",dpi=200)
E=pd.read_csv("ncdeg_out/ncdeg_metaprogram_enrichment.csv"); E=E[E.fdr<0.05].sort_values("fdr")
E.to_csv("ncdeg_out/ncdeg_metaprogram_enrichment_sig.csv",index=False); print(E[["niche","direction","meta_program","overlap","fdr"]].to_string())
