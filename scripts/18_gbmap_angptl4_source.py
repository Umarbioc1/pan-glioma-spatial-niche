"""18 - Single-cell source of ANGPTL4 and expression of its receptors in GBmap (core atlas reference subsample).

Input : gbmap_ref.h5ad (up to 400 cells per level-3 type from the GBmap core atlas, raw counts, 15 studies)
Output: per level-4 state mean log-normalised expression and % expressing (cells and study-balanced),
        MES-like hypoxia/MHC vs hypoxia-independent comparison (Mann-Whitney on cells; per-study direction),
        Figure S6 dot plot.
"""
import sys, os, numpy as np, pandas as pd, anndata as an, scanpy as sc, scipy.sparse as sp
from scipy.stats import mannwhitneyu, wilcoxon
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.linewidth": 0.6, "pdf.fonttype": 42})
REF = sys.argv[1] if len(sys.argv) > 1 else "/root/gbmap_ref.h5ad"
OUT = sys.argv[2] if len(sys.argv) > 2 else "gbmap_out"; FIG = sys.argv[3] if len(sys.argv) > 3 else "manuscript_figs"
os.makedirs(OUT, exist_ok=True)
R = an.read_h5ad(REF); R.var_names_make_unique()
R.layers["counts"] = R.X.copy()
sc.pp.normalize_total(R, target_sum=1e4); sc.pp.log1p(R)
LIG = ["ANGPTL4"]; REC = ["SDC2", "SDC3", "SDC4", "ITGA5", "ITGB1", "CD44", "CDH11"]
CTRL = ["VEGFA", "ADM", "SPP1"]; TIER = ["CAV1", "VMP1", "NRN1", "NDRG1", "NAMPT", "ARRDC3", "PDPN", "CA12"]
genes = [g for g in LIG + REC + CTRL + TIER if g in R.var_names]
L4 = R.obs.annotation_level_4.astype(str); L3 = R.obs.annotation_level_3.astype(str)
X = R[:, genes].X; X = X.toarray() if sp.issparse(X) else np.asarray(X)
C = R[:, genes].layers["counts"]; C = C.toarray() if sp.issparse(C) else np.asarray(C)
E = pd.DataFrame(X, columns=genes, index=R.obs_names); P = pd.DataFrame(C > 0, columns=genes, index=R.obs_names)
grp = L4.where(L4.map(L4.value_counts()) >= 25, L3)          # small level-4 groups fall back to level 3
mean = E.groupby(grp.values).mean(); pct = P.groupby(grp.values).mean() * 100; n = grp.value_counts()
# study-balanced mean (average of per-study means, studies with >=10 cells in the group)
sb = {}
for gname, idx in grp.groupby(grp).groups.items():
    a = R.obs.loc[idx, "author"]; ok = a.value_counts(); ok = ok[ok >= 10].index
    sb[gname] = E.loc[idx][a.isin(ok).values].groupby(a[a.isin(ok)].values).mean().mean() if len(ok) else E.loc[idx].mean()
sbm = pd.DataFrame(sb).T
mean.round(4).to_csv(f"{OUT}/gbmap_state_mean_lognorm.csv"); pct.round(2).to_csv(f"{OUT}/gbmap_state_pct_expressing.csv")
sbm.round(4).to_csv(f"{OUT}/gbmap_state_mean_study_balanced.csv"); n.to_csv(f"{OUT}/gbmap_state_n_cells.csv")
# MES-like hypoxia/MHC vs hypoxia-independent
a_, b_ = grp == "MES-like hypoxia/MHC", grp == "MES-like hypoxia independent"
rows = []
for g in genes:
    p = mannwhitneyu(E.loc[a_.values, g], E.loc[b_.values, g]).pvalue
    au = R.obs.author
    per = [(E.loc[(a_ & (au == s)).values, g].mean(), E.loc[(b_ & (au == s)).values, g].mean())
           for s in au.unique() if (a_ & (au == s)).sum() >= 10 and (b_ & (au == s)).sum() >= 10]
    rows.append(dict(gene=g, mean_hypoxiaMHC=E.loc[a_.values, g].mean(), mean_hypoxia_indep=E.loc[b_.values, g].mean(),
                     pct_hypoxiaMHC=P.loc[a_.values, g].mean() * 100, pct_hypoxia_indep=P.loc[b_.values, g].mean() * 100,
                     p_MWU=p, n_studies_both=len(per), n_studies_higher_in_hypoxia=sum(x > y for x, y in per)))
M = pd.DataFrame(rows); M.to_csv(f"{OUT}/gbmap_MES_hypoxia_vs_independent.csv", index=False)
# rank of each group for ANGPTL4
top = mean["ANGPTL4"].sort_values(ascending=False)
pd.set_option("display.width", 250)
print("ANGPTL4 by state (mean, % expressing, n):"); print(pd.concat([top, pct["ANGPTL4"], n], axis=1).head(12).round(3).to_string())
print(M.round(4).to_string())
for g in REC: print(g, "top:", mean[g].sort_values(ascending=False).head(4).round(2).to_dict())

# ---- Figure S6 dot plot ----
order = ["MES-like hypoxia/MHC", "MES-like hypoxia independent", "AC-like", "NPC-like OPC", "NPC-like neural", "OPC-like",
         "TAM-BDM hypoxia/MES", "TAM-BDM anti-infl", "TAM-BDM INF", "TAM-MG aging sig", "TAM-MG pro-infl I", "Mono hypoxia", "Mono naive",
         "Endo capilar", "Endo arterial", "Tip-like", "Pericyte", "SMC", "Oligodendrocyte", "OPC", "Astrocyte", "CD4/CD8"]
order = [o for o in order if o in mean.index]
gl = LIG + REC + CTRL + TIER; gl = [g for g in gl if g in genes]
Z = mean.loc[order, gl]; Zs = (Z - Z.min()) / (Z.max() - Z.min() + 1e-9)
fig, ax = plt.subplots(figsize=(7.2, 5.0))
yy, xx = np.meshgrid(range(len(order)), range(len(gl)), indexing="ij")
sca = ax.scatter(xx.ravel(), yy.ravel(), s=(pct.loc[order, gl].values.ravel() / 100) * 110 + 1,
                 c=Zs.values.ravel(), cmap="Reds", vmin=0, vmax=1, edgecolors="0.6", linewidths=0.3)
ax.set_xticks(range(len(gl))); ax.set_xticklabels(gl, rotation=60, ha="right")
ax.set_yticks(range(len(order))); ax.set_yticklabels([f"{o} ({n[o]})" for o in order]); ax.invert_yaxis()
for x in [0.5, 0.5 + len(REC), 0.5 + len(REC) + len(CTRL)]: ax.axvline(x, color="0.8", lw=0.6)
for y in [5.5, 12.5, 17.5]: ax.axhline(y, color="0.8", lw=0.6)
for x, t in [(0, "ligand"), (1 + (len(REC) - 1) / 2, "ANGPTL4 receptors"), (1 + len(REC) + 1, "hypoxia/MES controls"),
             (1 + len(REC) + len(CTRL) + (len(TIER) - 1) / 2, "tier-1 targets")]:
    ax.text(x, -1.1, t, ha="center", fontsize=6.5, style="italic")
ax.set_xlim(-0.6, len(gl) - 0.4); ax.set_ylim(len(order) - 0.4, -1.6)
cb = fig.colorbar(sca, ax=ax, fraction=0.025, pad=0.01); cb.set_label("mean expression (scaled per gene)", fontsize=6)
for pv in [10, 40, 80]: ax.scatter([], [], s=pv / 100 * 110 + 1, c="w", edgecolors="0.4", label=f"{pv}%")
ax.legend(title="% expressing", loc="upper left", bbox_to_anchor=(1.12, 1.0), frameon=False, fontsize=6, title_fontsize=6.5)
fig.savefig(f"{FIG}/Figure_S7.png", dpi=300, bbox_inches="tight"); fig.savefig(f"{FIG}/Figure_S7.pdf", bbox_inches="tight")
print("fig ok")
