"""14b_progression_scores.py -- GSE318238 patient BT2015050 (IDH-mut astrocytoma; 2015/2017/2019, Visium CytAssist).
Case study only: copy-number compartments fail (probe-based platform vs 3' cortex reference -> 100% 'malignant') and the
IDH-wt GBmap reference under-calls IDH-mut malignant cells, so niches are NOT projected. Instead: per-spot niche-signature
and meta-program scores (log1p CP10k, z across all spots of the three sections), compared between time points."""
import json, numpy as np, pandas as pd, scanpy as sc, anndata as an, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.stats import kruskal
S = json.load(open("sig_out/bulk_scoring_sets.json")); single = S.pop("_single")
A = an.concat([sc.read_h5ad(f"prog/qc/IDHm_GSE318238_{y}.h5ad") for y in ["2015", "2017", "2019"]], join="inner", index_unique="-")
sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A)
X = pd.DataFrame(A.X.toarray() if hasattr(A.X, "toarray") else A.X, index=A.obs_names, columns=A.var_names)
Z = (X - X.mean()) / X.std().replace(0, np.nan)
sc_ = pd.DataFrame({k: Z[[g for g in v if g in Z.columns]].mean(1) for k, v in S.items()})
cov = {k: f"{sum(g in Z.columns for g in v)}/{len(v)}" for k, v in S.items()}
for g in ["ANGPTL4", "VEGFA", "CA12", "NAMPT", "CD44"]:
    if g in X: sc_[g] = X[g]
sc_["year"] = A.obs["sample"].str[-4:].values
sc_[["array_row", "array_col"]] = A.obs[["array_row", "array_col"]].values
sc_.to_csv("prog/progression_spot_scores.csv.gz")
cols = ["HypoxicMES_N6", "Myeloid_N3", "Neuronal_N1", "Oligodendroglial_N10", "MP_Hypoxia", "MP_MES", "MP_AC", "MP_NPC", "MP_OPC", "ANGPTL4"]
T = sc_.groupby("year")[cols].mean().round(3)
T.loc["kruskal_p"] = [kruskal(*[sc_.loc[sc_.year == y, c] for y in ["2015", "2017", "2019"]]).pvalue for c in cols]
hi = sc_.HypoxicMES_N6 > sc_.HypoxicMES_N6.quantile(0.9)
T["pct_spots_top10pct_N6"] = list((hi.groupby(sc_.year).mean() * 100).round(1)) + [np.nan]
T.to_csv("prog/progression_score_summary.csv"); print(T.to_string()); print(cov)
D = pd.concat([pd.read_csv(f"prog/deconv/spots/IDHm_GSE318238_{y}_proportions.csv.gz", index_col=0).select_dtypes("number").drop(columns="fit_r").mean().rename(y) for y in ["2015", "2017", "2019"]], axis=1).round(3)
D.to_csv("prog/progression_mean_celltype_fractions.csv")
fig, ax = plt.subplots(3, 4, figsize=(15, 10.5))
for j, c in enumerate(["HypoxicMES_N6", "Myeloid_N3", "Neuronal_N1", "Oligodendroglial_N10"]):
    lim = np.nanpercentile(np.abs(sc_[c]), 98)
    for i, y in enumerate(["2015", "2017", "2019"]):
        d = sc_[sc_.year == y]; a = ax[i, j]
        im = a.scatter(d.array_col, -d.array_row, c=d[c], cmap="RdBu_r", vmin=-lim, vmax=lim, s=6 if y == "2017" else 14, lw=0)
        a.set_aspect("equal"); a.set_xticks([]); a.set_yticks([])
        a.set_title(f"{y}: {c} (mean {d[c].mean():+.2f})", fontsize=9)
plt.suptitle("GSE318238 patient BT2015050 (IDH-mut astrocytoma, 2015 → 2017 → 2019): niche-signature scores per spot", fontsize=11)
plt.tight_layout(); plt.savefig("prog/progression_niche_score_maps.png", dpi=170)
