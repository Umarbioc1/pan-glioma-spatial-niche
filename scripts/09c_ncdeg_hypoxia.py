"""
09c_ncdeg_hypoxia.py  --  Hypoxia-adjusted niche-dependent genes in malignant spots (sensitivity analysis)

Same design as 09_ncdeg.py (malignant focal spots: deconvolved malignant fraction >= 0.5 and CNV-malignant;
OLS per section on niche membership + other lineage fractions + log UMI; section = unit of replication),
with two extra models for niche 6 versus the rest of each section:
  H  : + per-spot hypoxia score (Gavish 'Cancer|MP6 Hypoxia', mean z of log-normalised expression within section)
  HS : + hypoxia score + malignant-state shares (MES/NPC/OPC-like within malignant; AC-like reference)
For genes that belong to the hypoxia program (e.g. ANGPTL4, NDRG1, NRN1, CAV1, CA12, VEGFA) the hypoxia score
is recomputed leaving that gene out, so no gene is adjusted for itself.
Note: hypoxia is part of what defines niche 6, so H/HS are deliberately conservative (collinear) tests.

usage: python3 09c_ncdeg_hypoxia.py <spot_LCP_niche.csv.gz> <deconv spots dir> <h5ad dirs comma> <reference_lineage_means.csv> <gavish_MPs.json> <outdir> [niche=6]
"""
import sys, os, json, numpy as np, pandas as pd, anndata as an, scanpy as sc, scipy.sparse as sp
from scipy.stats import ttest_1samp
from statsmodels.stats.multitest import multipletests

spot_f, dec_dir, h5dirs, gmeans_f, mp_f, out = sys.argv[1:7]
NICHE = "Niche_" + (sys.argv[7] if len(sys.argv) > 7 else "6")
os.makedirs(out, exist_ok=True)
GROUPS = {"Malignant": ["AC-like", "MES-like", "NPC-like", "OPC-like"],
          "Myeloid": ["TAM-MG", "TAM-BDM", "Mono", "DC"],
          "Lymphoid": ["CD4/CD8", "NK", "B cell", "Plasma B", "Mast"],
          "Vascular": ["Endothelial", "Mural cell"],
          "Glia": ["Oligodendrocyte", "OPC", "Astrocyte", "RG"],
          "Neuron": ["Neuron"]}
G = pd.read_csv(gmeans_f, index_col=0)
oth = G.drop(columns="Malignant").max(1)
drop = G.index.str.match(r"^(MT-|RPL|RPS|ENSG|MALAT1|NEAT1)")
genes = G.index[(G.Malignant >= 0.5) & (G.Malignant >= 0.5 * oth) & ~drop].tolist()
HYP = json.load(open(mp_f))["Cancer|MP6 Hypoxia"]

S = pd.read_csv(spot_f); S["niche"] = S.niche.astype(str)
h5 = {f[:-5]: os.path.join(d, f) for d in h5dirs.split(",") for f in os.listdir(d) if f.endswith(".h5ad")}

def ols_beta(X, Y):
    XtX_inv = np.linalg.pinv(X.T @ X); B = XtX_inv @ X.T @ Y
    res = Y - X @ B; s2 = (res ** 2).sum(0) / (X.shape[0] - X.shape[1])
    return B[1], B[1] / (np.sqrt(s2 * XtX_inv[1, 1]) + 1e-9)

rows = []
for s, g in S.groupby("sample"):
    pf = f"{dec_dir}/{s}_proportions.csv.gz"
    if s not in h5 or not os.path.exists(pf): continue
    P = pd.read_csv(pf, index_col=0)
    A = an.read_h5ad(h5[s]); A.var_names_make_unique()
    g = g.set_index("spot"); com = A.obs_names.intersection(g.index).intersection(P.index)
    A, g, P = A[com].copy(), g.loc[com], P.loc[com]
    umi = np.asarray(A.X.sum(1)).ravel()
    sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A)
    frac = pd.DataFrame({k: P[v].sum(1) for k, v in GROUPS.items()}, index=com)
    m = ((frac.Malignant >= 0.5) & (g.compartment2 == "Malignant")).values
    if m.sum() < 40: continue
    ind = (g.niche.values[m] == NICHE).astype(float)
    if ind.sum() < 15 or (1 - ind).sum() < 15: continue
    # hypoxia program z-scores within section (all spots), then restricted to focal spots
    hg = [x for x in HYP if x in A.var_names]
    H = A[:, hg].X; H = H.toarray() if sp.issparse(H) else np.asarray(H)
    Hz = (H - H.mean(0)) / (H.std(0) + 1e-9); Hz = Hz[m]
    hsum, hn = Hz.sum(1), len(hg)
    gs = [x for x in genes if x in A.var_names]
    Y = A[m, gs].X; Y = Y.toarray() if sp.issparse(Y) else np.asarray(Y)
    base = np.c_[frac.loc[m].drop(columns=["Malignant", "Lymphoid"]).values, np.log(umi[m])]
    mal = P.loc[m, GROUPS["Malignant"]].values; mal = mal / np.maximum(mal.sum(1, keepdims=True), 1e-9)
    state = mal[:, 1:]
    one = np.ones(len(ind))
    out_b = {}
    X0 = np.c_[one, ind, base]; out_b["unadj"] = ols_beta(X0, Y)
    hyp_all = hsum / hn
    for key, extra in [("H", []), ("HS", [state])]:
        Xh = np.c_[one, ind, base, hyp_all, *extra]
        b, t = ols_beta(Xh, Y)
        # leave-one-out hypoxia score for genes in the program
        for j, gg in enumerate(gs):
            if gg in hg:
                k = hg.index(gg); loo = (hsum - Hz[:, k]) / (hn - 1)
                bj, tj = ols_beta(np.c_[one, ind, base, loo, *extra], Y[:, [j]])
                b[j], t[j] = bj[0], tj[0]
        out_b[key] = (b, t)
    rows.append(pd.DataFrame({"sample": s, "cohort": g.cohort.iloc[0], "subtype": g.subtype.iloc[0], "gene": gs,
                              "n_in": int(ind.sum()), "n_out": int((1 - ind).sum()),
                              "r_niche_hypoxia": np.corrcoef(ind, hyp_all)[0, 1],
                              **{f"beta_{k}": v[0] for k, v in out_b.items()}}))
    print(s, int(ind.sum()), "niche spots; r(niche, hypoxia) =", round(np.corrcoef(ind, hyp_all)[0, 1], 2), flush=True)

D = pd.concat(rows, ignore_index=True)
D.to_csv(f"{out}/ncdeg_hypoxia_per_section.csv.gz", index=False)
res = []
for key in ["unadj", "H", "HS"]:
    def meta(x, c=f"beta_{key}"):
        b = x[c].values
        return pd.Series({"n_sections": len(b), "mean_beta": b.mean(), "p": ttest_1samp(b, 0).pvalue if len(b) >= 3 else np.nan,
                          "sign_consistency": max((b > 0).mean(), (b < 0).mean()), "n_up": int((b > 0).sum())})
    M = D.groupby("gene").apply(meta).reset_index(); M = M[M.n_sections >= 3]
    M["fdr"] = multipletests(M.p.fillna(1), method="fdr_bh")[1]; M["model"] = key
    M["NCDEG"] = (M.fdr < 0.05) & (M.mean_beta.abs() >= 0.2) & (M.sign_consistency >= 0.75)
    res.append(M)
R = pd.concat(res)
W = R.pivot(index="gene", columns="model", values=["mean_beta", "p", "fdr", "n_up", "n_sections"])
W.columns = [f"{a}_{b}" for a, b in W.columns]
W = W.sort_values("p_HS")
W.to_csv(f"{out}/ncdeg_hypoxia_adjusted_niche{NICHE}.csv")
print("NCDEGs per model:", R.groupby("model").NCDEG.sum().to_dict())
print("sections:", D["sample"].nunique(), "median r(niche, hypoxia):", round(D.groupby("sample").r_niche_hypoxia.first().median(), 2))
focus = ["ANGPTL4", "CAV1", "VMP1", "NRN1", "NDRG1", "NAMPT", "ARRDC3", "PDPN", "CA12", "CAMK2N1", "VEGFA", "ADM", "CD44", "CHI3L1"]
pd.set_option("display.width", 250)
print(W.loc[[x for x in focus if x in W.index], [c for c in W.columns if c.startswith(("mean_beta", "p_", "n_up"))]].round(4).to_string())
print(W.head(15)[["mean_beta_HS", "p_HS", "fdr_HS", "n_up_HS", "n_sections_HS"]].round(4).to_string())
