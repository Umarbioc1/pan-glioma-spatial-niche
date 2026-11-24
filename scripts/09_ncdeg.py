"""
09_ncdeg.py  --  Aim 4: niche-dependent, cell-type-resolved DEGs (NCDEGs)
C-SIDE-like approximation in Python (Li et al. 2026 used C-SIDE/spacexr in R).

For focal cell group F (Malignant = AC/MES/NPC/OPC-like; Myeloid = TAM-MG/TAM-BDM/Mono/DC):
  * focal spots: deconvolved fraction of F >= threshold (Malignant 0.5 + CNV-malignant; Myeloid 0.25)
  * genes: expressed by F in the GBmap reference (mean CP10k >= 0.5) and F-specific enough
           (mean in F >= 0.5 x max mean of any other lineage group) -> limits contamination
  * per section, per niche (>=15 focal spots in, >=15 out):
        log1p(CP10k) ~ 1 + I(niche) + other-lineage fractions + log(UMI)   (OLS, all genes at once)
    -> beta = niche effect within F-rich spots, adjusted for admixture
  * meta-analysis with SECTION as the unit (spots are spatially autocorrelated):
        mean beta, one-sample t-test across sections, sign consistency -> BH FDR per niche
  * annotation: hypergeometric overlap of NCDEG lists with Gavish 2023 meta-programs

usage: python3 09_ncdeg.py <spot_LCP_niche.csv.gz> <deconv spots dir> <h5ad dirs comma> <gbmap_ref.h5ad> <gavish_MPs.json> <outdir>
"""
import sys, os, json, numpy as np, pandas as pd, anndata as an, scanpy as sc, scipy.sparse as sp
from scipy.stats import ttest_1samp, hypergeom
from statsmodels.stats.multitest import multipletests

spot_f, dec_dir, h5dirs, ref_f, mp_f, out = sys.argv[1:7]
os.makedirs(out, exist_ok=True)
GROUPS = {"Malignant": ["AC-like", "MES-like", "NPC-like", "OPC-like"],
          "Myeloid": ["TAM-MG", "TAM-BDM", "Mono", "DC"],
          "Lymphoid": ["CD4/CD8", "NK", "B cell", "Plasma B", "Mast"],
          "Vascular": ["Endothelial", "Mural cell"],
          "Glia": ["Oligodendrocyte", "OPC", "Astrocyte", "RG"],
          "Neuron": ["Neuron"]}
FOCAL = {"Malignant": 0.5} if os.environ.get("STATE_ADJ") == "1" else {"Malignant": 0.5, "Myeloid": 0.25}

# ---- reference: lineage-group means -> gene filter per focal group ----
R = an.read_h5ad(ref_f); R.var_names_make_unique()
sc.pp.normalize_total(R, target_sum=1e4)
lvl = R.obs.annotation_level_3.astype(str).values
grp_of = {t: g for g, ts in GROUPS.items() for t in ts}
gl = np.array([grp_of.get(t, "Other") for t in lvl])
G = pd.DataFrame({g: np.asarray(R[gl == g].X.mean(0)).ravel() for g in GROUPS}, index=R.var_names)
drop = G.index.str.match(r"^(MT-|RPL|RPS|ENSG|MALAT1|NEAT1)")
genes_for = {}
for f in FOCAL:
    oth = G.drop(columns=f).max(1)
    keep = (G[f] >= 0.5) & (G[f] >= 0.5 * oth) & ~drop
    genes_for[f] = G.index[keep].tolist(); print(f, "genes tested:", len(genes_for[f]))
G.round(3).to_csv(f"{out}/reference_lineage_means.csv")

S = pd.read_csv(spot_f)
h5 = {f[:-5]: os.path.join(d, f) for d in h5dirs.split(",") for f in os.listdir(d) if f.endswith(".h5ad")}
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
    for f, thr in FOCAL.items():
        m = (frac[f] >= thr).values.copy()
        if f == "Malignant": m &= (g.compartment2 == "Malignant").values
        if m.sum() < 40: continue
        gs = [x for x in genes_for[f] if x in A.var_names]
        Y = A[m, gs].X; Y = Y.toarray() if sp.issparse(Y) else np.asarray(Y)
        cov = frac.loc[m].drop(columns=[f, "Lymphoid"]).values
        cov = np.c_[cov, np.log(umi[m])]
        if f == "Malignant" and os.environ.get("STATE_ADJ") == "1":
            mal = P.loc[m, GROUPS["Malignant"]].values; mal = mal / np.maximum(mal.sum(1, keepdims=True), 1e-9)
            cov = np.c_[cov, mal[:, 1:]]  # MES/NPC/OPC-like share within malignant (AC-like = reference)
        niche = g.niche.values[m]
        for n in np.unique(niche):
            ind = (niche == n).astype(float)
            if ind.sum() < 15 or (1 - ind).sum() < 15: continue
            X = np.c_[np.ones(len(ind)), ind, cov]
            XtX_inv = np.linalg.pinv(X.T @ X)
            B = XtX_inv @ X.T @ Y
            res = Y - X @ B; dfree = X.shape[0] - X.shape[1]
            s2 = (res ** 2).sum(0) / dfree
            se = np.sqrt(s2 * XtX_inv[1, 1]) + 1e-9
            rows.append(pd.DataFrame({"focal": f, "sample": s, "cohort": g.cohort.iloc[0], "subtype": g.subtype.iloc[0],
                                      "niche": n, "gene": gs, "beta": B[1], "t": B[1] / se,
                                      "n_in": int(ind.sum()), "n_out": int((1 - ind).sum()),
                                      "mean_in": Y[ind == 1].mean(0), "mean_out": Y[ind == 0].mean(0)}))
    print(s, "done", flush=True)

D = pd.concat(rows, ignore_index=True)
D.to_csv(f"{out}/ncdeg_per_section.csv.gz", index=False)

def meta(x):
    b = x.beta.values; k = len(b)
    p = ttest_1samp(b, 0).pvalue if k >= 3 else np.nan
    return pd.Series({"n_sections": k, "n_GBM": x.subtype.astype(str).str.contains("GBM").sum(),
                      "n_greenwald": (x.cohort == "Greenwald2024").sum(), "mean_beta": b.mean(), "median_beta": np.median(b),
                      "sign_consistency": max((b > 0).mean(), (b < 0).mean()), "p": p,
                      "mean_expr_in": x.mean_in.mean(), "mean_expr_out": x.mean_out.mean()})
M = D.groupby(["focal", "niche", "gene"]).apply(meta).reset_index()
M = M[M.n_sections >= 3].copy()
M["fdr"] = np.nan
for _, idx in M.groupby(["focal", "niche"]).groups.items():
    M.loc[idx, "fdr"] = multipletests(M.loc[idx, "p"].fillna(1), method="fdr_bh")[1]
M["direction"] = np.where(M.mean_beta > 0, "up", "down")
M["NCDEG"] = (M.fdr < 0.05) & (M.mean_beta.abs() >= 0.2) & (M.sign_consistency >= 0.75)
M = M.sort_values(["focal", "niche", "NCDEG", "mean_beta"], ascending=[True, True, False, False])
M.to_csv(f"{out}/ncdeg_meta.csv", index=False)
sig = M[M.NCDEG]; sig.to_csv(f"{out}/ncdeg_significant.csv", index=False)
print(sig.groupby(["focal", "niche", "direction"]).size().unstack(fill_value=0))

# ---- meta-program enrichment ----
MP = json.load(open(mp_f)); er = []
for (f, n, d), x in sig.groupby(["focal", "niche", "direction"]):
    U = set(M[(M.focal == f) & (M.niche == n)].gene); hits = set(x.gene)
    for name, gs in MP.items():
        gs = set(gs) & U
        if len(gs) < 5: continue
        k = len(hits & gs)
        if k < 3: continue
        p = hypergeom.sf(k - 1, len(U), len(gs), len(hits))
        er.append({"focal": f, "niche": n, "direction": d, "meta_program": name, "overlap": k, "mp_size_tested": len(gs),
                   "n_hits": len(hits), "p": p, "genes": ",".join(sorted(hits & gs))})
E = pd.DataFrame(er)
if len(E):
    E["fdr"] = multipletests(E.p, method="fdr_bh")[1]
    E = E.sort_values("p"); E.to_csv(f"{out}/ncdeg_metaprogram_enrichment.csv", index=False)
    print(E[E.fdr < 0.05].groupby(["focal", "niche", "direction"]).head(3)[["focal", "niche", "direction", "meta_program", "overlap", "fdr"]].to_string())
