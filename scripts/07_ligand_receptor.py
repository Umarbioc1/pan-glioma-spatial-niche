"""
07_ligand_receptor.py  --  Aim 3: niche-enriched ligand-receptor (LR) interactions
Adapted from Li et al. 2026 (Cell Rep Med) pan-cancer framework for glioma Visium.

Per section:
  * log-normalised expression (1e4), CellChatDB pairs (liana omni_resource, 'cellchatdb')
  * complexes: subunit expression = min over subunits (conservative)
  * spot score = 0.5*[sqrt(L_spot * R_nbr) + sqrt(R_spot * L_nbr)],
    nbr = mean over spot + first ring of neighbours (Visium hex, ~100-200 um)
  * per niche (>=20 spots in, >=20 out): Mann-Whitney niche vs rest of section -> AUC, z
Meta-analysis across sections (rank-based, as in Li et al.):
  * per section, pairs ranked by AUC within niche -> percentile rank
  * combine: mean percentile, median AUC, fraction of sections with AUC>0.6,
    Stouffer z -> p -> BH FDR
Sender/receiver annotation: mean expression of ligand/receptor per GBmap level-3 type.

usage: python3 07_ligand_receptor.py <spot_LCP_niche.csv.gz> <h5ad_dirs comma> <cellchatdb_lr.csv> <gbmap_ref.h5ad> <outdir>
"""
import sys, os, numpy as np, pandas as pd, anndata as an, scanpy as sc
from scipy.spatial import cKDTree
from scipy.stats import mannwhitneyu, norm, rankdata
import scipy.sparse as sp

spot_f, h5dirs, lr_f, ref_f, out = sys.argv[1:6]
os.makedirs(out, exist_ok=True)
S = pd.read_csv(spot_f)
LR = pd.read_csv(lr_f).drop_duplicates()
LR["L_sub"] = LR.ligand.str.split("_"); LR["R_sub"] = LR.receptor.str.split("_")
LR["pair"] = LR.ligand + "->" + LR.receptor
LR = LR.drop_duplicates("pair").reset_index(drop=True)
lr_genes = sorted(set(g for l in LR.L_sub for g in l) | set(g for l in LR.R_sub for g in l))

h5 = {}
for d in h5dirs.split(","):
    for f in os.listdir(d):
        if f.endswith(".h5ad"): h5[f[:-5]] = os.path.join(d, f)

def complex_expr(E, genes_idx, subs):
    """E spots x genes (dense); returns spots x pairs (min over subunits); NaN col if missing"""
    M = np.full((E.shape[0], len(subs)), np.nan, dtype=np.float32)
    for j, s in enumerate(subs):
        idx = [genes_idx.get(g) for g in s]
        if any(i is None for i in idx): continue
        M[:, j] = E[:, idx].min(1)
    return M

rows = []
for s, g in S.groupby("sample"):
    if s not in h5: print("no h5ad", s); continue
    A = an.read_h5ad(h5[s]); A.var_names_make_unique()
    g = g.set_index("spot"); common = A.obs_names.intersection(g.index)
    A = A[common].copy(); g = g.loc[common]
    sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A)
    keep = [x for x in lr_genes if x in A.var_names]
    X = A[:, keep].X; X = X.toarray() if sp.issparse(X) else np.asarray(X)
    gi = {x: i for i, x in enumerate(keep)}
    L = complex_expr(X, gi, LR.L_sub); R = complex_expr(X, gi, LR.R_sub)
    ok = ~(np.isnan(L).any(0) | np.isnan(R).any(0))
    L, R = L[:, ok], R[:, ok]; pairs = LR.pair.values[ok]
    # neighbours: first ring (hex: coords col, row*sqrt3; adjacent distance = 2)
    xy = np.c_[g.array_col.values, g.array_row.values * np.sqrt(3)]
    T = cKDTree(xy); nb = T.query_ball_point(xy, r=2.1)
    W = sp.lil_matrix((len(xy), len(xy)), dtype=np.float32)
    for i, n in enumerate(nb): W[i, n] = 1.0 / len(n)
    W = W.tocsr()
    Ln, Rn = W @ L, W @ R
    Sc = 0.5 * (np.sqrt(L * Rn) + np.sqrt(R * Ln))
    niche = g.niche.values
    for n in np.unique(niche[~pd.isna(niche)]):
        m = niche == n; n1, n2 = m.sum(), (~m).sum()
        if n1 < 20 or n2 < 20: continue
        a, b = Sc[m], Sc[~m]
        expr_frac = (a > 0).mean(0)
        U, p = mannwhitneyu(a, b, axis=0, alternative="two-sided", method="asymptotic")
        auc = U / (n1 * n2)
        # signed z from two-sided p
        z = norm.isf(np.clip(p, 1e-300, 1) / 2) * np.sign(auc - 0.5)
        pct = rankdata(np.nan_to_num(auc, nan=0.5)) / len(auc)
        rows.append(pd.DataFrame({"sample": s, "cohort": g.cohort.iloc[0], "subtype": g.subtype.iloc[0],
                                  "niche": str(n), "pair": pairs, "n_in": n1, "auc": auc, "z": z,
                                  "pct_rank": pct, "frac_in": expr_frac,
                                  "mean_in": a.mean(0), "mean_out": b.mean(0)}))
    print(s, A.n_obs, "spots,", ok.sum(), "pairs", flush=True)

D = pd.concat(rows, ignore_index=True)
D.to_csv(f"{out}/lr_per_section.csv.gz", index=False)

# ---- meta-analysis ----
D = D[D.frac_in >= 0.10]  # expressed in >=10% of niche spots in that section
def meta(x):
    k = len(x)
    zs = x.z.replace([np.inf, -np.inf], np.nan).dropna()
    Z = zs.sum() / np.sqrt(len(zs)) if len(zs) else np.nan
    return pd.Series({"n_sections": k, "n_greenwald": (x.cohort == "Greenwald2024").sum(),
                      "n_ravi": (x.cohort != "Greenwald2024").sum(),
                      "median_auc": x.auc.median(), "frac_auc_gt0.6": (x.auc > 0.6).mean(),
                      "mean_pct_rank": x.pct_rank.mean(), "stouffer_z": Z,
                      "log2FC_score": np.log2((x.mean_in.mean() + 1e-3) / (x.mean_out.mean() + 1e-3))})
M = D.groupby(["niche", "pair"]).apply(meta).reset_index()
M = M[M.n_sections >= 3].copy()
M["p"] = norm.sf(M.stouffer_z)  # one-sided: enriched in niche
from statsmodels.stats.multitest import multipletests
M["fdr"] = np.nan
for n, idx in M.groupby("niche").groups.items():
    M.loc[idx, "fdr"] = multipletests(M.loc[idx, "p"].fillna(1), method="fdr_bh")[1]
M["significant"] = (M.fdr < 0.05) & (M.median_auc > 0.6) & (M["frac_auc_gt0.6"] >= 0.6)

# ---- sender / receiver annotation from GBmap reference ----
Rf = an.read_h5ad(ref_f); Rf.var_names_make_unique()
sc.pp.normalize_total(Rf, target_sum=1e4); sc.pp.log1p(Rf)
keep = [x for x in lr_genes if x in Rf.var_names]
Xr = Rf[:, keep].X; Xr = Xr.toarray() if sp.issparse(Xr) else np.asarray(Xr)
types = Rf.obs.annotation_level_3.astype(str).values
TM = pd.DataFrame(Xr, columns=keep).groupby(types).mean()
gi = {x: i for i, x in enumerate(keep)}
def top_type(subs):
    if any(gg not in gi for gg in subs): return ("", np.nan)
    v = TM[subs].min(1); return (v.idxmax(), round(float(v.max()), 3))
LR2 = LR.set_index("pair")
ann = {}
for pr in M.pair.unique():
    lt, lv = top_type(LR2.loc[pr, "L_sub"]); rt, rv = top_type(LR2.loc[pr, "R_sub"])
    ann[pr] = (lt, lv, rt, rv)
M[["sender_top", "sender_expr", "receiver_top", "receiver_expr"]] = pd.DataFrame(
    [ann[p] for p in M.pair], index=M.index)
M = M.sort_values(["niche", "significant", "mean_pct_rank"], ascending=[True, False, False])
M.to_csv(f"{out}/lr_niche_meta.csv", index=False)
TM.T.to_csv(f"{out}/gbmap_LR_gene_mean_by_type.csv")
sig = M[M.significant]
print(sig.groupby("niche").size())
sig.to_csv(f"{out}/lr_niche_significant.csv", index=False)
