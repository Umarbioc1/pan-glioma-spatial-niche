"""
10_niche_signatures.py -- niche gene signatures for bulk scoring (GliomaNicheScore, Aim 5)
Per section: pseudobulk log1p CP10k per niche vs rest of section (>=30 spots each).
Per niche: mean logFC across sections, fraction of sections with logFC>0.25, t-test across sections.
Signature = top 50 genes by mean logFC with FDR<0.05, sign consistency >=0.75, logFC>0.3,
and not in the top-200 of more than 2 niches (specificity).
"""
import sys, os, numpy as np, pandas as pd, anndata as an, scanpy as sc, scipy.sparse as sp
from scipy.stats import ttest_1samp
from statsmodels.stats.multitest import multipletests
spot_f, h5dirs, out = sys.argv[1:4]; os.makedirs(out, exist_ok=True)
S = pd.read_csv(spot_f)
h5 = {f[:-5]: os.path.join(d, f) for d in h5dirs.split(",") for f in os.listdir(d) if f.endswith(".h5ad")}
rows = []
for s, g in S.groupby("sample"):
    if s not in h5: continue
    A = an.read_h5ad(h5[s]); A.var_names_make_unique()
    A = A[:, ~A.var_names.str.match(r"^(MT-|RPL|RPS|ENSG|MALAT1)")].copy()
    g = g.set_index("spot"); com = A.obs_names.intersection(g.index); A, g = A[com].copy(), g.loc[com]
    X = A.X.tocsr() if sp.issparse(A.X) else sp.csr_matrix(A.X)
    for n in g.niche.unique():
        m = (g.niche == n).values
        if m.sum() < 30 or (~m).sum() < 30: continue
        a = np.asarray(X[m].sum(0)).ravel(); b = np.asarray(X[~m].sum(0)).ravel()
        a = np.log1p(a / a.sum() * 1e4); b = np.log1p(b / b.sum() * 1e4)
        rows.append(pd.DataFrame({"sample": s, "niche": n, "gene": A.var_names, "lfc": a - b, "expr_in": a}))
    print(s, flush=True)
D = pd.concat(rows); D.to_csv(f"{out}/niche_pseudobulk_lfc_per_section.csv.gz", index=False)
def meta(x):
    v = x.lfc.values
    return pd.Series({"n_sections": len(v), "mean_lfc": v.mean(), "frac_up": (v > 0.25).mean(),
                      "p": ttest_1samp(v, 0).pvalue if len(v) >= 3 else np.nan, "mean_expr_in": x.expr_in.mean()})
M = D.groupby(["niche", "gene"]).apply(meta).reset_index()
M = M[(M.n_sections >= 3) & (M.mean_expr_in > 0.3)].copy()
M["fdr"] = M.groupby("niche").p.transform(lambda p: multipletests(p.fillna(1), method="fdr_bh")[1])
top200 = M.sort_values("mean_lfc", ascending=False).groupby("niche").head(200)
promisc = top200.gene.value_counts(); promisc = set(promisc[promisc > 2].index)
cand = M[(M.fdr < 0.05) & (M.frac_up >= 0.75) & (M.mean_lfc > 0.3) & ~M.gene.isin(promisc)]
sig = cand.sort_values("mean_lfc", ascending=False).groupby("niche").head(50)
M.to_csv(f"{out}/niche_marker_meta.csv", index=False); sig.to_csv(f"{out}/niche_signatures.csv", index=False)
print(sig.groupby("niche").size())
for n, x in sig.groupby("niche"): print(n, ", ".join(x.gene.head(15)))
