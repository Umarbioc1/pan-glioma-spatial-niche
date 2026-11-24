"""
Step 5 - Spot deconvolution with the GBmap reference (RCTD-like, Python)

Reference: Core GBmap (Ruiz-Moreno et al., Neuro-Oncology 2025; CELLxGENE), subsampled to <=400 cells per
annotation_level_3 type (00_make_gbmap_reference.py), raw counts.

Model (same principle as RCTD's 'full' mode, Cable et al. 2022):
  spot counts y_g  ~  N * sum_k w_k * gamma_g * mu_gk
  mu_gk   = mean normalised expression of gene g in cell type k (reference)
  gamma_g = gene-wise platform effect (scRNA -> Visium), estimated from all spots of a cohort (2 iterations)
  w       = non-negative cell-type weights, solved per spot by NNLS on variance-stabilised data, then scaled to sum 1
Genes: union of the top 40 marker genes per cell type (log fold-change vs all other types, min mean expression).

Usage: python 05_deconvolution.py <ref.h5ad> <h5ad_dirs comma-separated> <out_dir>
"""
import sys, os, glob, numpy as np, pandas as pd, scanpy as sc, anndata as an, warnings
from scipy.optimize import nnls
warnings.filterwarnings("ignore")
ref_f, dirs, out = sys.argv[1:4]
os.makedirs(f"{out}/spots", exist_ok=True)

# ---------- reference profiles ----------
R = an.read_h5ad(ref_f); R.var_names_make_unique()
R = R[:, ~R.var_names.str.match(r"^(MT-|RPL|RPS|ENSG)")].copy()
sc.pp.normalize_total(R, target_sum=1e4)
types = sorted(R.obs["annotation_level_3"].unique())
mu = pd.DataFrame({t: np.asarray(R[R.obs.annotation_level_3 == t].X.mean(0)).ravel() for t in types}, index=R.var_names)

# marker genes
markers = set()
lmu = np.log1p(mu)
for t in types:
    other = lmu.drop(columns=t).max(axis=1)
    lfc = lmu[t] - other
    ok = mu[t] > 0.5  # at least 0.5 counts per 10k on average
    markers |= set(lfc[ok].sort_values(ascending=False).index[:40])
markers = sorted(markers)
print(f"{len(types)} cell types, {len(markers)} marker genes")
mu.loc[markers].to_csv(f"{out}/reference_profiles_markers.csv")

def deconv(Y, M):
    """Y: spots x genes (counts), M: genes x types (normalised profiles). Returns weights spots x types."""
    N = Y.sum(1, keepdims=True); N[N == 0] = 1
    Yn = Y / N * 1e4                                 # same scale as profiles
    gamma = np.ones(M.shape[0])
    for it in range(3):
        Mg = M * gamma[:, None]
        sd = np.sqrt(Mg.mean(1) + 1e-3)               # variance stabilisation weights
        A = Mg / sd[:, None]
        W = np.vstack([nnls(A, Yn[i] / sd)[0] for i in range(Yn.shape[0])])
        pred = W @ Mg.T
        if it < 2:  # update platform effect: observed / predicted mean per gene (bounded)
            gamma = gamma * np.clip((Yn.mean(0) + 1e-3) / (pred.mean(0) + 1e-3), 0.1, 10)
    P = W / np.maximum(W.sum(1, keepdims=True), 1e-9)
    return P, gamma

summ = []
files = sorted(f for d in dirs.split(",") for f in glob.glob(f"{d}/*.h5ad"))
for f in files:
    s = os.path.basename(f)[:-5]
    a = sc.read_h5ad(f); a.var_names_make_unique()
    genes = [m for m in markers if m in a.var_names]           # per-section marker set (no cross-section intersection)
    M = mu.loc[genes].values
    Y = np.asarray(a[:, genes].X.todense())
    P, gamma = deconv(Y, M)                                     # platform effect estimated within the section
    Mg = M * gamma[:, None]; N = Y.sum(1, keepdims=True); N[N == 0] = 1; Yn = Y / N * 1e4
    pred = (P * 1.0) @ Mg.T
    r = np.array([np.corrcoef(Yn[i], pred[i])[0, 1] if Yn[i].std() > 0 else np.nan for i in range(len(Yn))])
    df = pd.DataFrame(P, index=a.obs_names, columns=types); df["fit_r"] = r
    df.insert(0, "sample", s); df.to_csv(f"{out}/spots/{s}_proportions.csv.gz")
    summ.append(dict(sample=s, genes=len(genes), spots=len(df), median_fit_r=round(float(np.nanmedian(r)), 3),
                     **{t: round(float(df[t].mean()), 4) for t in types}))
    print(s, len(genes), round(float(np.nanmedian(r)), 3), flush=True)
pd.DataFrame(summ).to_csv(f"{out}/mean_proportions_per_sample.csv", index=False)
