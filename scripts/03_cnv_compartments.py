"""
Step 3 - CNV inference against normal cortex and tissue compartments
(glioma adaptation of the Li et al. 2026 Benign/Stromal/Malignant pipeline)

Reference: Ravi 2022 normal-appearing cortex sections UKF256_C and UKF265_C (the same references
used by Greenwald et al. 2024, Module5_CNA.R).
Gene positions: hg38 table from github.com/jlaffy/infercna.

Per tumor section:
  1. sample spots + reference spots -> log-normalised -> infercnvpy (window 100, step 10)
  2. spot CNV score = mean squared CNV signal; smoothed per Leiden expression cluster
     (the paper denoises CNV per expression cluster)
  3. Malignant = cluster CNV score > Otsu threshold (computed on sample + reference spots)
  4. Non-malignant spots split into TME-rich vs Brain by microenvironment vs normal-brain program score
  5. Arm-level CNV per spot: chr7, chr10, 1p, 19q (GBM hallmark +7/-10; oligodendroglioma -1p/-19q)

Usage: python 03_cnv_compartments.py <qc_h5ad_dirs(comma-sep)> <ref_h5ads(comma-sep)> <genes_hg38.csv> <gene_sets.json> <out_dir> [sample filter substr]
"""
import sys, os, glob, json, warnings
import numpy as np, pandas as pd, scanpy as sc, anndata as an, infercnvpy as cnv
from skimage.filters import threshold_otsu
warnings.filterwarnings("ignore")
h5_dirs, refs, genes_csv, gs_json, out_dir = sys.argv[1:6]
only = sys.argv[6] if len(sys.argv) > 6 else ""
os.makedirs(f"{out_dir}/spots", exist_ok=True)
gene_sets = json.load(open(gs_json))
genes = pd.read_csv(genes_csv).drop_duplicates("symbol").set_index("symbol")
CEN = {"1": 123_400_000, "19": 26_200_000}  # hg38 centromere midpoints (approx.)

ref = an.concat([sc.read_h5ad(r) for r in refs.split(",")], join="inner", label="ref_sample",
                keys=[os.path.basename(r)[:-5] for r in refs.split(",")], index_unique="-")
ref.obs["src"] = "ref"

def zs(x): return (x - x.mean()) / (x.std() + 1e-9)

files = sorted(f for d in h5_dirs.split(",") for f in glob.glob(f"{d}/*.h5ad"))
files = [f for f in files if not os.path.basename(f).endswith("_C.h5ad") and only in f]
summary = []
for f in files:
    s = os.path.basename(f)[:-5]
    if os.path.exists(f"{out_dir}/spots/{s}_cnv.csv.gz"): continue
    ad = sc.read_h5ad(f); ad.obs["src"] = "sample"
    spatial = ad.obsm.get("spatial")
    both = an.concat([ad, ref], join="inner", index_unique=None)
    both.obs_names_make_unique()
    sc.pp.normalize_total(both, target_sum=1e4); sc.pp.log1p(both)
    v = both.var_names.intersection(genes.index)
    both = both[:, v].copy()
    both.var["chromosome"] = "chr" + genes.loc[v, "chromosome_name"].astype(str).values
    both.var["start"] = genes.loc[v, "start_position"].values.astype(int)
    both.var["end"] = genes.loc[v, "end_position"].values.astype(int)
    both = both[:, both.var.chromosome.isin([f"chr{i}" for i in range(1, 23)])].copy()
    cnv.tl.infercnv(both, reference_key="src", reference_cat=["ref"], window_size=100, step=10, n_jobs=1)
    X = both.obsm["X_cnv"]; X = X.toarray() if hasattr(X, "toarray") else X
    both.obs["cnv_score"] = (X ** 2).mean(axis=1)

    # arm / chromosome level signal: map CNV windows back to gene start positions
    chr_pos = both.uns["cnv"]["chr_pos"]; chroms = list(chr_pos); starts = list(chr_pos.values()) + [X.shape[1]]
    var_sorted = both.var.sort_values(["chromosome", "start"])
    seg = {}
    for i, c in enumerate(chroms):
        cols = np.arange(starts[i], starts[i + 1])
        g = var_sorted[var_sorted.chromosome == c].start.values
        pos = np.array([g[min(len(g) - 1, k * 10 + 50)] for k in range(len(cols))])
        cn = c.replace("chr", "")
        if cn in CEN:
            p, q = cols[pos < CEN[cn]], cols[pos >= CEN[cn]]
            if cn == "1": seg["1p"] = X[:, p].mean(1)
            if cn == "19": seg["19q"] = X[:, q].mean(1)
        if cn in ("7", "10"): seg[c] = X[:, cols].mean(1)

    # subtype-specific malignancy signature (hallmark arm-level CNVs)
    is_ref = both.obs.src.values == "ref"
    if "IDHm" in s and s.endswith("O"):
        sig_all, sig_name = -(seg["1p"] + seg["19q"]), "-(1p+19q)"          # oligodendroglioma codeletion
    elif "IDHm" in s:
        sig_all, sig_name = both.obs["cnv_score"].values, "CNV magnitude"     # astrocytoma: no single hallmark arm
    else:
        sig_all, sig_name = seg["chr7"] - seg["chr10"], "chr7-chr10"         # GBM hallmark +7/-10
    th = float(np.percentile(sig_all[is_ref], 95))                            # cut-off from normal cortex spots
    smp = both[~is_ref].copy()
    for k, vv in seg.items(): smp.obs[f"cnv_{k}"] = vv[~is_ref]
    smp.obs["malig_signature"] = sig_all[~is_ref]
    # expression clusters on the sample spots; signature smoothed per cluster (paper: per-cluster denoising)
    tmp = smp.copy(); sc.pp.highly_variable_genes(tmp, n_top_genes=2000)
    sc.pp.pca(tmp, n_comps=20, mask_var="highly_variable"); sc.pp.neighbors(tmp)
    sc.tl.leiden(tmp, resolution=1.0, key_added="cl", flavor="igraph", n_iterations=2)
    smp.obs["cluster"] = tmp.obs["cl"].values
    smp.obs["cnv_cluster"] = smp.obs.groupby("cluster", observed=True)["malig_signature"].transform("mean")
    malig = smp.obs["cnv_cluster"].values > th

    # programme scores for non-malignant split
    for k in ["Microglia", "Mono_Macrophage", "Endothelial", "T_cell", "Oligodendrocyte", "Neuron"]:
        gs = [g for g in gene_sets[k] if g in smp.var_names]
        if len(gs) >= 3: sc.tl.score_genes(smp, gs, score_name=k, random_state=0)
        else: smp.obs[k] = 0.0
    tme = np.nanmean(np.vstack([zs(smp.obs[k]) for k in ["Microglia", "Mono_Macrophage", "Endothelial", "T_cell"]]), 0)
    brain = np.nanmean(np.vstack([zs(smp.obs[k]) for k in ["Oligodendrocyte", "Neuron"]]), 0)
    smp.obs["compartment"] = np.where(malig, "Malignant", np.where(tme > brain, "TME-rich", "Brain"))

    o = smp.obs[["cnv_score", "malig_signature", "cnv_cluster", "compartment"] + [c for c in smp.obs if c.startswith("cnv_") and c not in ("cnv_score", "cnv_cluster", "malig_signature")]].copy()
    o.insert(0, "sample", s)
    if spatial is not None and len(spatial) == len(o):
        o["x"], o["y"] = spatial[:, 0], spatial[:, 1]
    o.to_csv(f"{out_dir}/spots/{s}_cnv.csv.gz")
    m = o.compartment == "Malignant"
    row = dict(sample=s, signature=sig_name, ref95_threshold=round(th, 4), spots=len(o), pct_malignant=round(100 * m.mean(), 1),
               pct_TME=round(100 * (o.compartment == "TME-rich").mean(), 1),
               pct_brain=round(100 * (o.compartment == "Brain").mean(), 1))
    for k in seg:
        if m.any() and (~m).any():
            row[f"{k}_malig_vs_nonmalig"] = round(float(o.loc[m, f"cnv_{k}"].mean() - o.loc[~m, f"cnv_{k}"].mean()), 4)
        row[f"{k}_malig_mean"] = round(float(o.loc[m, f"cnv_{k}"].mean()), 4) if m.any() else np.nan
    summary.append(row); print(row, flush=True)
    pd.DataFrame([row]).to_csv(f"{out_dir}/summary_parts_{s}.csv", index=False)
    del both, smp, tmp
