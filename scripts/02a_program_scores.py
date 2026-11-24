"""
Step 2a - Spot-level program scores (no CNV needed)
Gene sets: Gavish et al. 2023 3CA meta-programs (glioma MES/AC/NPC/OPC, hypoxia, cycling, IFN)
plus canonical microenvironment / normal-brain markers (see gene_sets_used.json).
Usage: python 02a_program_scores.py <h5ad_dir> <gene_sets.json> <out_dir>
"""
import sys, os, glob, json, warnings
import numpy as np, pandas as pd, scanpy as sc
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")
h5_dir, gs_json, out_dir = sys.argv[1:4]
os.makedirs(f"{out_dir}/spots", exist_ok=True)
gene_sets = json.load(open(gs_json))
show = ["MES_like", "AC_like", "NPC_like", "OPC_like", "Hypoxia", "Microglia", "Mono_Macrophage", "Endothelial", "T_cell", "Oligodendrocyte", "Neuron"]

rows, objs = [], {}
for f in sorted(glob.glob(f"{h5_dir}/*.h5ad")):
    s = os.path.basename(f)[:-5]
    ad = sc.read_h5ad(f); sc.pp.normalize_total(ad, target_sum=1e4); sc.pp.log1p(ad)
    for k, gs in gene_sets.items():
        gs = [g for g in gs if g in ad.var_names]
        if len(gs) >= 3: sc.tl.score_genes(ad, gs, score_name=k, random_state=0)
        else: ad.obs[k] = np.nan  # too few marker genes detected in this section
    obs = ad.obs[[c for c in ["array_row", "array_col", "pxl_row", "pxl_col", "total_counts"] + list(gene_sets) if c in ad.obs]].copy()
    obs.insert(0, "sample", s); obs.to_csv(f"{out_dir}/spots/{s}_program_scores.csv.gz")
    if "spatial" in ad.obsm and ad.obsm["spatial"].any(): objs[s] = (ad.obsm["spatial"], obs)
    rows.append(dict(sample=s, **{k: round(float(obs[k].mean()), 4) for k in gene_sets if k in obs}))
    print(s, flush=True)
pd.DataFrame(rows).to_csv(f"{out_dir}/program_means_per_sample.csv", index=False)

# Figure: one row per sample, one column per program (spatial maps)
samples = list(objs)
fig, axes = plt.subplots(len(samples), len(show), figsize=(1.9 * len(show), 1.9 * len(samples)))
for i, s in enumerate(samples):
    xy, obs = objs[s]
    for j, k in enumerate(show):
        ax = axes[i, j]; v = obs[k].values
        lo, hi = np.percentile(v, [2, 98])
        ax.scatter(xy[:, 0], -xy[:, 1], c=np.clip(v, lo, hi), s=0.6, cmap="viridis", rasterized=True)
        ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect("equal")
        for sp in ax.spines.values(): sp.set_visible(False)
        if i == 0: ax.set_title(k.replace("_", " "), fontsize=9)
        if j == 0: ax.set_ylabel(s.replace("GBM_", "").replace("IDHm_", "IDHm "), fontsize=8)
fig.tight_layout()
fig.savefig(f"{out_dir}/program_maps_all_samples.png", dpi=90)
