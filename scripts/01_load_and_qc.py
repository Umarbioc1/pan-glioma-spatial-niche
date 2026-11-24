"""
Step 1 - Load Greenwald 2024 Visium samples and run sample-level QC
Pan-glioma spatial atlas (framework: Li et al., Cell Rep Med 2026)

QC thresholds from the paper:
  spots under tissue > 500, median UMI per spot > 500, median mitochondrial fraction < 0.2
Gene filter for downstream use: genes detected in >= 3 spots; mitochondrial genes removed.

Usage: python 01_load_and_qc.py <raw_dir> <out_dir>
  raw_dir: folder containing one unpacked folder per sample (e.g. GBM_MGH258/)
"""
import sys, os, glob, json
import numpy as np, pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

raw_dir, out_dir = sys.argv[1], sys.argv[2]
os.makedirs(out_dir, exist_ok=True)
os.makedirs(os.path.join(out_dir, "h5ad"), exist_ok=True)

rows, adatas = [], {}
for sdir in sorted(glob.glob(os.path.join(raw_dir, "*/"))):
    sample = os.path.basename(os.path.normpath(sdir))
    h5 = glob.glob(os.path.join(sdir, "*filtered_feature_bc_matrix.h5"))
    if not h5:
        print(f"[skip] {sample}: no h5"); continue
    ad = sc.read_10x_h5(h5[0]); ad.var_names_make_unique()
    pos_f = glob.glob(os.path.join(sdir, "*tissue_positions_list.csv"))
    has_pos = bool(pos_f)
    if has_pos:
        pos = pd.read_csv(pos_f[0], header=None,
                          names=["barcode", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col"])
        pos = pos.set_index("barcode")
        pos = pos[pos.in_tissue == 1]
        keep = ad.obs_names.intersection(pos.index)
        ad = ad[keep].copy()
        ad.obs = ad.obs.join(pos.loc[keep])
        ad.obsm["spatial"] = ad.obs[["pxl_col", "pxl_row"]].to_numpy().astype(float)
        sf = json.load(open(glob.glob(os.path.join(sdir, "*scalefactors_json.json"))[0]))
        img = np.array(Image.open(glob.glob(os.path.join(sdir, "*tissue_lowres_image.png"))[0]))
        ad.uns["spatial"] = {sample: {"images": {"lowres": img}, "scalefactors": sf}}
    n_empty = int((np.asarray(ad.X.sum(1)).ravel() == 0).sum())
    sc.pp.filter_cells(ad, min_counts=1)  # drop spots with zero UMIs (would give NaN mito fraction)
    ad.var["mt"] = ad.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(ad, qc_vars=["mt"], inplace=True, percent_top=None)
    ad.obs["sample"] = sample
    n_spots = ad.n_obs
    med_umi = float(np.median(ad.obs.total_counts))
    med_genes = float(np.median(ad.obs.n_genes_by_counts))
    has_mt = bool(ad.var.mt.any())
    med_mt = float(np.median(ad.obs.pct_counts_mt)) / 100 if has_mt else np.nan
    # if the authors removed MT- genes from the matrix, mito QC cannot be applied -> judged on the other two criteria
    passed = (n_spots > 500) and (med_umi > 500) and ((not has_mt) or med_mt < 0.2)
    rows.append(dict(sample=sample,
                     subtype=("IDH-mutant glioma" if "IDHm" in sample
                              else "Cortex (Ravi normal-appearing control)" if sample.endswith("_C")
                              else "GBM (IDH-wt)"),
                     spots_under_tissue=n_spots, median_UMI=round(med_umi), median_genes=round(med_genes),
                     median_mito_frac=round(med_mt, 4), genes_detected=int((ad.var.n_cells_by_counts > 0).sum()),
                     spatial_coords=has_pos, QC_pass=passed,
                     note=("; ".join(filter(None, [
                         "" if has_mt else "MT- genes absent from matrix (mito QC not applicable)",
                         f"{n_empty} zero-UMI spot(s) removed" if n_empty else "",
                         "" if has_pos else "tissue_positions_list.csv missing: no spatial coordinates, cannot be used for niche analysis"])))))
    adatas[sample] = ad
    print(f"{sample}: spots={n_spots} medUMI={med_umi:.0f} medMT={med_mt:.3f} pass={passed}")

qc = pd.DataFrame(rows)
qc.to_csv(os.path.join(out_dir, "QC_summary.csv"), index=False)

# Spatial QC panels: total UMI per spot on the H&E (lowres)
withpos = [s for s in adatas if "spatial" in adatas[s].uns]
n = len(withpos); ncol = 5; nrow = int(np.ceil(n / ncol))
fig, axes = plt.subplots(nrow, ncol, figsize=(4 * ncol, 4 * nrow))
for ax, s in zip(axes.flat, withpos):
    ad = adatas[s]; sf = ad.uns["spatial"][s]["scalefactors"]["tissue_lowres_scalef"]
    ax.imshow(ad.uns["spatial"][s]["images"]["lowres"])
    xy = ad.obsm["spatial"] * sf
    p = ax.scatter(xy[:, 0], xy[:, 1], c=np.log10(ad.obs.total_counts + 1), s=2, cmap="magma")
    r = qc.set_index("sample").loc[s]
    ax.set_title(f"{s}\n{r.spots_under_tissue} spots, med UMI {r.median_UMI}", fontsize=9)
    ax.axis("off")
for ax in list(axes.flat)[n:]:
    ax.axis("off")
fig.colorbar(p, ax=axes, shrink=0.4, label="log10 UMI per spot")
fig.savefig(os.path.join(out_dir, "QC_spatial_UMI.png"), dpi=110, bbox_inches="tight")

# Save per-sample objects (QC-passed), mito genes removed, genes in >= 3 spots kept
for s, ad in adatas.items():
    if qc.set_index("sample").loc[s, "QC_pass"]:
        ad = ad[:, ~ad.var.mt].copy()
        sc.pp.filter_genes(ad, min_cells=3)
        ad.uns.pop("spatial", None)  # keep files small; images stay in raw folders
        ad.write_h5ad(os.path.join(out_dir, "h5ad", f"{s}.h5ad"), compression="gzip")
print(qc.to_string(index=False))
