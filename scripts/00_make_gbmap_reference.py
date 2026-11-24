"""
Make a small deconvolution reference from Core GBmap (CELLxGENE, 338,564 cells, 8.1 GB h5ad)
Runs on your PC; reads the big file in 'backed' mode, so it needs little RAM (~4-8 GB is enough).

Output: GBmap_reference_subsample.h5ad (~100-200 MB) with RAW COUNTS, up to 400 cells per cell type,
plus all annotation columns. Send/leave this small file in 02_metadata so Claude can upload it.

Usage (in Anaconda Prompt):
  python 00_make_gbmap_reference.py "D:\\GLIOMA SPATIAL PROJECT\\01_raw_data\\GBmap\\861acfd8-25f0-418b-a445-aa96da232827.h5ad"
"""
import sys, os
import numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp

src = sys.argv[1]
out = os.path.join(r"D:\GLIOMA SPATIAL PROJECT\02_metadata", "GBmap_reference_subsample.h5ad")
per_type = 400
rng = np.random.default_rng(0)

A = ad.read_h5ad(src, backed="r")
print("cells x genes:", A.shape)
print("annotation columns:", [c for c in A.obs.columns if "annot" in c.lower() or "cell_type" in c.lower()])

# choose the finest stable annotation available
for col in ["annotation_level_3", "annotation_level_2", "cell_type"]:
    if col in A.obs.columns:
        label = col; break
print("subsampling by:", label)
print(A.obs[label].value_counts().to_string())

idx = []
for ct, rows in A.obs.groupby(label, observed=True).indices.items():
    take = rows if len(rows) <= per_type else rng.choice(rows, per_type, replace=False)
    idx.extend(sorted(take))
idx = np.array(sorted(idx))
print("selected cells:", len(idx))

# raw counts: CELLxGENE stores them in .raw.X (X is normalised)
if A.raw is not None:
    X = A.raw.X[idx]
    var = A.raw.var.copy()
else:
    X = A.X[idx]; var = A.var.copy()
X = sp.csr_matrix(X)
obs = A.obs.iloc[idx].copy()
# use gene symbols as names when available (CELLxGENE: var['feature_name'])
if "feature_name" in var.columns:
    var.index = var["feature_name"].astype(str).values
var = var[[]]
B = ad.AnnData(X=X, obs=obs, var=var)
B.var_names_make_unique()
B.write_h5ad(out, compression="gzip")
print("written:", out, round(os.path.getsize(out) / 1e6, 1), "MB")
