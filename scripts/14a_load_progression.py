"""14a_load_progression.py -- load GSE318238 Visium CytAssist (patient BT2015050; 2015/2017/2019) into the atlas h5ad format.
QC as atlas: in_tissue spots, >=1 count, >=200 genes; MT genes flagged (probe set may lack them)."""
import sys, os, pandas as pd, numpy as np, scanpy as sc
src, out = sys.argv[1:3]; os.makedirs(out, exist_ok=True)
S = {"2015": "GSM9490659", "2017": "GSM9490660", "2019": "GSM9490661"}
rows = []
for yr, gsm in S.items():
    h5 = f"{src}/{gsm}_filtered_feature_bc_matrix.h5" if yr == "2015" else f"{src}/{gsm}_{yr}_filtered_feature_bc_matrix.h5"
    a = sc.read_10x_h5(h5); a.var_names_make_unique()
    p = pd.read_csv(f"{src}/{gsm}_{yr}_spatial_tissue_positions.csv.gz")
    p.columns = ["barcode", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col"][:p.shape[1]] if p.columns[0] != "barcode" else p.columns
    p = p.set_index("barcode")
    a = a[a.obs_names.isin(p.index[p.in_tissue == 1])].copy()
    a.obs = a.obs.join(p.loc[a.obs_names, ["in_tissue", "array_row", "array_col", "pxl_row_in_fullres" if "pxl_row_in_fullres" in p else "pxl_row", "pxl_col_in_fullres" if "pxl_col_in_fullres" in p else "pxl_col"]])
    a.obs.columns = ["in_tissue", "array_row", "array_col", "pxl_row", "pxl_col"]
    a.obsm["spatial"] = a.obs[["pxl_col", "pxl_row"]].values.astype(float)
    a.var["mt"] = a.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(a, qc_vars=["mt"], inplace=True)
    n0 = a.n_obs
    sc.pp.filter_cells(a, min_counts=1); sc.pp.filter_cells(a, min_genes=200)
    s = f"IDHm_GSE318238_{yr}"; a.obs["sample"] = s
    a.write_h5ad(f"{out}/{s}.h5ad")
    rows.append({"sample": s, "gsm": gsm, "year": yr, "spots_in_tissue": n0, "spots_after_QC": a.n_obs, "genes": a.n_vars,
                 "median_UMI": float(np.median(a.obs.total_counts)), "median_genes": float(np.median(a.obs.n_genes_by_counts)),
                 "pct_mt_median": float(np.median(a.obs.pct_counts_mt))})
pd.DataFrame(rows).to_csv(f"{out}/qc_summary.csv", index=False); print(pd.DataFrame(rows).to_string())
