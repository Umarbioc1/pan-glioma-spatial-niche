"""
Step 2 - Spot-level cell-program scores and tissue compartments (Greenwald 2024 Visium)

Interim, reference-free version (GBmap/RCTD deconvolution to follow once the reference is downloaded).
 * Malignant-state and microenvironment programs: Gavish et al. 2023 (Nature) 3CA meta-programs
   (github.com/tiroshlab/3ca) + canonical cell-type markers.
 * CNV: infercnvpy, hg38 gene positions from github.com/jlaffy/infercna.
   Reference spots = within-sample spots richest in normal brain (oligodendrocyte + neuron score)
   and poorest in malignant-program score (the paper's "top 1% non-parenchymal" rule, adapted to brain).
 * Compartments (glioma adaptation of Li et al. 2026):
     Malignant  = CNV score above the sample's Otsu threshold
     TME-rich   = low CNV, microenvironment (myeloid+vascular+lymphoid) score > brain score
     Brain      = low CNV, normal-brain (oligo+neuron+astro) score >= microenvironment score

Usage: python 02_programs_and_compartments.py <h5ad_dir> <genes_hg38.csv> <gavish_MPs.json> <out_dir>
"""
import sys, os, glob, json
import numpy as np, pandas as pd, scanpy as sc, infercnvpy as cnv
from skimage.filters import threshold_otsu
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

h5_dir, genes_csv, mp_json, out_dir = sys.argv[1:5]
os.makedirs(out_dir, exist_ok=True); os.makedirs(f"{out_dir}/h5ad", exist_ok=True)

mps = json.load(open(mp_json))
G = lambda k: mps[k]
gene_sets = {
    # malignant states (Gavish 2023 cancer MPs; glioma-specific ones map to Neftel states)
    "MES_like": G("Cancer|MP16 MES (glioma)"),
    "AC_like": G("Cancer|MP25 Astrocytes"),
    "NPC_like": G("Cancer|MP26 NPC Glioma"),
    "OPC_like": G("Cancer|MP27 Oligo Progenitor"),
    "Hypoxia": G("Cancer|MP6 Hypoxia"),
    "Cycling": G("Cancer|MP1  Cell Cycle - G2/M") + G("Cancer|MP2  Cell Cycle - G1/S"),
    "Interferon_MHC": G("Cancer|MP17 Interferon/MHC-II (I)"),
    # microenvironment / normal brain (canonical markers)
    "Microglia": ["P2RY12", "TMEM119", "CX3CR1", "SALL1", "SIGLEC8", "CSF1R", "GPR34", "OLFML3", "P2RY13", "SELPLG"],
    "Mono_Macrophage": ["CD163", "F13A1", "LYZ", "FCN1", "VCAN", "TGFBI", "CD14", "MRC1", "LYVE1", "S100A11"],
    "Neutrophil": ["S100A8", "S100A9", "FCGR3B", "CSF3R", "CXCR2", "FPR1", "MNDA"],
    "T_cell": ["CD3D", "CD3E", "CD2", "TRAC", "CD247", "IL7R"],
    "CD8_T_NK": ["CD8A", "CD8B", "NKG7", "GZMK", "GZMA", "GNLY", "KLRD1", "PRF1"],
    "B_cell": ["MS4A1", "CD79A", "CD79B", "CD19", "BANK1"],
    "Plasma": ["JCHAIN", "MZB1", "IGKC", "IGHG1", "IGHA1", "XBP1"],
    "DC": ["CLEC9A", "XCR1", "CD1C", "FCER1A", "CLEC10A", "LAMP3"],
    "Endothelial": ["PECAM1", "CLDN5", "VWF", "FLT1", "CDH5", "ESAM", "KDR"],
    "Mural_pericyte": ["PDGFRB", "RGS5", "ACTA2", "KCNJ8", "MYH11", "NOTCH3"],
    "Fibroblast": ["COL1A1", "COL1A2", "DCN", "LUM", "COL3A1"],
    "Oligodendrocyte": ["MBP", "PLP1", "MOG", "MOBP", "MAG", "CNP", "ERMN", "OPALIN", "TF", "CLDN11"],
    "Neuron": ["SNAP25", "SYT1", "RBFOX3", "STMN2", "GAD1", "GAD2", "SLC17A7", "NRGN", "VSNL1", "CALM1"],
    "Astrocyte_normal": ["AQP4", "GJA1", "SLC1A2", "SLC1A3", "ALDOC", "GLUL", "SPARCL1", "MLC1"],
}
json.dump(gene_sets, open(f"{out_dir}/gene_sets_used.json", "w"), indent=1)

genes = pd.read_csv(genes_csv).drop_duplicates("symbol").set_index("symbol")
chrom_ok = [str(c) for c in list(range(1, 23))]  # autosomes only

summary = []
for f in sorted(glob.glob(f"{h5_dir}/*.h5ad")):
    s = os.path.basename(f)[:-5]
    ad = sc.read_h5ad(f)
    ad.layers["counts"] = ad.X.copy()
    sc.pp.normalize_total(ad, target_sum=1e4); sc.pp.log1p(ad)
    for k, gs in gene_sets.items():
        gs = [g for g in gs if g in ad.var_names]
        sc.tl.score_genes(ad, gs, score_name=k, random_state=0) if len(gs) >= 3 else ad.obs.__setitem__(k, np.nan)
    z = lambda c: (ad.obs[c] - ad.obs[c].mean()) / (ad.obs[c].std() + 1e-9)
    ad.obs["malig_prog"] = ad.obs[["MES_like", "AC_like", "NPC_like", "OPC_like"]].max(axis=1)
    ad.obs["brain_score"] = (z("Oligodendrocyte") + z("Neuron")) / 2
    ad.obs["tme_score"] = (z("Microglia") + z("Mono_Macrophage") + z("Endothelial") + z("T_cell")) / 4

    # --- CNV ---
    v = ad.var_names.intersection(genes.index)
    ad = ad[:, v].copy()
    ad.var["chromosome"] = "chr" + genes.loc[v, "chromosome_name"].astype(str).values
    ad.var["start"] = genes.loc[v, "start_position"].values.astype(int)
    ad.var["end"] = genes.loc[v, "end_position"].values.astype(int)
    ad = ad[:, ad.var.chromosome.isin(["chr" + c for c in chrom_ok])].copy()
    ref_rank = (ad.obs["brain_score"].rank(pct=True) - ad.obs["malig_prog"].rank(pct=True))
    n_ref = max(30, int(0.05 * ad.n_obs))
    ref_idx = ref_rank.sort_values(ascending=False).index[:n_ref]
    ad.obs["cnv_ref"] = np.where(ad.obs_names.isin(ref_idx), "ref", "query")
    cnv.tl.infercnv(ad, reference_key="cnv_ref", reference_cat=["ref"], window_size=100, step=10, n_jobs=1)
    X = ad.obsm["X_cnv"]; X = X.toarray() if hasattr(X, "toarray") else X
    ad.obs["cnv_score"] = np.sqrt((X ** 2).mean(axis=1))
    # smooth CNV score over Leiden clusters (paper denoises per expression cluster)
    sc.pp.highly_variable_genes(ad, n_top_genes=2000); sc.pp.pca(ad, n_comps=20, mask_var="highly_variable")
    sc.pp.neighbors(ad); sc.tl.leiden(ad, resolution=0.8, key_added="cluster", flavor="igraph", n_iterations=2)
    ad.obs["cnv_cluster"] = ad.obs.groupby("cluster", observed=True)["cnv_score"].transform("mean")
    th = threshold_otsu(ad.obs["cnv_cluster"].values)
    malig = ad.obs["cnv_cluster"] > th
    ad.obs["compartment"] = np.where(malig, "Malignant",
                              np.where(ad.obs["tme_score"] > ad.obs["brain_score"], "TME-rich", "Brain"))
    # chr7 gain / chr10 loss check (GBM hallmark) in malignant vs other spots
    chr_pos = ad.uns["cnv"]["chr_pos"]; chroms = list(chr_pos.keys()); starts = list(chr_pos.values()) + [X.shape[1]]
    seg = {c: X[:, starts[i]:starts[i + 1]].mean(axis=1) for i, c in enumerate(chroms)}
    for c in ["chr7", "chr10", "chr1", "chr19"]:
        if c in seg: ad.obs[f"cnv_{c}"] = seg[c]
    row = dict(sample=s, spots=ad.n_obs, otsu=float(th),
               pct_malignant=round(100 * malig.mean(), 1),
               pct_TME=round(100 * (ad.obs.compartment == "TME-rich").mean(), 1),
               pct_brain=round(100 * (ad.obs.compartment == "Brain").mean(), 1))
    for c in ["chr7", "chr10", "chr1", "chr19"]:
        if f"cnv_{c}" in ad.obs:
            row[f"{c}_malig_minus_other"] = round(float(ad.obs.loc[malig, f"cnv_{c}"].mean() - ad.obs.loc[~malig, f"cnv_{c}"].mean()), 4)
    summary.append(row); print(row, flush=True)
    keep_obs = ["sample", "array_row", "array_col", "pxl_row", "pxl_col", "total_counts", "cnv_score", "cnv_cluster",
                "compartment", "malig_prog", "brain_score", "tme_score"] + list(gene_sets) + [c for c in ad.obs if c.startswith("cnv_chr")]
    out = sc.AnnData(obs=ad.obs[[c for c in keep_obs if c in ad.obs]].copy())
    out.obsm["spatial"] = ad.obsm.get("spatial", np.zeros((ad.n_obs, 2)))
    out.write_h5ad(f"{out_dir}/h5ad/{s}_spots.h5ad", compression="gzip")
    del ad

pd.DataFrame(summary).to_csv(f"{out_dir}/compartment_summary.csv", index=False)
