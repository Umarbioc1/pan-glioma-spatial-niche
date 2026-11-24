"""
Step 4 - Local cellular programs (LCPs) and spatial niches  [INTERIM: program scores as composition proxy]

Follows Li et al. 2026 (Cell Rep Med) Methods:
  * LCPs: k-means on per-spot cell composition, separately in Malignant and non-malignant compartments
    (MLCP / TLCP). Composition proxy here = program scores (Gavish 2023 MPs + markers), z-scored within
    cohort, clipped at 3. To be replaced by RCTD cell-type proportions once the GBmap reference is available.
  * Niches: LCP composition of neighbouring spots within two rings (Visium hex grid), z-scored, clipped at 3,
    k-means; niche = recurrent multicellular neighbourhood.
Sections with <5% malignant spots (no confident CNV call) contribute non-malignant spots only.

Usage: python 04_lcp_niches.py <out_dir>
"""
import sys, os, glob, numpy as np, pandas as pd
from sklearn.cluster import KMeans
from scipy.spatial import cKDTree
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
out = sys.argv[1]; os.makedirs(out, exist_ok=True)
rng = 0
MAL = ["MES_like", "AC_like", "NPC_like", "OPC_like", "Hypoxia", "Cycling", "Interferon_MHC"]
TME = ["Microglia", "Mono_Macrophage", "Neutrophil", "T_cell", "CD8_T_NK", "B_cell", "Plasma", "DC",
       "Endothelial", "Mural_pericyte", "Fibroblast"]
BRAIN = ["Oligodendrocyte", "Neuron", "Astrocyte_normal"]
FEAT = MAL + TME + BRAIN

# ---------- assemble spots ----------
cs = pd.read_csv("cnv_out/compartment_summary.csv").set_index("sample")
rows = []
for f in sorted(glob.glob("cnv_out/spots/*_cnv.csv.gz")):
    c = pd.read_csv(f, index_col=0); s = c["sample"].iloc[0]
    pdir = "step2a_ravi" if s.startswith("Ravi_") else "step2a_out"
    p = pd.read_csv(f"{pdir}/spots/{s}_program_scores.csv.gz", index_col=0)
    d = c.join(p.drop(columns=["sample"]), how="inner")
    if "array_row" not in d or d["array_row"].isna().all() or "x" not in d: continue  # no spatial coordinates
    d["cohort"] = "Ravi2022" if s.startswith("Ravi_") else "Greenwald2024"
    d["subtype"] = "IDH-mut" if "IDHm" in s else "GBM"
    d["confident_cnv"] = cs.loc[s, "pct_malignant"] >= 5
    d["spot"] = d.index; rows.append(d)
D = pd.concat(rows, ignore_index=True)
D[FEAT] = D[FEAT].fillna(0)
Z = D.groupby("cohort")[FEAT].transform(lambda x: (x - x.mean()) / (x.std() + 1e-9)).clip(-3, 3)
# batch correction across cohorts (Harmony on the programme matrix; corrected values used for clustering)
import harmonypy
ho = harmonypy.run_harmony(Z[FEAT].values, D[["cohort"]], ["cohort"], max_iter_harmony=20, random_state=0)
Zc = ho.Z_corr.T if ho.Z_corr.shape[0] == len(FEAT) else ho.Z_corr
Z = pd.DataFrame(np.asarray(Zc), columns=FEAT, index=D.index).clip(-3, 3)
D["compartment2"] = np.where((D.compartment == "Malignant") & D.confident_cnv, "Malignant", "NonMalignant")
D.loc[(D.compartment == "Malignant") & ~D.confident_cnv, "compartment2"] = "NonMalignant"
print(D.groupby(["subtype", "compartment2"]).size())

# ---------- LCPs ----------
def name_cluster(center, cols, prefix):
    top = pd.Series(center, index=cols).sort_values(ascending=False)
    keep = [k for k in top.index[:2] if top[k] > 0.5] or [top.index[0]]
    return prefix + "_" + "/".join(k.replace("_like", "").replace("Mono_Macrophage", "Macro").replace("Oligodendrocyte", "Oligo")
                                    .replace("Astrocyte_normal", "Astro").replace("Mural_pericyte", "Mural").replace("Interferon_MHC", "IFN")
                                    .replace("Endothelial", "Endo").replace("CD8_T_NK", "CD8T/NK") for k in keep)

D["LCP"] = ""
for comp, prefix, k in [("Malignant", "MLCP", 14), ("NonMalignant", "TLCP", 14)]:
    idx = D.compartment2 == comp
    km = KMeans(n_clusters=k, n_init=20, random_state=rng).fit(Z.loc[idx, FEAT].values)
    names = [name_cluster(c, FEAT, prefix) for c in km.cluster_centers_]
    # make names unique
    seen = {}; uniq = []
    for n in names:
        seen[n] = seen.get(n, 0) + 1; uniq.append(n if seen[n] == 1 else f"{n}#{seen[n]}")
    D.loc[idx, "LCP"] = np.array(uniq)[km.labels_]
    pd.DataFrame(km.cluster_centers_, index=uniq, columns=FEAT).round(3).to_csv(f"{out}/{prefix}_centers.csv")
lcps = sorted(D.LCP.unique())

# ---------- niches: neighbourhood LCP composition within 2 rings ----------
NB = np.zeros((len(D), len(lcps)), dtype=np.float32); li = {l: i for i, l in enumerate(lcps)}
lab = D.LCP.map(li).values
for s, g in D.groupby("sample"):
    xy = np.c_[g.array_col.values, g.array_row.values * np.sqrt(3)]  # hex grid -> nearest neighbour distance 2
    tree = cKDTree(xy); nbrs = tree.query_ball_point(xy, r=4.1)        # 2 rings = up to 18 neighbours + self
    gi = g.index.values
    for j, nb in enumerate(nbrs):
        cnt = np.bincount(lab[gi[nb]], minlength=len(lcps)); NB[gi[j]] = cnt / cnt.sum()
NBz = pd.DataFrame(NB, columns=lcps)
NBz = NBz.groupby(D.subtype.values).transform(lambda x: (x - x.mean()) / (x.std() + 1e-9)).clip(-3, 3)
K = 10
km = KMeans(n_clusters=K, n_init=20, random_state=rng).fit(NBz.values)
D["niche"] = [f"Niche_{i+1}" for i in km.labels_]
# characterise niches by mean neighbourhood LCP composition and programme z-scores
comp = pd.DataFrame(NB, columns=lcps).groupby(D.niche.values).mean()
prog = Z[FEAT].groupby(D.niche.values).mean()
comp.round(3).to_csv(f"{out}/niche_LCP_composition.csv"); prog.round(3).to_csv(f"{out}/niche_program_profile.csv")
freq = pd.crosstab(D["sample"], D.niche, normalize="index")
freq.round(4).to_csv(f"{out}/niche_frequency_per_sample.csv")
sub = pd.crosstab(D.subtype, D.niche, normalize="index"); sub.round(4).to_csv(f"{out}/niche_frequency_by_subtype.csv")
D[["sample", "spot", "cohort", "subtype", "compartment2", "LCP", "niche", "x", "y", "array_row", "array_col"]].to_csv(f"{out}/spot_LCP_niche.csv.gz", index=False)

# ---------- figures ----------
order = prog.index.tolist()
fig, ax = plt.subplots(1, 2, figsize=(16, 5), gridspec_kw={"width_ratios": [2, 1.1]})
im = ax[0].imshow(prog.loc[order].values, cmap="RdBu_r", vmin=-1.5, vmax=1.5, aspect="auto")
ax[0].set_xticks(range(len(FEAT))); ax[0].set_xticklabels([f.replace("_", " ") for f in FEAT], rotation=60, ha="right", fontsize=8)
ax[0].set_yticks(range(K)); ax[0].set_yticklabels(order, fontsize=9); ax[0].set_title("Niche programme profile (mean z-score)")
fig.colorbar(im, ax=ax[0], shrink=0.6)
sub[order].T.plot.barh(ax=ax[1], color=["#c0392b", "#2e86c1"]); ax[1].invert_yaxis()
ax[1].set_xlabel("Fraction of spots"); ax[1].set_title("Niche frequency: GBM vs IDH-mutant")
fig.tight_layout(); fig.savefig(f"{out}/niche_profiles.png", dpi=110)

cmap = plt.get_cmap("tab10"); col = {f"Niche_{i+1}": cmap(i) for i in range(K)}
ss = sorted(D["sample"].unique()); nc = 7; nr = int(np.ceil(len(ss) / nc))
fig, axs = plt.subplots(nr, nc, figsize=(2.3 * nc, 2.5 * nr))
for a, s in zip(axs.flat, ss):
    g = D[D["sample"] == s]; a.scatter(g.x, -g.y, c=g.niche.map(col), s=0.8, rasterized=True)
    a.set_title(s.replace("GBM_", "").replace("Ravi_", "R_").replace("IDHm_", "IDHm "), fontsize=8); a.set_aspect("equal"); a.axis("off")
for a in list(axs.flat)[len(ss):]: a.axis("off")
fig.legend(handles=[plt.Line2D([], [], marker="o", ls="", color=c, label=k) for k, c in col.items()], loc="lower right", fontsize=9, frameon=False)
fig.tight_layout(); fig.savefig(f"{out}/niche_maps.png", dpi=95)
print(prog.round(2).to_string()); print(sub.round(3).to_string())
