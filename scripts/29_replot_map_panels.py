"""29_replot_map_panels.py -- re-draw the spatial-map and niche-profile panels from saved tables, as PNG plus SVG with live text.
Panels: compartment maps (Fig 1C), niche maps (Fig 2C), niche profiles from program scores (Fig S2A, S2B),
program maps per cohort (Fig S1A, S1B). Plotting code follows 03/04/06/02a; spots are rasterized in the SVG, text stays editable.
Usage: python3 29_replot_map_panels.py <outdir>
"""
import sys, os, glob, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
OUT = sys.argv[1] if len(sys.argv) > 1 else "edit_figs/panels"; os.makedirs(OUT, exist_ok=True)

def save(fig, name, dpi):
    fig.savefig(f"{OUT}/{name}.png", dpi=dpi)
    with matplotlib.rc_context({"svg.fonttype": "none"}):
        fig.savefig(f"{OUT}/{name}.svg", dpi=dpi)
    plt.close(fig)

def short(s): return s.replace("GBM_", "").replace("Ravi_", "R_").replace("IDHm_", "IDHm ")

def map_grid(D, colour_of, legend, name):
    ss = sorted(D["sample"].unique()); nc = 7; nr = int(np.ceil(len(ss) / nc))
    fig, axs = plt.subplots(nr, nc, figsize=(2.3 * nc, 2.5 * nr))
    for a, s in zip(axs.flat, ss):
        g = D[D["sample"] == s]; a.scatter(g.x, -g.y, c=colour_of(g), s=0.8, rasterized=True)
        a.set_title(short(s), fontsize=8); a.set_aspect("equal"); a.axis("off")
    for a in list(axs.flat)[len(ss):]: a.axis("off")
    fig.legend(handles=[plt.Line2D([], [], marker="o", ls="", color=c, label=k) for k, c in legend.items()], loc="lower right", fontsize=9, frameon=False)
    fig.tight_layout(); save(fig, name, 95)

# Fig 1C: copy-number compartments
C = pd.concat([pd.read_csv(f) for f in sorted(glob.glob("cnv_out/spots/*_cnv.csv.gz"))])
C = C[C["sample"].isin(pd.read_csv("niche_deconv/spot_LCP_niche.csv.gz", usecols=["sample"])["sample"].unique())]  # the 36 niche-analysis sections, as in the original panel
ccol = {"Malignant": "#c0392b", "TME-rich": "#2e86c1", "Brain": "#a9cce3"}
map_grid(C, lambda g: g.compartment.map(ccol).fillna("#dddddd"), ccol, "compartment_maps")

# Fig 2C: niches from deconvolution
N = pd.read_csv("niche_deconv/spot_LCP_niche.csv.gz"); K = N.niche.nunique()
cm = plt.get_cmap("tab10"); ncol = {f"Niche_{i+1}": cm(i) for i in range(K)}
map_grid(N, lambda g: g.niche.map(ncol), ncol, "niche_maps")

# Fig S2: niche profiles from program scores (with / without Harmony)
for src, name in [("niche_out", "niche_profiles_programscore"), ("niche_out_h", "niche_profiles_programscore_harmony")]:
    prog = pd.read_csv(f"{src}/niche_program_profile.csv", index_col=0); sub = pd.read_csv(f"{src}/niche_frequency_by_subtype.csv", index_col=0)
    order = prog.index.tolist(); FEAT = prog.columns.tolist()
    fig, ax = plt.subplots(1, 2, figsize=(16, 5), gridspec_kw={"width_ratios": [2, 1.1]})
    im = ax[0].imshow(prog.loc[order].values, cmap="RdBu_r", vmin=-1.5, vmax=1.5, aspect="auto")
    ax[0].set_xticks(range(len(FEAT))); ax[0].set_xticklabels([f.replace("_", " ") for f in FEAT], rotation=60, ha="right", fontsize=8)
    ax[0].set_yticks(range(len(order))); ax[0].set_yticklabels(order, fontsize=9); ax[0].set_title("Niche programme profile (mean z-score)")
    fig.colorbar(im, ax=ax[0], shrink=0.6)
    sub[order].T.plot.barh(ax=ax[1], color=["#c0392b", "#2e86c1"]); ax[1].invert_yaxis()
    ax[1].set_xlabel("Fraction of spots"); ax[1].set_ylabel("niche"); ax[1].set_title("Niche frequency: GBM vs IDH-mutant")
    fig.tight_layout(); save(fig, name, 110)

# Fig S1: program maps per cohort
show = ["MES_like", "AC_like", "NPC_like", "OPC_like", "Hypoxia", "Microglia", "Mono_Macrophage", "Endothelial", "T_cell", "Oligodendrocyte", "Neuron"]
for src, name in [("step2a_out", "program_maps_greenwald"), ("step2a_ravi", "program_maps_ravi")]:
    objs = []
    for f in sorted(glob.glob(f"{src}/spots/*_program_scores.csv.gz")):
        o = pd.read_csv(f, index_col=0)
        if {"pxl_row", "pxl_col"} <= set(o.columns) and o[["pxl_row", "pxl_col"]].notna().all().all() and o[["pxl_row", "pxl_col"]].abs().sum().sum() > 0:
            objs.append((o["sample"].iloc[0], o))
    fig, axes = plt.subplots(len(objs), len(show), figsize=(1.9 * len(show), 1.9 * len(objs)), squeeze=False)
    for i, (s, o) in enumerate(objs):
        for j, k in enumerate(show):
            ax = axes[i, j]
            if k in o and o[k].notna().any():
                v = o[k].values; lo, hi = np.nanpercentile(v, [2, 98])
                ax.scatter(o.pxl_col, -o.pxl_row, c=np.clip(v, lo, hi), s=0.6, cmap="viridis", rasterized=True)
            ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect("equal")
            for sp in ax.spines.values(): sp.set_visible(False)
            if i == 0: ax.set_title(k.replace("_", " "), fontsize=9)
            if j == 0: ax.set_ylabel(s.replace("GBM_", "").replace("IDHm_", "IDHm "), fontsize=8)
    fig.tight_layout(); save(fig, name, 90)
print("ok", sorted(os.listdir(OUT)))
