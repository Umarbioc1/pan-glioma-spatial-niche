"""17b - Figure S5: IvyGAP anatomical validation (style matches 16_article_figures.py)."""
import sys, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7, "axes.titlesize": 7.5, "axes.labelsize": 7,
                     "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6.5, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.linewidth": 0.6, "pdf.fonttype": 42})
IN = sys.argv[1] if len(sys.argv) > 1 else "ivy_out"; OUT = sys.argv[2] if len(sys.argv) > 2 else "manuscript_figs"
df = pd.read_csv(f"{IN}/ivygap_sample_scores.csv", index_col=0)
means = pd.read_csv(f"{IN}/ivygap_region_means.csv", index_col=0)
res = pd.read_csv(f"{IN}/ivygap_mixedlm_region_contrasts.csv")
order = ["LE", "IT", "CT", "CTmvp", "CTpan"]
lab = {"LE": "LE", "IT": "IT", "CT": "CT", "CTmvp": "MVP", "CTpan": "PAN"}
col = {"LE": "#8fb3d9", "IT": "#6a9fcf", "CT": "#9a9a9a", "CTmvp": "#d99a3d", "CTpan": "#c0392b"}

def L(ax, s, dx=-0.12, dy=1.04):
    ax.text(dx, dy, s, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom", ha="right")

fig = plt.figure(figsize=(7.2, 5.6))
gs = fig.add_gridspec(2, 1, height_ratios=[1, 1.05], hspace=0.75)
top = gs[0].subgridspec(1, 6, wspace=0.45)
panels = [("HypoxicMES_N6", "N6 hypoxic MES"), ("MP_Hypoxia", "Hypoxia MP"), ("Myeloid_N3", "N3 myeloid"),
          ("ANGPTL4", "ANGPTL4"), ("Neuronal_N1", "N1 neuronal"), ("Oligodendroglial_N10", "N10 oligodendroglial")]
rng = np.random.default_rng(0)
for i, (c, t) in enumerate(panels):
    ax = fig.add_subplot(top[i])
    data = [df.loc[df.region == r, c].values for r in order]
    bp = ax.boxplot(data, widths=0.6, patch_artist=True, showfliers=False, medianprops=dict(color="k", lw=0.8),
                    whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6), boxprops=dict(lw=0.6))
    for p, r in zip(bp["boxes"], order): p.set_facecolor(col[r]); p.set_alpha(0.55)
    for j, (r, d) in enumerate(zip(order, data)):
        ax.scatter(j + 1 + rng.uniform(-0.18, 0.18, len(d)), d, s=3, color=col[r], lw=0, alpha=0.9)
    ax.set_xticks(range(1, 6)); ax.set_xticklabels([lab[r] for r in order], rotation=90)
    rr = res[(res.score == c) & (res.contrast == "CTpan vs CT") & (res.test == "OLS+tumour FE")].iloc[0]
    ax.set_title(f"{t}\n" + f"PAN vs CT β {rr.beta:+.2f}\nFDR {rr.fdr:.0e}".replace("e-0", "e-"), fontsize=6.3)
    if i == 0: ax.set_ylabel("score (mean z)"); L(ax, "A", -0.35, 1.18)

genes = ["HypoxicMES_N6", "MP_Hypoxia", "ANGPTL4", "SDC4", "CD44", "SDC2", "ITGA5", "ITGB1",
         "CAV1", "VMP1", "NRN1", "NDRG1", "NAMPT", "ARRDC3", "PDPN", "CA12", "CAMK2N1"]
names = {"HypoxicMES_N6": "N6 signature", "MP_Hypoxia": "Hypoxia MP"}
M = means.loc[genes, order].astype(float)
ax = fig.add_subplot(gs[1])
im = ax.imshow(M.T.values, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
ax.set_xticks(range(len(genes))); ax.set_xticklabels([names.get(g, g) for g in genes], rotation=55, ha="right")
ax.set_yticks(range(5)); ax.set_yticklabels(["Leading edge (LE)", "Infiltrating (IT)", "Cellular tumour (CT)",
                                             "Microvascular prolif. (MVP)", "Pseudopalisading (PAN)"])
for j, g in enumerate(genes):
    i = order.index(M.loc[g].idxmax()); ax.add_patch(plt.Rectangle((j - .5, i - .5), 1, 1, fill=False, lw=1, ec="k"))
for x in [1.5, 2.5, 7.5]: ax.axvline(x, color="w", lw=1.5)
for x, t in [(0.5, "signatures"), (2, "ligand"), (5, "ANGPTL4 receptors"), (12, "tier-1 targets")]:
    ax.text(x, -0.62, t, ha="center", va="bottom", fontsize=6.3, style="italic")
for k in range(8, 17):
    ax.get_xticklabels()[k].set_color("#c0392b" if genes[k] != "CAMK2N1" else "#d4880f")
    ax.get_xticklabels()[k].set_fontweight("bold")
cb = fig.colorbar(im, ax=ax, fraction=0.02, pad=0.01); cb.set_label("mean z", fontsize=6); cb.ax.tick_params(labelsize=6)
L(ax, "B", -0.2, 1.1)
fig.savefig(f"{OUT}/Figure_S5.png", dpi=300, bbox_inches="tight"); fig.savefig(f"{OUT}/Figure_S5.pdf", bbox_inches="tight")
print("ok")
