"""21 - Main Figure 5: the ANGPTL4 niche axis (hypoxia adjustment, anatomy, cellular source, receptors, protein).
Inputs: ncdeg_hypoxia/ (09c), ivy_out/ (17), gbmap_out/ (18), cptac_out/cptac_results.json (19)."""
import sys, glob, json, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.linewidth": 0.6, "pdf.fonttype": 42, "axes.titlesize": 7})
FIG = sys.argv[1] if len(sys.argv) > 1 else "manuscript_figs"
RED, ORG, GREY, BLUE = "#c0392b", "#e08a3c", "#9a9a9a", "#2A6FB0"
def L(ax, s, dx=-0.16, dy=1.08): ax.text(dx, dy, s, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom")
fig = plt.figure(figsize=(7.2, 8.4))
gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 1.05], hspace=0.75, wspace=0.95)

# A hypoxia adjustment
W = pd.read_csv(glob.glob("ncdeg_hypoxia/ncdeg_hypoxia_adjusted_*.csv")[0], index_col=0)
genes = ["ANGPTL4", "VMP1", "NAMPT", "PDPN", "CAV1", "VEGFA", "ADM", "PGK1", "LDHA"]
ax = fig.add_subplot(gs[0, 0]); x = np.arange(len(genes)); w = 0.27
for i, (m, c, lb) in enumerate([("unadj", "#bdbdbd", "unadjusted"), ("H", ORG, "+ hypoxia"), ("HS", RED, "+ hypoxia + state")]):
    ax.bar(x + (i - 1) * w, W.loc[genes, f"mean_beta_{m}"], w, color=c, label=lb)
    for j, g in enumerate(genes):
        if m != "unadj" and W.loc[g, f"p_{m}"] < 0.05: ax.text(x[j] + (i - 1) * w, W.loc[g, f"mean_beta_{m}"] + 0.01, "*", ha="center", fontsize=6)
ax.axhline(0, color="k", lw=0.5); ax.axvline(4.5, color="0.75", ls="--", lw=0.6)
ax.set_xticks(x); ax.set_xticklabels(genes, rotation=50, ha="right")
for k in range(5, len(genes)): ax.get_xticklabels()[k].set_color("0.45")
ax.set_ylabel("niche-6 effect (mean β, 15 sections)"); ax.set_ylim(0, 0.68)
ax.legend(frameon=False, fontsize=5.5, ncol=3, loc="lower center", bbox_to_anchor=(0.5, 1.0))
ax.text(2, 0.62, "niche genes", ha="center", style="italic", fontsize=6); ax.text(6.5, 0.62, "canonical hypoxia", ha="center", style="italic", fontsize=6)
L(ax, "A", dy=1.12)

# B IvyGAP region boxplots
ivy = pd.read_csv("ivy_out/ivygap_sample_scores.csv", index_col=0)
order = ["LE", "IT", "CT", "CTmvp", "CTpan"]; lab = ["LE", "IT", "CT", "MVP", "PAN"]
col = {"LE": "#8fb3d9", "IT": "#6a9fcf", "CT": "#9a9a9a", "CTmvp": ORG, "CTpan": RED}
sub = gs[0, 1].subgridspec(1, 2, wspace=0.45)
rng = np.random.default_rng(0)
for k, (c, t) in enumerate([("HypoxicMES_N6", "N6 signature"), ("ANGPTL4", "ANGPTL4")]):
    ax = fig.add_subplot(sub[k]); data = [ivy.loc[ivy.region == r, c].values for r in order]
    bp = ax.boxplot(data, widths=0.6, patch_artist=True, showfliers=False, medianprops=dict(color="k", lw=0.8))
    for p, r in zip(bp["boxes"], order): p.set_facecolor(col[r]); p.set_alpha(0.55)
    for j, (r, d) in enumerate(zip(order, data)): ax.scatter(j + 1 + rng.uniform(-0.18, 0.18, len(d)), d, s=2.5, color=col[r], lw=0)
    ax.set_xticks(range(1, 6)); ax.set_xticklabels(lab, rotation=90, fontsize=6); ax.set_title(t)
    if k == 0: ax.set_ylabel("score (z), IvyGAP"); L(ax, "B", dx=-0.5, dy=1.12)

# C IvyGAP receptors heatmap
means = pd.read_csv("ivy_out/ivygap_region_means.csv", index_col=0)
rg = ["ANGPTL4", "SDC4", "CD44", "SDC2", "ITGA5", "ITGB1"]
M = means.loc[rg, order].astype(float)
ax = fig.add_subplot(gs[1, 0]); im = ax.imshow(M.values, cmap="RdBu_r", vmin=-1.5, vmax=1.5, aspect="auto")
ax.set_yticks(range(len(rg))); ax.set_yticklabels(rg); ax.set_xticks(range(5)); ax.set_xticklabels(["LE", "IT", "CT", "MVP", "PAN"], fontsize=6.5)
for i, g in enumerate(rg):
    j = order.index(M.loc[g].idxmax()); ax.add_patch(plt.Rectangle((j - .5, i - .5), 1, 1, fill=False, lw=1, ec="k"))
ax.axhline(0.5, color="w", lw=1.5)
cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02); cb.set_label("mean z", fontsize=6); cb.ax.tick_params(labelsize=5.5)
ax.set_title("IvyGAP: ligand vs receptors by region", loc="left"); L(ax, "C", dx=-0.3)

# D GBmap ANGPTL4 % expressing
pct = pd.read_csv("gbmap_out/gbmap_state_pct_expressing.csv", index_col=0)
states = [("MES-like hypoxia/MHC", "MES-like, hypoxic"), ("MES-like hypoxia independent", "MES-like, non-hypoxic"), ("AC-like", "AC-like"),
          ("OPC-like", "OPC-like"), ("NPC-like neural", "NPC-like"), ("TAM-BDM hypoxia/MES", "Macrophage, hypoxic"), ("Mono hypoxia", "Monocyte, hypoxic"),
          ("Astrocyte", "Astrocyte"), ("Pericyte", "Pericyte"), ("SMC", "Smooth muscle"), ("Endo arterial", "Endothelial, arterial")]
states = [(a, b) for a, b in states if a in pct.index]
ax = fig.add_subplot(gs[1, 1]); v = [pct.loc[a, "ANGPTL4"] for a, _ in states]
cl = [RED if a.startswith("MES") else (ORG if a in ("TAM-BDM hypoxia/MES", "Mono hypoxia") else ("#6a9fcf" if a in ("Pericyte", "SMC", "Endo arterial") else GREY)) for a, _ in states]
ax.barh(range(len(states)), v, color=cl, height=0.65); ax.set_yticks(range(len(states))); ax.set_yticklabels([b for _, b in states], fontsize=6); ax.invert_yaxis()
for i, x_ in enumerate(v): ax.text(x_ + 1, i, f"{x_:.0f}%", va="center", fontsize=5.5)
ax.set_xlabel("% cells expressing ANGPTL4 (GBmap)"); ax.set_title("Cellular source of ANGPTL4", loc="left"); L(ax, "D", dx=-0.55)

# E GBmap receptors by cell class
mean = pd.read_csv("gbmap_out/gbmap_state_mean_lognorm.csv", index_col=0)
cls = [("MES-like hypoxia/MHC", "MES-like hyp."), ("MES-like hypoxia independent", "MES-like indep."), ("TAM-BDM hypoxia/MES", "Macrophage hyp."),
       ("Endo arterial", "Endo arterial"), ("Tip-like", "Endo tip"), ("Pericyte", "Pericyte"), ("SMC", "SMC"), ("Astrocyte", "Astrocyte")]
cls = [(a, b) for a, b in cls if a in mean.index]
rec = ["SDC2", "SDC4", "ITGA5", "ITGB1", "CD44", "CDH11"]
Z = mean.loc[[a for a, _ in cls], rec]; Zs = (Z - Z.min()) / (Z.max() - Z.min() + 1e-9); Pp = pct.loc[[a for a, _ in cls], rec]
ax = fig.add_subplot(gs[2, 0]); yy, xx = np.meshgrid(range(len(cls)), range(len(rec)), indexing="ij")
sc = ax.scatter(xx.ravel(), yy.ravel(), s=Pp.values.ravel() / 100 * 90 + 1, c=Zs.values.ravel(), cmap="Reds", vmin=0, vmax=1, edgecolors="0.6", linewidths=0.3)
ax.set_xticks(range(len(rec))); ax.set_xticklabels(rec, rotation=45, ha="right"); ax.set_yticks(range(len(cls))); ax.set_yticklabels([b for _, b in cls], fontsize=6); ax.invert_yaxis()
ax.set_xlim(-0.6, len(rec) - 0.4); ax.set_ylim(len(cls) - 0.4, -0.6)
cb = fig.colorbar(sc, ax=ax, fraction=0.04, pad=0.02); cb.set_label("scaled mean", fontsize=6); cb.ax.tick_params(labelsize=5.5)
ax.set_title("ANGPTL4 receptors by cell type (size = % expressing)", loc="left"); L(ax, "E", dx=-0.42)

# F CPTAC protein
R = json.load(open("cptac_out/cptac_results.json"))
items = [("Hypoxia MP†", R["ANGPTL4p_vs"]["HYPnoA_p"]), ("N6 signature†", R["ANGPTL4p_vs"]["N6noA_p"]), ("N3 myeloid", R["ANGPTL4p_vs"]["N3_p"]),
         ("MES MP", R["ANGPTL4p_vs"]["MES_p"]), ("ITGA5", R["ANGPTL4p_receptors"]["ITGA5"]), ("ITGB1", R["ANGPTL4p_receptors"]["ITGB1"]),
         ("CD44", R["ANGPTL4p_receptors"]["CD44"]), ("SDC2", R["ANGPTL4p_receptors"]["SDC2"]), ("SDC4", R["ANGPTL4p_receptors"]["SDC4"])]
ax = fig.add_subplot(gs[2, 1]); v = [x["rho"] for _, x in items]
ax.barh(range(len(items)), v, color=[BLUE] * 4 + [RED] * 5, height=0.65); ax.axvline(0, color="k", lw=0.5); ax.axhline(3.5, color="0.8", lw=0.6)
ax.set_yticks(range(len(items))); ax.set_yticklabels([n for n, _ in items], fontsize=6); ax.invert_yaxis(); ax.set_xlim(-0.1, 0.95)
for i, (_, x_) in enumerate(items): ax.text(max(x_["rho"], 0) + 0.02, i, "***" if x_["p"] < 1e-3 else "**" if x_["p"] < 0.01 else "*" if x_["p"] < 0.05 else "", va="center", fontsize=6)
ax.set_xlabel("Spearman ρ with ANGPTL4 protein (CPTAC, n = 99)"); ax.set_title("Protein level", loc="left"); L(ax, "F", dx=-0.42)
fig.savefig(f"{FIG}/Figure_5_new.png", dpi=300, bbox_inches="tight"); fig.savefig(f"{FIG}/Figure_5_new.pdf", bbox_inches="tight"); print("ok")
