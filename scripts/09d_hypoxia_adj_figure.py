"""09d - Figure S6: niche-6 effects in malignant spots before and after hypoxia (H) and hypoxia + state (HS) adjustment."""
import sys, glob, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.linewidth": 0.6, "pdf.fonttype": 42})
IN = sys.argv[1] if len(sys.argv) > 1 else "ncdeg_hypoxia"; FIG = sys.argv[2] if len(sys.argv) > 2 else "manuscript_figs"
W = pd.read_csv(glob.glob(f"{IN}/ncdeg_hypoxia_adjusted_*.csv")[0], index_col=0)
genes = ["ANGPTL4", "NAMPT", "PDPN", "VMP1", "CAV1", "CA12", "NDRG1", "NRN1", "ARRDC3", "VEGFA", "ADM", "PGK1", "LDHA"]
genes = [g for g in genes if g in W.index]
fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.9), gridspec_kw=dict(width_ratios=[1.7, 1], wspace=0.35))
ax = axs[0]; x = np.arange(len(genes)); w = 0.27
cols = {"unadj": "#bdbdbd", "H": "#e08a3c", "HS": "#c0392b"}
labs = {"unadj": "unadjusted", "H": "+ hypoxia", "HS": "+ hypoxia + state"}
for i, m in enumerate(["unadj", "H", "HS"]):
    ax.bar(x + (i - 1) * w, W.loc[genes, f"mean_beta_{m}"], w, color=cols[m], label=labs[m])
    for j, g in enumerate(genes):
        if W.loc[g, f"p_{m}"] < 0.05 and m != "unadj":
            ax.text(x[j] + (i - 1) * w, W.loc[g, f"mean_beta_{m}"] + 0.015, "*", ha="center", fontsize=7)
ax.axhline(0, color="k", lw=0.5); ax.set_xticks(x); ax.set_xticklabels(genes, rotation=55, ha="right")
for k in range(9, len(genes)): ax.get_xticklabels()[k].set_color("0.45")
ax.axvline(8.5, color="0.8", lw=0.6, ls="--")
ax.set_ylim(top=0.66); ax.text(4, 0.635, "tier-1 targets", ha="center", style="italic", fontsize=6.5)
ax.text(10.5, 0.635, "canonical hypoxia genes", ha="center", style="italic", fontsize=6.5)
ax.set_ylabel("niche-6 effect, mean β across sections"); ax.legend(frameon=False, fontsize=6, loc="lower center", bbox_to_anchor=(0.5, 1.02), ncol=3)
ax.text(-0.1, 1.1, "A", transform=ax.transAxes, fontsize=10, fontweight="bold", ha="right")
ax = axs[1]
up = W[W.mean_beta_HS > 0].sort_values("p_HS")
r = -np.log10(up.p_HS.values)
ax.scatter(np.arange(1, len(r) + 1), r, s=3, color="0.6", lw=0)
for g in ["ANGPTL4", "VMP1", "CAV1"]:
    k = list(up.index).index(g)
    ax.scatter(k + 1, r[k], s=14, color="#c0392b", zorder=3); ax.text((k + 1) * 1.6, r[k], f"{g} (rank {k + 1})", va="center", fontsize=6.5, color="#c0392b")
ax.set_xscale("log"); ax.set_xlabel(f"rank among {len(up):,} genes higher in niche 6\n(hypoxia + state model)"); ax.set_ylabel("−log10 p")
ax.text(-0.18, 1.1, "B", transform=ax.transAxes, fontsize=10, fontweight="bold", ha="right")
fig.savefig(f"{FIG}/Figure_S6.png", dpi=300, bbox_inches="tight"); fig.savefig(f"{FIG}/Figure_S6.pdf", bbox_inches="tight"); print("ok")
