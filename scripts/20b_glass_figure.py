"""20b - Figure S9: niche signatures through recurrence in GLASS."""
import sys, numpy as np, pandas as pd, matplotlib
from scipy.stats import kruskal
matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7, "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.6, "pdf.fonttype": 42})
IN = sys.argv[1] if len(sys.argv) > 1 else "glass_out"; FIG = sys.argv[2] if len(sys.argv) > 2 else "manuscript_figs"
P = pd.read_csv(f"{IN}/glass_pairs_scores.csv"); R = pd.read_csv(f"{IN}/glass_paired_changes.csv")
SUB = {"IDHmut-codel": "Oligo-\ndendroglioma", "IDHmut-noncodel": "Astrocytoma\nIDH-mut", "IDHwt": "IDH-wildtype"}
COL = {"IDHmut-codel": "#6a9fcf", "IDHmut-noncodel": "#e08a3c", "IDHwt": "#c0392b"}
fig = plt.figure(figsize=(7.2, 6.0)); gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.05], width_ratios=[1, 1.2], wspace=0.35, hspace=0.75)
ax = fig.add_subplot(gs[0, 0]); order = ["IDHmut-codel", "IDHmut-noncodel", "IDHwt"]
data = [P.loc[P.subtype == s, "HypoxicMES_N6_a"].dropna() for s in order]
bp = ax.boxplot(data, widths=0.6, patch_artist=True, showfliers=False, medianprops=dict(color="k", lw=0.8))
for p_, s in zip(bp["boxes"], order): p_.set_facecolor(COL[s]); p_.set_alpha(0.5)
rng = np.random.default_rng(0)
for i, (s, d) in enumerate(zip(order, data)): ax.scatter(i + 1 + rng.uniform(-0.18, 0.18, len(d)), d, s=3, color=COL[s], lw=0)
ax.set_xticks([1, 2, 3]); ax.set_xticklabels([SUB[s] for s in order], fontsize=6); ax.set_ylabel("N6 score, earlier tumour")
ax.set_title(f"Type ordering (Kruskal p = {kruskal(*data).pvalue:.0e})".replace("e-0", "e-"), loc="left", fontsize=6.5)
ax.text(-0.22, 1.06, "A", transform=ax.transAxes, fontsize=10, fontweight="bold")
ax = fig.add_subplot(gs[1, :])
scs = ["HypoxicMES_N6", "N6_resid_hypoxia", "ANGPTL4", "ANGPTL4_resid_hypoxia", "Myeloid_N3", "MP_Hypoxia", "MP_MES", "Neuronal_N1", "Oligodendroglial_N10"]
lab = ["N6", "N6 | hypoxia", "ANGPTL4", "ANGPTL4 | hypoxia", "N3 myeloid", "Hypoxia MP", "MES MP", "N1 neuronal", "N10 oligodendr."]
grp = ["IDHmut-codel", "IDHmut-noncodel", "IDHmut-noncodel, grade increase", "IDHwt"]
glab = ["Oligodendroglioma", "Astrocytoma", "Astro, grade ↑", "IDH-wildtype"]
M = np.full((len(grp), len(scs)), np.nan); Pv = M.copy(); F = M.copy(); N = []
for i, g in enumerate(grp):
    for j, c in enumerate(scs):
        r = R[(R.group == g) & (R.score == c)]
        if len(r): M[i, j], Pv[i, j], F[i, j] = r.median_delta.iloc[0], r.p.iloc[0], r.fdr.iloc[0]
    N.append(int(R[(R.group == g)].n_pairs.iloc[0]))
im = ax.imshow(M, cmap="RdBu_r", vmin=-0.4, vmax=0.4, aspect="auto")
for i in range(len(grp)):
    for j in range(len(scs)):
        t = "**" if F[i, j] < 0.05 else ("*" if Pv[i, j] < 0.05 else "")
        ax.text(j, i, t, ha="center", va="center", fontsize=7)
ax.set_xticks(range(len(scs))); ax.set_xticklabels(lab, rotation=55, ha="right", fontsize=6)
ax.set_yticks(range(len(grp))); ax.set_yticklabels([f"{g} ({n} pairs)" for g, n in zip(glab, N)], fontsize=6.5); ax.set_xticklabels(lab, rotation=40, ha="right", fontsize=6.5)
cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01); cb.set_label("median change, later − earlier", fontsize=6); cb.ax.tick_params(labelsize=5.5)
ax.set_title("Paired change at recurrence (* p<0.05, ** FDR<0.05)", loc="left", fontsize=6.5)
ax.text(-0.16, 1.06, "C", transform=ax.transAxes, fontsize=10, fontweight="bold")
ax = fig.add_subplot(gs[0, 1])
A = P[(P.subtype == "IDHmut-noncodel") & P.grade_change.isin(["increase", "stable"])]
for k, (gc, x0, col) in enumerate([("stable", 0, "#9a9a9a"), ("increase", 2, "#e08a3c")]):
    d = A[A.grade_change == gc]
    for _, r in d.iterrows(): ax.plot([x0, x0 + 1], [r.ANGPTL4_a, r.ANGPTL4_b], color=col, lw=0.6, alpha=0.6)
    ax.plot([x0, x0 + 1], [d.ANGPTL4_a.median(), d.ANGPTL4_b.median()], color="k", lw=1.6)
    rr = R[(R.group == f"IDHmut-noncodel, grade {gc}") & (R.score == "ANGPTL4")].iloc[0]
    ax.text(x0 + 0.5, ax.get_ylim()[1] if False else 2.6, f"{gc} (n={len(d)})\n{rr.frac_up*100:.0f}% up, p={rr.p:.2f}", ha="center", fontsize=6)
ax.set_xticks([0, 1, 2, 3]); ax.set_xticklabels(["earlier", "later", "earlier", "later"], fontsize=6); ax.set_ylim(-2.6, 3.2)
ax.set_ylabel("ANGPTL4 (z)"); ax.set_title("Astrocytoma: ANGPTL4 by grade change", loc="left", fontsize=6.5)
ax.text(-0.18, 1.06, "B", transform=ax.transAxes, fontsize=10, fontweight="bold")
fig.savefig(f"{FIG}/Figure_S9.png", dpi=300, bbox_inches="tight"); fig.savefig(f"{FIG}/Figure_S9.pdf", bbox_inches="tight"); print("ok")
