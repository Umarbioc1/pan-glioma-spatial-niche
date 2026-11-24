"""34_replot_bulk_panels.py -- re-draw the bulk, survival, progression and immunotherapy panels
from the saved analysis tables, for the publication restyle. No analysis is repeated: every value
is read from the tables written by scripts 11, 13, 14b and 15. Descriptive panel headings are
dropped; identity labels and statistics are kept (log-rank P moved inside the axes).
Usage: python3 33_pubstyle.py long 34_replot_bulk_panels.py <outdir>
"""
import sys, os, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from scipy.stats import mannwhitneyu

OUT = sys.argv[1] if len(sys.argv) > 1 else "pub_panels"; os.makedirs(OUT, exist_ok=True)
SIG = ["HypoxicMES_N6", "Myeloid_N3", "Neuronal_N1", "Oligodendroglial_N10"]
RED, TEAL = "#B5172F", "#2A8C82"
def save(fig, name, dpi=600):
    fig.savefig(f"{OUT}/{name}.png", dpi=dpi, bbox_inches="tight"); plt.close(fig); print("saved", name)

# --- TCGA scores by type and grade (script 11) ---
D = pd.read_csv("bulk_out/tcga_analysis_table.csv")
A = D[D.who == "Astrocytoma, IDH-mut"]
cats = [("Oligo", D.who == "Oligodendroglioma"), ("Astro G2", (A.grade == 2).reindex(D.index, fill_value=False)),
        ("Astro G3", (A.grade == 3).reindex(D.index, fill_value=False)), ("Astro G4", (A.grade == 4).reindex(D.index, fill_value=False)),
        ("IDH-wt LGG", D.who == "IDH-wt lower-grade"), ("GBM IDH-wt", D.who == "GBM, IDH-wt")]
fig, ax = plt.subplots(1, 4, figsize=(16, 4.2))
for a, s in zip(ax, SIG):
    a.boxplot([D.loc[m, s] for _, m in cats], tick_labels=[c for c, _ in cats], showfliers=False)
    a.set_title(s, fontsize=10); a.tick_params(axis="x", rotation=45)
fig.tight_layout(); save(fig, "tcga_scores_by_type_grade")

# --- CGGA scores by type and grade (script 13) ---
C = pd.read_csv("cgga_out/cgga_analysis_table.csv"); P = C[C.PRS_type == "Primary"]
cats = [("Oligo", P.who == "Oligodendroglioma"), ("Astro G2", (P.who == "Astrocytoma, IDH-mut") & (P.grade == 2)),
        ("Astro G3", (P.who == "Astrocytoma, IDH-mut") & (P.grade == 3)), ("Astro G4", (P.who == "Astrocytoma, IDH-mut") & (P.grade == 4)),
        ("IDH-wt LGG", P.who == "IDH-wt lower-grade"), ("GBM IDH-wt", P.who == "GBM, IDH-wt")]
fig, ax = plt.subplots(1, 4, figsize=(16, 4.2))
for a, s in zip(ax, SIG):
    a.boxplot([P.loc[m, s] for _, m in cats], tick_labels=[f"{c}\n(n={m.sum()})" for c, m in cats], showfliers=False)
    a.set_title(s, fontsize=10); a.tick_params(axis="x", rotation=45, labelsize=8)
fig.tight_layout(); save(fig, "cgga_scores_by_type_grade")

# --- CGGA Kaplan-Meier within WHO types (script 13) ---
fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
for a, g in zip(ax, ["Astrocytoma, IDH-mut", "Oligodendroglioma", "GBM, IDH-wt"]):
    df = P[P.who == g]; hi = df.HypoxicMES_N6 > df.groupby("cohort").HypoxicMES_N6.transform("median")
    for lab, m, col in [("High", hi, RED), ("Low", ~hi, TEAL)]:
        KaplanMeierFitter().fit(df.OS_MONTHS[m], df.OS_event[m], label=f"{lab} (n={m.sum()})").plot_survival_function(ax=a, ci_show=False, color=col)
    p = logrank_test(df.OS_MONTHS[hi], df.OS_MONTHS[~hi], df.OS_event[hi], df.OS_event[~hi]).p_value
    a.text(0.97, 0.93, f"{g}\nlog-rank P = {p:.2g}", transform=a.transAxes, ha="right", va="top", fontsize=9)
    a.set_xlabel("Overall survival (months)"); a.set_ylabel("Survival probability")
fig.tight_layout(); save(fig, "cgga_KM_HypoxicMES_N6")

# --- GSE121810 arms (script 15) ---
sc = pd.read_csv("icb_out/gse121810_scores.csv", index_col=0)
R = pd.read_csv("icb_out/gse121810_arm_comparison.csv").set_index("score")
show = ["IFNG_Tcell_control", "CellCycle_control", "HypoxicMES_N6", "Myeloid_N3", "MP_MES", "MP_Hypoxia", "ANGPTL4", "Neuronal_N1"]
fig, ax = plt.subplots(1, len(show), figsize=(2.1 * len(show), 3.8))
for a, c in zip(ax, show):
    d = [sc.loc[sc.arm == x, c] for x in ["adjuvant", "neoadjuvant"]]
    a.boxplot(d, tick_labels=["adj.", "neo-adj."], showfliers=False)
    for i, v in enumerate(d): a.scatter(np.random.default_rng(0).normal(i + 1, 0.06, len(v)), v, s=10, c=RED if i else TEAL, zorder=3)
    a.set_title(f"{c}\nP = {R.loc[c, 'p_MWU']:.2g}", fontsize=8)
fig.tight_layout(); save(fig, "gse121810_arm_scores")

# --- GSE318238 progression maps (script 14b) ---
sc_ = pd.read_csv("prog/progression_spot_scores.csv.gz")
sc_["year"] = sc_["year"].astype(str)
fig, ax = plt.subplots(3, 4, figsize=(15, 10.5))
for j, c in enumerate(SIG):
    lim = np.nanpercentile(np.abs(sc_[c]), 98)
    for i, y in enumerate(["2015", "2017", "2019"]):
        d = sc_[sc_.year == y]; a = ax[i, j]
        a.scatter(d.array_col, -d.array_row, c=d[c], cmap="RdBu_r", vmin=-lim, vmax=lim, s=6 if y == "2017" else 14, lw=0, rasterized=True)
        a.set_aspect("equal"); a.set_xticks([]); a.set_yticks([])
        a.set_title(f"{y} · {c}", fontsize=9)
        a.text(0.5, -0.03, f"mean {d[c].mean():+.2f}", transform=a.transAxes, ha="center", va="top", fontsize=8)
fig.tight_layout(); save(fig, "progression_niche_score_maps", dpi=300)
