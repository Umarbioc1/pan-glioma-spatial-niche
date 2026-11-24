"""
16_article_figures.py -- assemble journal-style multi-panel figures (Fig 1-6, Fig S1-S4) from result tables and maps.
Output: 06_manuscript/figures/{Figure_N}.png (300 dpi) and .pdf. Width 183 mm (double column).
"""
import os, json, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from adjustText import adjust_text
from matplotlib.gridspec import GridSpec
from PIL import Image
Image.MAX_IMAGE_PIXELS = None
OUT = "manuscript_figs"; os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7, "axes.titlesize": 7.5, "axes.labelsize": 7,
                     "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6.5, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.linewidth": 0.6, "pdf.fonttype": 42})
GBM, IDH, GREY, INK = "#B5172F", "#2A8C82", "#9AA3AE", "#1E1E24"
W = 7.2  # inches

def letter(ax, L, dx=-0.02, dy=1.02):
    ax.text(dx, dy, L, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom", ha="right")

def img(ax, path, crop=None, blank=None):
    im = Image.open(path).convert("RGB")
    if crop: im = im.crop(crop)
    if blank:
        from PIL import ImageDraw
        if blank == "br":  # bottom-right legend block, as a fraction of the panel
            w, h = im.size; blank = (round(w * 1330 / 1529), round(h * 1180 / 1425), w, h)
        ImageDraw.Draw(im).rectangle(blank, fill="white")
    ax.imshow(im); ax.axis("off")

def save(fig, name):
    fig.savefig(f"{OUT}/{name}.png", dpi=300, bbox_inches="tight"); fig.savefig(f"{OUT}/{name}.pdf", bbox_inches="tight"); plt.close(fig)
    print("saved", name)

N_ORDER = [f"Niche_{i}" for i in range(1, 11)]
NLAB = {"Niche_1": "N1 Neuronal", "Niche_2": "N2 Vascular", "Niche_3": "N3 Myeloid", "Niche_4": "N4 OPC-like",
        "Niche_5": "N5 Astrocytic", "Niche_6": "N6 Hypoxic MES", "Niche_7": "N7 Brain/infiltration",
        "Niche_8": "N8 NPC-like", "Niche_9": "N9 AC-like", "Niche_10": "N10 Oligodendroglial"}

# ---------------- Figure 1 ----------------
fig = plt.figure(figsize=(W, 7.4)); gs = GridSpec(2, 2, figure=fig, height_ratios=[1, 1.6], width_ratios=[1, 1], hspace=0.3, wspace=0.6)
ax = fig.add_subplot(gs[0, 0])
cats = ["GBM IDH-wt", "Astrocytoma IDH-mut", "Oligodendroglioma", "IDH-mut (Ravi)", "Cortex control"]
g = [13, 3, 3, 0, 0]; r = [16, 0, 0, 2, 3]
y = np.arange(len(cats))[::-1]
ax.barh(y + 0.18, g, 0.36, color="#5B6C8F", label="Greenwald 2024 (19)"); ax.barh(y - 0.18, r, 0.36, color="#D9A441", label="Ravi 2022 (21)")
for yi, a, b in zip(y, g, r):
    if a: ax.text(a + 0.3, yi + 0.18, str(a), va="center", fontsize=6.5)
    if b: ax.text(b + 0.3, yi - 0.18, str(b), va="center", fontsize=6.5)
ax.set_yticks(y); ax.set_yticklabels(cats); ax.set_xlabel("Visium sections passing QC"); ax.legend(frameon=False, loc="lower right")
ax.set_title("Atlas: 40 sections, 109,415 spots", loc="left")
letter(ax, "A", dx=-0.42)
ax = fig.add_subplot(gs[0, 1])
import glob as _g
_J = []
for _f in _g.glob("deconv_out/spots/*_proportions.csv.gz"):
    _s = os.path.basename(_f).replace("_proportions.csv.gz", "")
    _pf = [q for q in (f"step2a_out/spots/{_s}_program_scores.csv.gz", f"step2a_ravi/spots/{_s}_program_scores.csv.gz") if os.path.exists(q)][0]
    _d = pd.read_csv(_f, index_col=0); _p = pd.read_csv(_pf, index_col=0)
    _J.append(_d.join(_p[["Neuron", "Oligodendrocyte", "Microglia", "Mono_Macrophage", "MES_like"]].add_suffix("_p"), how="inner"))
_J = pd.concat(_J)  # pooled Spearman rho, deconvolved fraction vs marker-program score, all spots
_J["my"] = _J[["TAM-BDM", "TAM-MG", "Mono", "DC"]].sum(1); _J["myp"] = _J[["Microglia_p", "Mono_Macrophage_p"]].mean(1)
types = ["Neuron", "Oligodendrocyte", "All myeloid", "MES-like"]
rv = [round(_J[a].corr(_J[b], method="spearman"), 2) for a, b in [("Neuron", "Neuron_p"), ("Oligodendrocyte", "Oligodendrocyte_p"), ("my", "myp"), ("MES-like", "MES_like_p")]]
print("Fig 1B rho:", rv)
ax.barh(types[::-1], rv[::-1], color=IDH, height=0.55)
for i, v in enumerate(rv[::-1]): ax.text(v + 0.01, i, f"{v:.2f}", va="center", fontsize=6.5)
ax.set_xlim(0, 1); ax.set_xlabel("Spearman ρ, deconvolved fraction vs marker score")
ax.set_title("Deconvolution validation (median fit r = 0.90)", loc="left")
letter(ax, "B", dx=-0.32)
ax = fig.add_subplot(gs[1, :]); img(ax, "cnv_out/compartment_maps.png")
ax.set_title("Copy-number-defined compartments per section (Malignant / TME-rich / Brain)", loc="left"); letter(ax, "C", dx=0.0)
save(fig, "Figure_1")

# ---------------- Figure 2 ----------------
P = pd.read_csv("niche_deconv/niche_program_profile.csv", index_col=0).loc[N_ORDER]
F = pd.read_csv("niche_deconv/niche_frequency_by_subtype.csv", index_col=0)[N_ORDER] * 100
fig = plt.figure(figsize=(W, 9.0)); gs = GridSpec(2, 2, figure=fig, height_ratios=[1, 1.45], width_ratios=[1.75, 1], hspace=0.3, wspace=0.08)
ax = fig.add_subplot(gs[0, 0])
im = ax.imshow(P.values, cmap="RdBu_r", vmin=-1.5, vmax=1.5, aspect="auto")
ax.set_yticks(range(10)); ax.set_yticklabels([NLAB[n] for n in N_ORDER]); ax.set_xticks(range(P.shape[1]))
ax.set_xticklabels(P.columns, rotation=90); ax.tick_params(length=0)
for s in ax.spines.values(): s.set_visible(False)
cb = fig.colorbar(im, ax=ax, orientation="horizontal", fraction=0.04, pad=0.33, aspect=40); cb.set_label("mean z (deconvolved fraction)", fontsize=6); cb.ax.tick_params(labelsize=6)
ax.set_title("Cell-type composition of the 10 niches", loc="left"); letter(ax, "A", dx=-0.3)
ax = fig.add_subplot(gs[0, 1])
yy = np.arange(10)[::-1]
ax.barh(yy + 0.2, F.loc["GBM"].values, 0.4, color=GBM, label="GBM IDH-wt"); ax.barh(yy - 0.2, F.loc["IDH-mut"].values, 0.4, color=IDH, label="IDH-mutant")
ax.set_yticks(yy); ax.set_yticklabels([f"N{i}" for i in range(1, 11)]); ax.set_xlabel("% of spots"); ax.legend(frameon=False, loc="lower right")
ax.text(F.loc["GBM", "Niche_6"] + 0.8, yy[5] + 0.2, "10.4% vs 1.2%", fontsize=6.5, color=GBM, va="center", fontweight="bold")
ax.set_title("Niche frequency by subtype", loc="left"); letter(ax, "B", dx=-0.04)
ax = fig.add_subplot(gs[1, :]); img(ax, "niche_deconv/niche_maps.png", blank="br"); ax.set_title("Niche maps across tumour sections", loc="left"); letter(ax, "C", dx=0.0)
tab = plt.get_cmap("tab10")
from matplotlib.lines import Line2D
ax.legend(handles=[Line2D([], [], marker="o", ls="", color=tab(i), label=NLAB[f"Niche_{i+1}"], markersize=5) for i in range(10)],
          loc="center left", bbox_to_anchor=(1.0, 0.5), frameon=False, fontsize=6.5, title="Niche", title_fontsize=7)
save(fig, "Figure_2")

# ---------------- Figure 3 ----------------
B = pd.read_csv("lr_out/lr_significant_by_subtype.csv"); S = pd.read_csv("lr_out/lr_niche_significant.csv")
b6 = B[B.niche == "Niche_6"].merge(S[S.niche == "Niche_6"][["pair", "median_auc"]], on="pair").sort_values("median_auc_GBM")
def grp(p):
    l = p.split("->")[0]
    return "Hypoxia" if l in ("ANGPTL4", "NAMPT") else ("Macrophage" if l in ("SPP1", "MIF") else "Matrix")
gc = {"Hypoxia": "#B5172F", "Matrix": "#5B6C8F", "Macrophage": "#D9A441"}
fig = plt.figure(figsize=(W, 6.2)); gs = GridSpec(1, 2, figure=fig, width_ratios=[1.25, 1], wspace=0.55)
ax = fig.add_subplot(gs[0, 0]); img(ax, "lr_out/lr_niche_dotplot.png"); letter(ax, "A", dx=0.02, dy=1.0)
ax = fig.add_subplot(gs[0, 1]); yy = np.arange(len(b6))
for i, (_, r) in enumerate(b6.iterrows()):
    if pd.notna(r["median_auc_IDH-mut"]): ax.plot([r["median_auc_IDH-mut"], r.median_auc_GBM], [i, i], c="#D5D9DF", lw=1.2, zorder=1)
ax.scatter(b6.median_auc_GBM, yy, c=GBM, s=16, zorder=3, label="GBM sections (n = 13-14)")
ax.scatter(b6["median_auc_IDH-mut"], yy, c=IDH, s=16, zorder=3, label="IDH-mut sections (n = 1-2)")
ax.set_yticks(yy); ax.set_yticklabels([p.replace("->", " → ").replace("_", "–") for p in b6.pair], fontsize=6)
for t, p in zip(ax.get_yticklabels(), b6.pair): t.set_color(gc[grp(p)])
ax.axvline(0.5, c=GREY, ls="--", lw=0.6); ax.set_xlabel("median AUC, niche 6 vs rest of section"); ax.set_xlim(0.45, 0.8)
ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(0.0, -0.08), fontsize=6, ncol=2)
ax.set_title("22 niche-6 pairs: GBM vs IDH-mutant\n(label colour: hypoxia / matrix / macrophage ligand)", loc="left"); letter(ax, "B", dx=-0.45, dy=1.0)
save(fig, "Figure_3")

# ---------------- Figure 4 ----------------
J = pd.read_csv("ncdeg_out/niche6_malignant_unadjusted_vs_stateadjusted.csv")
E = pd.read_csv("ncdeg_out/ncdeg_metaprogram_enrichment_sig.csv")
E = E[(E.niche == "Niche_6") & E.meta_program.str.startswith("Cancer|")].sort_values("fdr").head(7)
fig = plt.figure(figsize=(W, 5.6)); gs = GridSpec(2, 2, figure=fig, height_ratios=[1.25, 0.8], hspace=0.45, wspace=0.35)
ax = fig.add_subplot(gs[0, 0])
c = np.where(~J.NCDEG, "#D5D9DF", np.where(J.mean_beta > 0, GBM, "#2A6FB0"))
ax.scatter(J.mean_beta, -np.log10(J.fdr), s=5, c=c, lw=0)
lab4 = J[J.NCDEG & J.gene.isin(["ANGPTL4", "CHI3L1", "NDRG1", "VEGFA", "ADM", "CA12", "CD44", "PDPN", "LDHA", "BCAN", "OLIG1", "NCAN", "CSPG5"])]
tx = [ax.text(r.mean_beta, -np.log10(r.fdr), r.gene, fontsize=5.5, fontweight="bold" if r.gene == "ANGPTL4" else "normal") for _, r in lab4.iterrows()]
adjust_text(tx, ax=ax, arrowprops=dict(arrowstyle="-", color=GREY, lw=0.4),
            expand_text=(1.25, 1.35), expand_points=(1.3, 1.4), force_text=(0.4, 0.8),
            force_points=(0.3, 0.6), max_move=40, time_lim=12)
ax.axhline(-np.log10(0.05), c=GREY, ls="--", lw=0.6); [ax.axvline(v, c=GREY, ls=":", lw=0.6) for v in (-0.2, 0.2)]
ax.set_xlabel("niche-6 effect in malignant-rich spots (mean β)"); ax.set_ylabel("−log10 FDR (section-level)")
ax.set_title("47 malignant NCDEGs in niche 6 (32 up, 15 down)", loc="left"); letter(ax, "A", dx=-0.15)
ax = fig.add_subplot(gs[0, 1]); S6 = J[J.NCDEG]
ax.scatter(S6.mean_beta, S6.mean_beta_adj, c=np.where(S6.mean_beta > 0, GBM, "#2A6FB0"), s=14)
ax.margins(0.16)
OFF = {"ANGPTL4": (-26, 20, "right"), "NDRG1": (-8, 24, "right"), "VEGFA": (24, 10, "left"),
       "CHI3L1": (-10, -20, "right"), "BCAN": (-6, -20, "right")}  # fixed label offsets (points)
for _, r in S6.iterrows():
    if r.gene not in OFF: continue
    dx, dy, ha = OFF[r.gene]
    ax.annotate(r.gene, (r.mean_beta, r.mean_beta_adj), textcoords="offset points", xytext=(dx, dy),
                ha=ha, va="center", fontsize=5.8, fontweight="bold" if r.gene == "ANGPTL4" else "normal",
                arrowprops=dict(arrowstyle="-", color=GREY, lw=0.4, shrinkA=0, shrinkB=2))
ax.plot([-0.6, 0.8], [-0.6, 0.8], c=GREY, ls="--", lw=0.6); ax.axhline(0, c=INK, lw=0.4); ax.axvline(0, c=INK, lw=0.4)
ax.set_xlabel("effect, admixture-adjusted"); ax.set_ylabel("effect, + malignant-state adjusted")
ax.set_title("Most effects shrink after state adjustment;\nANGPTL4 stays up (14/15 sections, p = 1×10⁻⁴)", loc="left"); letter(ax, "B", dx=-0.17)
ax = fig.add_subplot(gs[1, :])
lab = [m.split("|")[1].split(" ", 1)[1] + f"  ({d})" for m, d in zip(E.meta_program, E.direction)]
ax.barh(lab[::-1], -np.log10(E.fdr.values[::-1]), color=np.where(E.direction.values[::-1] == "up", GBM, "#2A6FB0"), height=0.6)
ax.set_xlabel("−log10 FDR (hypergeometric overlap with Gavish cancer meta-programs)")
ax.set_title("Niche-6 NCDEGs map to hypoxia and MES programs (up) and astrocyte/OPC programs (down)", loc="left"); letter(ax, "C", dx=-0.28)
save(fig, "Figure_4")

# ---------------- Figure 5 ----------------
fig = plt.figure(figsize=(W, 8.0)); gs = GridSpec(4, 1, figure=fig, height_ratios=[0.9, 0.9, 0.95, 1.0], hspace=0.28)
for i, (p, L, t) in enumerate([("bulk_out/scores_by_type_grade.png", "A", "TCGA LGG + GBM (n = 618 classified tumours)"),
                               ("cgga_out/cgga_scores_by_type_grade.png", "B", "CGGA-693 + CGGA-325 primary tumours"),
                               ("cgga_out/meta_forest.png", "C", "Overall survival meta-analysis (1,159 patients, 521 deaths)"),
                               ("cgga_out/cgga_KM_HypoxicMES_N6.png", "D", "CGGA: niche-6 score within WHO types")]):
    ax = fig.add_subplot(gs[i, 0]); img(ax, p); ax.set_title(t, loc="left"); letter(ax, L, dx=0.0)
save(fig, "Figure_5")

# ---------------- Figure 6 ----------------
T = pd.read_csv("targets_out/niche6_target_prioritisation.csv")
sel = T[T.tier.astype(str).str.startswith("1") | T.gene.isin(["VEGFA", "LGALS3", "MIF", "CD44", "TNC", "ADM"])].copy()
sel = sel[sel.stateadj_beta.notna()].sort_values("stateadj_beta")
fig, axs = plt.subplots(1, 3, figsize=(W, 4.0), sharey=True, gridspec_kw={"width_ratios": [1.3, 0.8, 0.8], "wspace": 0.18})
yy = np.arange(len(sel)); tier1 = (sel.tier.astype(str) == "1").values
col = np.where(tier1, GBM, np.where(sel.tier.astype(str).str.startswith("1b").values, "#D9A441", GREY))
ax = axs[0]; ax.hlines(yy, 0, sel.stateadj_beta, color=col, lw=1.2); ax.scatter(sel.stateadj_beta, yy, c=col, s=22, zorder=3)
for yi, (b, p) in enumerate(zip(sel.stateadj_beta, sel.stateadj_p)): ax.text(max(b, 0) + 0.008, yi, f"$P$ = {p:.2g}", va="center", fontsize=5.5)
ax.axvline(0, c=INK, lw=0.4); ax.set_xlim(-0.01, 0.27); ax.set_yticks(yy); ax.set_yticklabels(sel.gene)
for t, k in zip(ax.get_yticklabels(), tier1): t.set_fontweight("bold" if k else "normal")
ax.set_xlabel("niche-6 effect after state adjustment (β)"); ax.set_title("Niche specificity", loc="left"); letter(ax, "A", dx=-0.18)
ax = axs[1]; ax.barh(yy, sel.OT_GBM_assoc, color="#5B6C8F", height=0.55); ax.set_xlabel("Open Targets GBM score")
ax.set_title("Prior GBM evidence", loc="left"); letter(ax, "B", dx=-0.02)
ax = axs[2]; ax.scatter(sel.DepMap_glioma_median_effect, yy, c=np.where(sel.DepMap_glioma_median_effect < -0.5, GBM, GREY), s=18)
ax.axvline(-1, c=GREY, ls="--", lw=0.6); ax.set_xlabel("DepMap glioma (Chronos)")
stage = {4: "approved", 3: "Ph3", 2.5: "Ph2/3", 2: "Ph2", 1.5: "Ph1/2", 1: "Ph1"}
for yi, (_, r) in enumerate(sel.iterrows()):
    if r.max_stage > 0 and r.gene not in ("CD74", "SDC4"): ax.text(0.55, yi, stage.get(r.max_stage, ""), va="center", fontsize=5.5, color=GBM)
ax.set_xlim(-1.3, 0.85); ax.set_title("Dependency · drug stage", loc="left"); letter(ax, "C", dx=-0.02)
fig.suptitle("Niche-resolved targets. Red/bold: tier 1 (niche-specific after state adjustment, tractable); amber: tier 1b (not tractable); grey: tier 2", x=0.01, ha="left", fontsize=7.5)
save(fig, "Figure_6")

# ---------------- Supplementary ----------------
def stack(paths, titles, name, h):
    fig = plt.figure(figsize=(W, h)); gs = GridSpec(len(paths), 1, figure=fig, hspace=0.18)
    for i, (p, t) in enumerate(zip(paths, titles)):
        ax = fig.add_subplot(gs[i, 0]); img(ax, p); ax.set_title(t, loc="left"); letter(ax, "ABCD"[i], dx=0.0)
    save(fig, name)
fig = plt.figure(figsize=(W, 5.2)); gs = GridSpec(1, 2, figure=fig, wspace=0.05)
for i, (p, t) in enumerate([("step2a_out/program_maps_all_samples.png", "Greenwald 2024"), ("step2a_ravi/program_maps_all_samples.png", "Ravi 2022")]):
    ax = fig.add_subplot(gs[0, i]); img(ax, p); ax.set_title(t, loc="left"); letter(ax, "AB"[i], dx=0.0)
save(fig, "Figure_S1")
stack(["niche_out/niche_profiles.png", "niche_out_h/niche_profiles.png"], ["Niches from program scores (no cohort correction)", "Niches from program scores with Harmony"], "Figure_S2", 5.6)
stack(["prog/progression_niche_score_maps.png"], ["GSE318238: one IDH-mutant astrocytoma, 2015 → 2017 → 2019"], "Figure_S3", 5.4)
stack(["icb_out/gse121810_arm_scores.png"], ["GSE121810: adjuvant (n = 10) vs neoadjuvant (n = 12) pembrolizumab"], "Figure_S4", 2.6)
