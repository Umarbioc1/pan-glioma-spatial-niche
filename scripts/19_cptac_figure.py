"""19 - Figure S8: protein-level check of niche signatures and ANGPTL4 in CPTAC GBM (99 tumours).
Input: cptac_out/cptac_results.json (computed in-browser from cBioPortal gbm_cptac_2021 by 19a_fetch_cptac_cbioportal.js)."""
import sys, json, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.linewidth": 0.6, "pdf.fonttype": 42})
R = json.load(open(sys.argv[1] if len(sys.argv) > 1 else "cptac_out/cptac_results.json")); FIG = sys.argv[2] if len(sys.argv) > 2 else "manuscript_figs"
RED, GREY, BLUE = "#c0392b", "#9a9a9a", "#2A6FB0"
fig, axs = plt.subplots(1, 3, figsize=(7.6, 2.9), gridspec_kw=dict(width_ratios=[1.3, 1.15, 1.1], wspace=0.75))
ax = axs[0]; pv = R["prot_vs_rna"]
keys = ["N6", "N3", "N1", "N10", "HYP", "MES", "ANGPTL4", "CAV1", "CA12", "NDRG1", "PDPN", "NAMPT", "VMP1", "NRN1"]
lab = {"N6": "N6 signature", "N3": "N3 signature", "N1": "N1 signature", "N10": "N10 signature", "HYP": "Hypoxia MP", "MES": "MES MP"}
v = [pv[k]["rho"] for k in keys]
ax.barh(range(len(keys)), v, color=[BLUE if k in lab else (RED if k == "ANGPTL4" else GREY) for k in keys], height=0.65)
ax.set_yticks(range(len(keys))); ax.set_yticklabels([lab.get(k, k) for k in keys]); ax.invert_yaxis(); ax.set_xlim(0, 1)
for i, x in enumerate(v): ax.text(x + 0.02, i, f"{x:.2f}", va="center", fontsize=6)
ax.set_xlabel("Spearman ρ, protein vs mRNA"); ax.set_title("Protein–mRNA agreement", loc="left", fontsize=7)
ax.text(-0.5, 1.03, "A", transform=ax.transAxes, fontsize=10, fontweight="bold")
ax = axs[1]; a = R["ANGPTL4p_vs"]; rc = R["ANGPTL4p_receptors"]
items = [("Hypoxia MP†", a["HYPnoA_p"]), ("N6 signature†", a["N6noA_p"]), ("N3 signature", a["N3_p"]), ("MES MP", a["MES_p"]), ("N1 signature", a["N1_p"]),
         ("ITGA5", rc["ITGA5"]), ("ITGB1", rc["ITGB1"]), ("CD44", rc["CD44"]), ("SDC2", rc["SDC2"]), ("SDC4", rc["SDC4"]), ("CDH11", rc["CDH11"])]
v = [x["rho"] for _, x in items]
ax.barh(range(len(items)), v, color=[BLUE] * 5 + [RED] * 6, height=0.65); ax.axvline(0, color="k", lw=0.5)
ax.set_yticks(range(len(items))); ax.set_yticklabels([n for n, _ in items]); ax.invert_yaxis(); ax.set_xlim(-0.35, 0.95)
for i, (n, x) in enumerate(items):
    ax.text(max(x["rho"], 0) + 0.02, i, ("***" if x["p"] < 1e-3 else "**" if x["p"] < 0.01 else "*" if x["p"] < 0.05 else ""), va="center", fontsize=6)
ax.axhline(4.5, color="0.8", lw=0.6)
ax.set_xlabel("Spearman ρ with ANGPTL4 protein"); ax.set_title("ANGPTL4 protein correlates", loc="left", fontsize=7)
ax.text(-0.62, 1.03, "B", transform=ax.transAxes, fontsize=10, fontweight="bold")
ax = axs[2]
mods = [("protein", "model_p_all", "p_"), ("protein\n+ MES", "model_p_mes", "p_"), ("mRNA", "model_r_all", "r_"), ("mRNA\n+ MES", "model_r_mes", "r_")]
x = np.arange(len(mods)); w = 0.36
for j, (var, col, lb) in enumerate([("N6noA", RED, "N6 signature†"), ("HYPnoA", BLUE, "Hypoxia MP†")]):
    b = [R[m][p + var]["beta"] for _, m, p in mods]; pp = [R[m][p + var]["p"] for _, m, p in mods]
    ax.bar(x + (j - 0.5) * w, b, w, color=col, label=lb)
    for k in range(len(mods)): ax.text(x[k] + (j - 0.5) * w, b[k] + 0.015, "*" if pp[k] < 0.05 else "ns", ha="center", fontsize=5.5)
ax.set_xticks(x); ax.set_xticklabels(["prot.", "prot.\n+MES", "mRNA", "mRNA\n+MES"], fontsize=6); ax.set_ylabel("standardised β for ANGPTL4")
ax.legend(frameon=False, fontsize=5.5, loc="upper left"); ax.set_ylim(0, 0.95); ax.set_title("Model for ANGPTL4", loc="left", fontsize=7)
ax.text(-0.38, 1.03, "C", transform=ax.transAxes, fontsize=10, fontweight="bold")
fig.savefig(f"{FIG}/Figure_S8.png", dpi=300, bbox_inches="tight"); fig.savefig(f"{FIG}/Figure_S8.pdf", bbox_inches="tight"); print("ok")
