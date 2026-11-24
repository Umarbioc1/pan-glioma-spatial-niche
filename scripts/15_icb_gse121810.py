"""15_icb_gse121810.py -- niche signatures in recurrent GBM treated with pembrolizumab (Cloughesy 2019, GSE121810).
NCBI-generated TPM (GRCh38.p13) + NCBI annotation. GEO records only the arm: neoadjuvant (n=14, tumour resected
after one pre-surgical dose) vs adjuvant (n=15, resected before any anti-PD-1). Survival/response are not in GEO.
Scores: mean z of log2(TPM+1) per signature (z across the 29 tumours). Arm comparison: two-sided Mann-Whitney, BH FDR.
Positive controls from the original paper: IFN-gamma/T-cell (Ayers 2017 6-gene) up and cell-cycle down after neoadjuvant."""
import json, numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
src = "/mnt/user-data/uploads/GLIOMA SPATIAL PROJECT/01_raw_data"
T = pd.read_csv(f"{src}/GSE121810_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz", sep="\t", index_col=0)
A = pd.read_csv(f"{src}/Human.GRCh38.p13.annot.tsv.gz", sep="\t", index_col=0, usecols=["GeneID", "Symbol"])
T.index = A.reindex(T.index).Symbol.values; T = T[T.index.notna()].groupby(level=0).mean()
E = np.log2(T + 1); E = E[(T > 1).sum(1) >= 3]
Z = E.sub(E.mean(1), axis=0).div(E.std(1), axis=0)
S = json.load(open("sig_out/bulk_scoring_sets.json")); single = S.pop("_single")
S["IFNG_Tcell_control"] = ["IDO1", "CXCL10", "CXCL9", "HLA-DRA", "STAT1", "IFNG"]
S["CellCycle_control"] = ["MKI67", "TOP2A", "CCNB1", "CDK1", "BUB1", "AURKA", "CCNA2", "CDC20", "PLK1", "TPX2"]
sc = pd.DataFrame({k: Z.reindex([g for g in v if g in Z.index]).mean() for k, v in S.items()})
cov = {k: f"{sum(g in Z.index for g in v)}/{len(v)}" for k, v in S.items()}
for g in ["ANGPTL4", "CA12", "NAMPT", "SPP1", "CD44"]:
    if g in E.index: sc[g] = E.loc[g]
gsm = sc.index.str.extract(r"(\d+)$")[0].astype(int).values
sc["arm"] = np.where(gsm <= 3447021, "neoadjuvant", "adjuvant")
sc.to_csv("icb_out/gse121810_scores.csv")
rows = []
for c in [c for c in sc.columns if c != "arm"]:
    a, b = sc.loc[sc.arm == "neoadjuvant", c], sc.loc[sc.arm == "adjuvant", c]
    rows.append({"score": c, "mean_neoadjuvant": a.mean(), "mean_adjuvant": b.mean(), "diff": a.mean() - b.mean(),
                 "p_MWU": mannwhitneyu(a, b).pvalue})
R = pd.DataFrame(rows); R["fdr"] = multipletests(R.p_MWU, method="fdr_bh")[1]
R.round(4).to_csv("icb_out/gse121810_arm_comparison.csv", index=False)
print(sc.arm.value_counts().to_dict(), cov); print(R.round(4).to_string())
show = ["IFNG_Tcell_control", "CellCycle_control", "HypoxicMES_N6", "Myeloid_N3", "MP_MES", "MP_Hypoxia", "ANGPTL4", "Neuronal_N1"]
fig, ax = plt.subplots(1, len(show), figsize=(2.1 * len(show), 3.8))
for a, c in zip(ax, show):
    d = [sc.loc[sc.arm == x, c] for x in ["adjuvant", "neoadjuvant"]]
    a.boxplot(d, tick_labels=["adj.", "neo-adj."], showfliers=False)
    for i, v in enumerate(d): a.scatter(np.random.default_rng(0).normal(i + 1, 0.06, len(v)), v, s=10, c="#B5172F" if i else "#2A8C82", zorder=3)
    p = R.set_index("score").loc[c, "p_MWU"]; a.set_title(f"{c}\np = {p:.2g}", fontsize=8)
plt.suptitle("GSE121810 recurrent GBM: adjuvant (n=15) vs neoadjuvant (n=14) pembrolizumab", fontsize=10)
plt.tight_layout(); plt.savefig("icb_out/gse121810_arm_scores.png", dpi=200)
