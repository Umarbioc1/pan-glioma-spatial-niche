"""17 - Anatomical validation of niche signatures in the Ivy Glioblastoma Atlas (IvyGAP).

Input  : fpkm_table.csv, rows-genes.csv, columns-samples.csv (gene_expression_matrix_2014-11-25.zip)
         rna_seq_samples_details.csv (sample -> donor, structure)
         sig_out/bulk_scoring_sets.json (niche signatures, Gavish MPs, single genes)
Primary: 'Anatomic Structures RNA Seq' samples dissected by reference histology
         (LE, IT, CT, CTmvp, CTpan), 10 tumours.
Scoring: log2(FPKM+1), z per gene across primary samples, mean z per signature (as for bulk cohorts).
Stats  : OLS score ~ region (reference CT) + tumour fixed effect (within-tumour contrasts; a random-intercept
         mixed model was singular for most scores because between-tumour variance was ~0 after z-scoring);
         paired Wilcoxon on per-tumour region means (CTpan vs CT, CTpan vs CTmvp); BH across all tests.
         Specificity: N6 adjusted for the Gavish hypoxia MP (residual) by region.
"""
import json, sys, os
import numpy as np, pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
import warnings; warnings.filterwarnings("ignore")

IVY = sys.argv[1] if len(sys.argv) > 1 else "ivy"
SETS = sys.argv[2] if len(sys.argv) > 2 else "sig_out/bulk_scoring_sets.json"
OUT = sys.argv[3] if len(sys.argv) > 3 else "ivy_out"
os.makedirs(OUT, exist_ok=True)

genes = pd.read_csv(f"{IVY}/rows-genes.csv")
fpkm = pd.read_csv(f"{IVY}/fpkm_table.csv", index_col=0)
fpkm.columns = fpkm.columns.astype(int)
sym = genes.set_index("gene_id")["gene_symbol"]
fpkm.index = sym.reindex(fpkm.index).values
fpkm = fpkm[~fpkm.index.isna()].groupby(level=0).mean()

det = pd.read_csv(f"{IVY}/rna_seq_samples_details.csv")
prim = det[det.study_name == "Anatomic Structures RNA Seq"].copy()
prim = prim[prim.sample_id.isin(fpkm.columns)]
REG = {"LE": "Leading edge", "IT": "Infiltrating tumour", "CT": "Cellular tumour",
       "CTmvp": "Microvascular proliferation", "CTpan": "Pseudopalisading (peri-necrotic)"}
prim["region"] = prim.structure_parent_acronym
prim = prim[prim.region.isin(REG)]
print("primary samples", len(prim), "tumours", prim.tumor_name.nunique(), prim.region.value_counts().to_dict())

X = np.log2(fpkm[prim.sample_id] + 1)
X = X[(X > 1).sum(axis=1) >= 5]                       # expressed in >=5 samples
Z = X.sub(X.mean(1), axis=0).div(X.std(1).replace(0, np.nan), axis=0)

sets = json.load(open(SETS))
single = sets.pop("_single") + ["NRN1", "ARRDC3", "CAMK2N1"]
cov, scores = {}, {}
for k, g in sets.items():
    gg = [x for x in g if x in Z.index]; cov[k] = f"{len(gg)}/{len(g)}"
    scores[k] = Z.loc[gg].mean()
for g in single:
    if g in Z.index: scores[g] = Z.loc[g]
S = pd.DataFrame(scores)
S.index = S.index.astype(int)
df = prim.set_index("sample_id")[["tumor_name", "donor_id", "region"]].join(S)

# N6 beyond hypoxia: residual of N6 on the Gavish hypoxia MP
b = np.polyfit(df.MP_Hypoxia, df.HypoxicMES_N6, 1)
df["N6_resid_hypoxia"] = df.HypoxicMES_N6 - np.polyval(b, df.MP_Hypoxia)
df.to_csv(f"{OUT}/ivygap_sample_scores.csv")

order = ["LE", "IT", "CT", "CTmvp", "CTpan"]
score_cols = [c for c in S.columns] + ["N6_resid_hypoxia"]
from scipy.stats import wilcoxon
pt = df.groupby(["tumor_name", "region"])[score_cols].mean()
rows = []
for c in score_cols:
    d = df[["tumor_name", "region", c]].rename(columns={c: "y"}).copy()
    d["region"] = pd.Categorical(d.region, categories=["CT", "LE", "IT", "CTmvp", "CTpan"])
    m = smf.ols("y ~ C(region) + C(tumor_name)", d).fit()
    for r in ["LE", "IT", "CTmvp", "CTpan"]:
        t = f"C(region)[T.{r}]"
        rows.append(dict(score=c, test="OLS+tumour FE", contrast=f"{r} vs CT", beta=m.params[t], se=m.bse[t], p=m.pvalues[t], n_tumours=np.nan))
    w = pt[c].unstack()
    for a, bb in [("CTpan", "CT"), ("CTpan", "CTmvp")]:
        ww = w[[a, bb]].dropna()
        st = wilcoxon(ww[a], ww[bb])
        rows.append(dict(score=c, test="paired Wilcoxon (tumour means)", contrast=f"{a} vs {bb}", beta=float((ww[a]-ww[bb]).median()), se=np.nan, p=st.pvalue, n_tumours=len(ww)))
res = pd.DataFrame(rows)
res["fdr"] = multipletests(res.p, method="fdr_bh")[1]
means = df.groupby("region")[score_cols].mean().loc[order].T
means["top_region"] = means[order].idxmax(axis=1)
# within-tumour: fraction of tumours where CTpan is the highest region (tumours with CTpan sampled)
frac = {}
for c in score_cols:
    w = pt[c].unstack()
    w = w[w["CTpan"].notna()]
    frac[c] = f"{int((w.idxmax(axis=1) == 'CTpan').sum())}/{len(w)}"
means["CTpan_highest_in_tumours"] = pd.Series(frac)
res.to_csv(f"{OUT}/ivygap_mixedlm_region_contrasts.csv", index=False)
means.to_csv(f"{OUT}/ivygap_region_means.csv")
json.dump(dict(coverage=cov, n_samples=len(df), n_tumours=int(df.tumor_name.nunique()),
               n_per_region=df.region.value_counts().to_dict()), open(f"{OUT}/ivygap_summary.json", "w"), indent=1)
pd.set_option("display.width", 200)
print(means.round(2))
print(res[res.score.isin(["HypoxicMES_N6","Myeloid_N3","Neuronal_N1","Oligodendroglial_N10","MP_Hypoxia","MP_MES","ANGPTL4","N6_resid_hypoxia","CA12","CAV1","NDRG1","NAMPT","VMP1","PDPN"]) & res.contrast.isin(["CTpan vs CT","CTpan vs CTmvp","LE vs CT"])].to_string())
