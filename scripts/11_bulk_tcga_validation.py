"""
11_bulk_tcga_validation.py -- Aim 5: niche signatures in TCGA LGG+GBM (PanCancer Atlas, cBioPortal API)
Input: tcga_glioma_niche_scores.csv (per-sample signature scores = mean z of log2(RSEM+1), z across all 674 tumours;
       single genes log2(RSEM+1); clinical from cBioPortal). Fetched via 03_scripts/11a_fetch_tcga_cbioportal.js
WHO 2021 grouping from PanCancer SUBTYPE:
  GBM_IDHwt -> GBM, IDH-wt ; LGG_IDHwt -> IDH-wt lower-grade (molecular GBM-like) ; *_IDHmut-non-codel -> Astrocytoma, IDH-mut ;
  LGG_IDHmut-codel -> Oligodendroglioma
Analyses: (1) scores by type/grade, (2) Cox OS per SD within type, adjusted for age (+grade),
(3) independence of HypoxicMES_N6 from Buffa hypoxia score and MES program, (4) ANGPTL4 alone, (5) KM median split.
"""
import sys, os, numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
from scipy.stats import kruskal, spearmanr
inp, out = sys.argv[1:3]; os.makedirs(out, exist_ok=True)
D = pd.read_csv(inp)
D = D[D.SAMPLE_TYPE.fillna("Primary") == "Primary"].copy()
D["OS_event"] = D.OS_STATUS.astype(str).str.startswith("1").astype(int)
D["PFS_event"] = D.PFS_STATUS.astype(str).str.startswith("1").astype(int)
m = {"GBM_IDHwt": "GBM, IDH-wt", "LGG_IDHwt": "IDH-wt lower-grade", "LGG_IDHmut-non-codel": "Astrocytoma, IDH-mut",
     "GBM_IDHmut-non-codel": "Astrocytoma, IDH-mut", "LGG_IDHmut-codel": "Oligodendroglioma"}
D["who"] = D.SUBTYPE.map(m)
D["grade"] = D.GRADE.map({"G2": 2, "G3": 3}); D.loc[D.study.str.startswith("gbm"), "grade"] = 4
D = D[D.who.notna() & D.OS_MONTHS.notna() & (D.OS_MONTHS > 0)].copy()
SIG = ["HypoxicMES_N6", "Myeloid_N3", "Neuronal_N1", "Oligodendroglial_N10"]
REF = ["MP_MES", "MP_Hypoxia", "MP_AC", "MP_NPC", "MP_OPC"]
D.to_csv(f"{out}/tcga_analysis_table.csv", index=False)
print(D.who.value_counts()); print(pd.crosstab(D.who, D.grade))

# (1) scores by type and grade
order = ["Oligodendroglioma", "Astrocytoma, IDH-mut", "IDH-wt lower-grade", "GBM, IDH-wt"]
S = D.groupby("who")[SIG + REF].mean().loc[order].round(3); S.to_csv(f"{out}/score_means_by_type.csv")
kw = {s: kruskal(*[D.loc[D.who == w, s] for w in order]).pvalue for s in SIG}
A = D[D.who == "Astrocytoma, IDH-mut"]
gr = {s: spearmanr(A.grade, A[s]) for s in SIG}
pd.DataFrame({"kruskal_p_type": kw, "astro_grade_rho": {k: v[0] for k, v in gr.items()},
              "astro_grade_p": {k: v[1] for k, v in gr.items()}}).to_csv(f"{out}/score_type_grade_tests.csv")
GS = A.groupby("grade")[SIG].mean().round(3); GS.to_csv(f"{out}/astro_IDHmut_scores_by_grade.csv")
print(S[SIG]); print(kw); print(GS)

# (2)-(4) Cox
def cox(df, var, covs, label, endpoint="OS"):
    x = df[[f"{endpoint}_MONTHS", f"{endpoint}_event", var] + covs].dropna().copy()
    if x[f"{endpoint}_event"].sum() < 10: return None
    x[var] = (x[var] - x[var].mean()) / x[var].std()
    for c in covs:
        if x[c].nunique() > 2: x[c] = (x[c] - x[c].mean()) / x[c].std()
    c = CoxPHFitter(penalizer=0.0).fit(x, f"{endpoint}_MONTHS", f"{endpoint}_event")
    r = c.summary.loc[var]
    return {"group": label, "variable": var, "endpoint": endpoint, "adjusted_for": "+".join(covs) or "none",
            "n": len(x), "events": int(x[f"{endpoint}_event"].sum()), "HR_per_SD": r["exp(coef)"],
            "CI_low": r["exp(coef) lower 95%"], "CI_high": r["exp(coef) upper 95%"], "p": r["p"]}
res = []
D["is_GBMwt"] = (D.who == "GBM, IDH-wt").astype(int); D["is_astro"] = (D.who == "Astrocytoma, IDH-mut").astype(int)
D["is_IDHwtLGG"] = (D.who == "IDH-wt lower-grade").astype(int)
groups = {"All gliomas": (D, ["AGE", "grade", "is_GBMwt", "is_astro", "is_IDHwtLGG"]),
          "GBM, IDH-wt": (D[D.who == "GBM, IDH-wt"], ["AGE"]),
          "Astrocytoma, IDH-mut": (A, ["AGE", "grade"]),
          "Oligodendroglioma": (D[D.who == "Oligodendroglioma"], ["AGE", "grade"]),
          "IDH-wt lower-grade": (D[D.who == "IDH-wt lower-grade"], ["AGE", "grade"])}
for g, (df, cv) in groups.items():
    for v in SIG + ["ANGPTL4"]:
        for ep in ["OS", "PFS"]:
            r = cox(df, v, cv, g, ep); r and res.append(r)
    # independence of HypoxicMES_N6 from established hypoxia / MES scores
    for extra in [["BUFFA_HYPOXIA_SCORE"], ["MP_MES"], ["MP_Hypoxia"]]:
        r = cox(df, "HypoxicMES_N6", cv + extra, g); r and res.append(r)
        r = cox(df, "ANGPTL4", cv + extra, g); r and res.append(r)
R = pd.DataFrame(res); R.to_csv(f"{out}/cox_results.csv", index=False)
pd.set_option("display.width", 250)
print(R[(R.endpoint == "OS")].round(4).to_string())

# (5) KM: HypoxicMES_N6 median split
fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
for a, (g, df) in zip(ax, [("All gliomas", D), ("Astrocytoma, IDH-mut", A), ("GBM, IDH-wt", D[D.who == "GBM, IDH-wt"])]):
    hi = df.HypoxicMES_N6 > df.HypoxicMES_N6.median()
    for lab, mm, col in [("High", hi, "#B5172F"), ("Low", ~hi, "#2A8C82")]:
        KaplanMeierFitter().fit(df.OS_MONTHS[mm], df.OS_event[mm], label=f"{lab} (n={mm.sum()})").plot_survival_function(ax=a, ci_show=False, color=col)
    p = logrank_test(df.OS_MONTHS[hi], df.OS_MONTHS[~hi], df.OS_event[hi], df.OS_event[~hi]).p_value
    a.set_title(f"{g}: hypoxic-MES niche score\nlog-rank p = {p:.2g}", fontsize=10); a.set_xlabel("Overall survival (months)"); a.set_ylabel("Survival probability")
plt.tight_layout(); plt.savefig(f"{out}/KM_HypoxicMES_N6.png", dpi=200)
# score by type/grade plot
fig, ax = plt.subplots(1, 4, figsize=(16, 4.2))
cats = [("Oligo", D.who == "Oligodendroglioma"), ("Astro G2", (A.grade == 2).reindex(D.index, fill_value=False)),
        ("Astro G3", (A.grade == 3).reindex(D.index, fill_value=False)), ("Astro G4", (A.grade == 4).reindex(D.index, fill_value=False)),
        ("IDH-wt LGG", D.who == "IDH-wt lower-grade"), ("GBM IDH-wt", D.who == "GBM, IDH-wt")]
for a, s in zip(ax, SIG):
    a.boxplot([D.loc[mm, s] for _, mm in cats], tick_labels=[c for c, _ in cats], showfliers=False)
    a.set_title(s, fontsize=10); a.tick_params(axis="x", rotation=45)
plt.tight_layout(); plt.savefig(f"{out}/scores_by_type_grade.png", dpi=200)
