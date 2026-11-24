"""
13_bulk_cgga_validation.py -- Aim 5 replication in CGGA mRNAseq_693 and mRNAseq_325 (RSEM FPKM, cgga.org.cn, 2020-05-06 release)
Same signatures and scoring as TCGA (11a/11): log2(FPKM+1), per-gene z within cohort, signature = mean z.
WHO 2021 groups from IDH and 1p/19q status (+ grade): GBM IDH-wt (grade IV, IDH-wt); IDH-wt lower-grade (II/III, IDH-wt);
Astrocytoma IDH-mut (IDH-mut, non-codel, II-IV); Oligodendroglioma (IDH-mut, codel).
Main survival analyses use PRIMARY tumours; recurrent/secondary kept for a sensitivity analysis.
Then fixed-effect inverse-variance meta-analysis of log-HR across TCGA + CGGA693 + CGGA325.
usage: python3 13_bulk_cgga_validation.py <cgga_raw_dir> <signatures+MP json> <tcga cox_results.csv> <outdir>
"""
import sys, os, json, numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
from scipy.stats import kruskal, norm
raw, setf, tcga_cox, out = sys.argv[1:5]; os.makedirs(out, exist_ok=True)
SETS = json.load(open(setf)); SINGLE = SETS.pop("_single")
SIG = ["HypoxicMES_N6", "Myeloid_N3", "Neuronal_N1", "Oligodendroglial_N10"]

def load(c):
    E = pd.read_csv(f"{raw}/CGGA.mRNAseq_{c}.RSEM-genes.20200506.txt", sep="\t", index_col=0)
    E = np.log2(E.groupby(level=0).mean() + 1)
    txt = open(f"{raw}/CGGA.mRNAseq_{c}_clinical.20200506.txt", newline="").read().replace("\r\n", "\n").replace("\r", "\n")
    from io import StringIO
    C = pd.read_csv(StringIO(txt), sep="\t").set_index("CGGA_ID")
    C.columns = ["PRS_type", "Histology", "Grade", "Gender", "Age", "OS_days", "OS_event", "Radio", "TMZ", "IDH", "codel", "MGMT"]
    Z = E.sub(E.mean(1), axis=0).div(E.std(1).replace(0, np.nan), axis=0)
    S = pd.DataFrame({k: Z.reindex([g for g in v if g in Z.index]).mean() for k, v in SETS.items()})
    for g in SINGLE:
        if g in E.index: S[g] = E.loc[g]
    D = C.join(S, how="inner"); D["cohort"] = f"CGGA{c}"
    cov = {k: f"{sum(g in Z.index for g in v)}/{len(v)}" for k, v in SETS.items()}
    return D, cov

frames, covs = [], {}
for c in ["693", "325"]:
    D, cv = load(c); frames.append(D); covs[c] = cv
D = pd.concat(frames)
D["grade"] = D.Grade.map({"WHO II": 2, "WHO III": 3, "WHO IV": 4})
D["OS_MONTHS"] = D.OS_days / 30.44
idh, cod = D.IDH.astype(str), D.codel.astype(str)
D["who"] = np.select([(idh == "Wildtype") & (D.grade == 4), (idh == "Wildtype") & D.grade.isin([2, 3]),
                      (idh == "Mutant") & (cod == "Non-codel"), (idh == "Mutant") & (cod == "Codel")],
                     ["GBM, IDH-wt", "IDH-wt lower-grade", "Astrocytoma, IDH-mut", "Oligodendroglioma"], default="")
D["primary"] = (D.PRS_type == "Primary").astype(int)
D = D[(D.who != "") & D.OS_MONTHS.notna() & (D.OS_MONTHS > 0) & D.OS_event.notna()].copy()
D["OS_event"] = D.OS_event.astype(int)
D.to_csv(f"{out}/cgga_analysis_table.csv")
json.dump(covs, open(f"{out}/signature_gene_coverage.json", "w"), indent=1)
print(pd.crosstab([D.cohort, D.who], D.grade)); print(covs)

P = D[D.primary == 1]
order = ["Oligodendroglioma", "Astrocytoma, IDH-mut", "IDH-wt lower-grade", "GBM, IDH-wt"]
M = P.groupby(["cohort", "who"])[SIG].mean().round(3); M.to_csv(f"{out}/score_means_by_type.csv"); print(M)
A = P[P.who == "Astrocytoma, IDH-mut"]
print(A.groupby(["cohort", "grade"])[SIG].mean().round(3))
A.groupby(["cohort", "grade"])[SIG].mean().round(3).to_csv(f"{out}/astro_IDHmut_scores_by_grade.csv")

def cox(df, var, covs, label, cohort):
    x = df[["OS_MONTHS", "OS_event", var] + covs].dropna().copy()
    if x.OS_event.sum() < 10 or x[var].std() == 0: return None
    x[var] = (x[var] - x[var].mean()) / x[var].std()
    for c in covs:
        if x[c].nunique() > 2: x[c] = (x[c] - x[c].mean()) / x[c].std()
    x = x.loc[:, x.std() > 0] if len(covs) else x
    try: c = CoxPHFitter(penalizer=0.0).fit(x, "OS_MONTHS", "OS_event")
    except Exception: return None
    r = c.summary.loc[var]
    return {"cohort": cohort, "group": label, "variable": var, "adjusted_for": "+".join(covs) or "none", "n": len(x),
            "events": int(x.OS_event.sum()), "HR_per_SD": r["exp(coef)"], "CI_low": r["exp(coef) lower 95%"],
            "CI_high": r["exp(coef) upper 95%"], "p": r["p"], "logHR": r["coef"], "se": r["se(coef)"]}
for k in ["is_GBMwt", "is_astro", "is_IDHwtLGG"]: pass
D["is_GBMwt"] = (D.who == "GBM, IDH-wt").astype(int); D["is_astro"] = (D.who == "Astrocytoma, IDH-mut").astype(int)
D["is_IDHwtLGG"] = (D.who == "IDH-wt lower-grade").astype(int)
res = []
for coh in ["CGGA693", "CGGA325"]:
    for prim_only in [True, False]:
        X = D[(D.cohort == coh) & ((D.primary == 1) if prim_only else True)]
        extra = [] if prim_only else ["primary"]
        tag = "" if prim_only else " (incl. recurrent)"
        groups = {"All gliomas": (X, ["Age", "grade", "is_GBMwt", "is_astro", "is_IDHwtLGG"] + extra)}
        for w in order: groups[w] = (X[X.who == w], ["Age"] + ([] if w == "GBM, IDH-wt" else ["grade"]) + extra)
        for g, (df, cv) in groups.items():
            for v in SIG + ["ANGPTL4"]:
                r = cox(df, v, cv, g + tag, coh); r and res.append(r)
            for v in ["HypoxicMES_N6", "ANGPTL4"]:
                r = cox(df, v, cv + ["MP_Hypoxia"], g + tag, coh); r and res.append(r)
                r = cox(df, v, cv + ["MP_MES"], g + tag, coh); r and res.append(r)
R = pd.DataFrame(res); R.to_csv(f"{out}/cox_results_cgga.csv", index=False)

# meta-analysis with TCGA (primary-tumour models, same adjustment sets)
T = pd.read_csv(tcga_cox); T = T[T.endpoint == "OS"].copy(); T["cohort"] = "TCGA"
T["logHR"] = np.log(T.HR_per_SD); T["se"] = (np.log(T.CI_high) - np.log(T.CI_low)) / (2 * 1.959964)
T["adjusted_for"] = T.adjusted_for.str.replace("AGE", "Age")
R2 = R[~R.group.str.contains("incl.")]
def key(a): return a.str.replace("BUFFA_HYPOXIA_SCORE", "BUFFA").str.contains("MP_Hypoxia").map({True: "+MP_Hypoxia", False: ""}) + \
                   a.str.contains("MP_MES").map({True: "+MP_MES", False: ""})
T["adj"] = key(T.adjusted_for); R2 = R2.assign(adj=key(R2.adjusted_for))
T = T[~T.adjusted_for.str.contains("BUFFA")]
ALL = pd.concat([T[R2.columns.intersection(T.columns).tolist()], R2[R2.columns.intersection(T.columns).tolist()]])
meta = []
for (g, v, a), x in ALL.groupby(["group", "variable", "adj"]):
    if len(x) < 2: continue
    w = 1 / x.se ** 2; b = (w * x.logHR).sum() / w.sum(); se = np.sqrt(1 / w.sum())
    Q = (w * (x.logHR - b) ** 2).sum(); I2 = max(0, (Q - (len(x) - 1)) / Q) if Q > 0 else 0
    meta.append({"group": g, "variable": v, "extra_adjustment": a or "age (+grade, type)", "cohorts": ",".join(x.cohort),
                 "n_total": int(x.n.sum()), "events": int(x.events.sum()), "HR_meta": np.exp(b), "CI_low": np.exp(b - 1.96 * se),
                 "CI_high": np.exp(b + 1.96 * se), "p_meta": 2 * norm.sf(abs(b / se)), "I2": round(I2, 2),
                 "per_cohort_HR": "; ".join(f"{c}:{h:.2f}" for c, h in zip(x.cohort, np.exp(x.logHR)))})
MT = pd.DataFrame(meta).sort_values(["variable", "extra_adjustment", "group"]); MT.to_csv(f"{out}/meta_TCGA_CGGA.csv", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 60)
print(MT.round(4).to_string())

# figures: scores by type/grade (primary), KM within types (pooled CGGA, primary)
fig, ax = plt.subplots(1, 4, figsize=(16, 4.2))
cats = [("Oligo", P.who == "Oligodendroglioma"), ("Astro G2", (P.who == "Astrocytoma, IDH-mut") & (P.grade == 2)),
        ("Astro G3", (P.who == "Astrocytoma, IDH-mut") & (P.grade == 3)), ("Astro G4", (P.who == "Astrocytoma, IDH-mut") & (P.grade == 4)),
        ("IDH-wt LGG", P.who == "IDH-wt lower-grade"), ("GBM IDH-wt", P.who == "GBM, IDH-wt")]
for a, s in zip(ax, SIG):
    a.boxplot([P.loc[m, s] for _, m in cats], tick_labels=[f"{c}\n(n={m.sum()})" for c, m in cats], showfliers=False)
    a.set_title(f"{s} (CGGA primary)", fontsize=10); a.tick_params(axis="x", rotation=45, labelsize=8)
plt.tight_layout(); plt.savefig(f"{out}/cgga_scores_by_type_grade.png", dpi=200)
fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
for a, g in zip(ax, ["Astrocytoma, IDH-mut", "Oligodendroglioma", "GBM, IDH-wt"]):
    df = P[P.who == g]; hi = df.HypoxicMES_N6 > df.groupby("cohort").HypoxicMES_N6.transform("median")
    for lab, m, col in [("High", hi, "#B5172F"), ("Low", ~hi, "#2A8C82")]:
        KaplanMeierFitter().fit(df.OS_MONTHS[m], df.OS_event[m], label=f"{lab} (n={m.sum()})").plot_survival_function(ax=a, ci_show=False, color=col)
    p = logrank_test(df.OS_MONTHS[hi], df.OS_MONTHS[~hi], df.OS_event[hi], df.OS_event[~hi]).p_value
    a.set_title(f"CGGA {g}: hypoxic-MES (N6) score\nlog-rank p = {p:.2g} (median split within cohort)", fontsize=10)
    a.set_xlabel("Overall survival (months)"); a.set_ylabel("Survival probability")
plt.tight_layout(); plt.savefig(f"{out}/cgga_KM_HypoxicMES_N6.png", dpi=200)
