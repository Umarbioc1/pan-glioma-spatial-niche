"""20 - Niche signatures through glioma recurrence in GLASS (Varn et al. 2022; Synapse syn17038081).

Inputs (Synapse exports, renamed in 01_raw_data/GLASS):
  gene_tpm_matrix_all_samples.tsv (syn69961520), clinical_surgeries (syn69931127), clinical_cases (syn69931132),
  analysis_rna_silver_set (syn69930927; one quality-checked RNA pair per patient: tumour a = earlier, b = later)
Scoring : log2(TPM+1), z per gene across all RNA samples, mean z per signature (as for TCGA/CGGA).
Tests   : paired change b - a within WHO subtype (IDH-wt; IDH-mut non-codeleted = astrocytoma; codeleted = oligodendroglioma),
          Wilcoxon signed-rank; within astrocytoma, change compared between pairs with and without grade increase (Mann-Whitney);
          hypoxia-adjusted N6 (residual after regression on the Gavish hypoxia program across samples);
          Cox models for overall survival (from diagnosis; clinical_cases) on the score of the earlier tumour,
          adjusted for age and subtype (lifelines). BH correction across the paired tests.
usage: python3 20_glass_longitudinal.py <glass_dir> <bulk_scoring_sets.json> <outdir>
"""
import sys, os, json, numpy as np, pandas as pd
from scipy.stats import wilcoxon, mannwhitneyu
from statsmodels.stats.multitest import multipletests
from lifelines import CoxPHFitter
import warnings; warnings.filterwarnings("ignore")

GD, SETS_F, OUT = sys.argv[1:4]
os.makedirs(OUT, exist_ok=True)
tpm = pd.read_csv(f"{GD}/gene_tpm_matrix_all_samples.tsv", sep="\t", index_col=0)
tpm.columns = [c.replace(".", "-") for c in tpm.columns]
tpm = tpm.groupby(level=0).mean()
X = np.log2(tpm + 1)
X = X[(X > 1).sum(1) >= 20]
Z = X.sub(X.mean(1), axis=0).div(X.std(1), axis=0)

sets = json.load(open(SETS_F)); single = sets.pop("_single")
keep = ["HypoxicMES_N6", "Myeloid_N3", "Neuronal_N1", "Oligodendroglial_N10", "MP_Hypoxia", "MP_MES"]
S, cov = {}, {}
for k in keep:
    g = [x for x in sets[k] if x in Z.index]; cov[k] = f"{len(g)}/{len(sets[k])}"; S[k] = Z.loc[g].mean()
for g in ["ANGPTL4", "VEGFA", "CA12", "CAV1", "NDRG1", "PDPN", "NAMPT"]:
    if g in Z.index: S[g] = Z.loc[g]
S = pd.DataFrame(S)
b = np.polyfit(S.MP_Hypoxia, S.HypoxicMES_N6, 1); S["N6_resid_hypoxia"] = S.HypoxicMES_N6 - np.polyval(b, S.MP_Hypoxia)
b = np.polyfit(S.MP_Hypoxia, S.ANGPTL4, 1); S["ANGPTL4_resid_hypoxia"] = S.ANGPTL4 - np.polyval(b, S.MP_Hypoxia)
S.index.name = "aliquot_barcode"; S["sample_barcode"] = [x[:15] for x in S.index]
S.to_csv(f"{OUT}/glass_sample_scores.csv")

su = pd.read_csv(f"{GD}/clinical_surgeries.csv"); ca = pd.read_csv(f"{GD}/clinical_cases.csv"); si = pd.read_csv(f"{GD}/analysis_rna_silver_set.csv")
su = su.dropna(subset=["sample_barcode"]).drop_duplicates("sample_barcode").set_index("sample_barcode")
G = {"II": 2, "III": 3, "IV": 4}
P = si[["case_barcode", "tumor_barcode_a", "tumor_barcode_b"]].copy()
P = P[P.tumor_barcode_a.isin(S.index) & P.tumor_barcode_b.isin(S.index)]
P["sa"], P["sb"] = P.tumor_barcode_a.str[:15], P.tumor_barcode_b.str[:15]
for t in ["a", "b"]:
    m = su.reindex(P[f"s{t}"])
    P[f"grade_{t}"] = m.grade.map(G).values; P[f"subtype_{t}"] = m.idh_codel_subtype.values
    P[f"interval_{t}"] = m.surgical_interval_mo.values
    P[f"tmz_{t}"] = m.treatment_tmz.values; P[f"rt_{t}"] = m.treatment_radiotherapy.values
P["subtype"] = P.subtype_a.fillna(P.subtype_b)
P["grade_change"] = np.where(P.grade_b > P.grade_a, "increase", np.where(P.grade_b.notna() & P.grade_a.notna(), "stable", "unknown"))
P["is_primary_a"] = P.tumor_barcode_a.str[13:15] == "TP"
scs = ["HypoxicMES_N6", "N6_resid_hypoxia", "Myeloid_N3", "Neuronal_N1", "Oligodendroglial_N10", "MP_Hypoxia", "MP_MES", "ANGPTL4", "ANGPTL4_resid_hypoxia"]
for c in scs:
    P[f"{c}_a"] = S.loc[P.tumor_barcode_a, c].values; P[f"{c}_b"] = S.loc[P.tumor_barcode_b, c].values
    P[f"{c}_delta"] = P[f"{c}_b"] - P[f"{c}_a"]
P.to_csv(f"{OUT}/glass_pairs_scores.csv", index=False)

rows = []
groups = [("IDHwt", P.subtype == "IDHwt"), ("IDHmut-noncodel", P.subtype == "IDHmut-noncodel"), ("IDHmut-codel", P.subtype == "IDHmut-codel"),
          ("IDHmut-noncodel, grade increase", (P.subtype == "IDHmut-noncodel") & (P.grade_change == "increase")),
          ("IDHmut-noncodel, grade stable", (P.subtype == "IDHmut-noncodel") & (P.grade_change == "stable"))]
for gname, mask in groups:
    d = P[mask]
    for c in scs:
        x = d[f"{c}_delta"].dropna()
        if len(x) < 6: continue
        rows.append(dict(group=gname, score=c, n_pairs=len(x), mean_a=d[f"{c}_a"].mean(), mean_b=d[f"{c}_b"].mean(),
                         median_delta=x.median(), frac_up=(x > 0).mean(), p=wilcoxon(x).pvalue))
R = pd.DataFrame(rows); R["fdr"] = multipletests(R.p, method="fdr_bh")[1]
R.to_csv(f"{OUT}/glass_paired_changes.csv", index=False)
# grade increase vs stable within astrocytoma
A = P[(P.subtype == "IDHmut-noncodel") & P.grade_change.isin(["increase", "stable"])]
gi = []
for c in scs:
    u, s_ = A.loc[A.grade_change == "increase", f"{c}_delta"].dropna(), A.loc[A.grade_change == "stable", f"{c}_delta"].dropna()
    if len(u) >= 5 and len(s_) >= 5:
        gi.append(dict(score=c, n_increase=len(u), n_stable=len(s_), median_delta_increase=u.median(), median_delta_stable=s_.median(),
                       p_MWU=mannwhitneyu(u, s_).pvalue))
GI = pd.DataFrame(gi); GI.to_csv(f"{OUT}/glass_astro_grade_increase_vs_stable.csv", index=False)
# survival: OS from diagnosis vs score of earlier tumour (primary only)
C = ca.set_index("case_barcode")[["case_age_diagnosis_years", "case_vital_status", "case_overall_survival_mo"]]
Q = P[P.is_primary_a].join(C, on="case_barcode").dropna(subset=["case_overall_survival_mo", "case_age_diagnosis_years"])
Q["event"] = (Q.case_vital_status == "dead").astype(int); Q = Q[Q.case_overall_survival_mo > 0]
cox = []
for sub in ["all", "IDHwt", "IDHmut-noncodel"]:
    q = Q if sub == "all" else Q[Q.subtype == sub]
    for c in ["HypoxicMES_N6", "N6_resid_hypoxia", "ANGPTL4", "Myeloid_N3", "Neuronal_N1"]:
        dd = pd.DataFrame({"T": q.case_overall_survival_mo, "E": q.event, "x": (q[f"{c}_a"] - q[f"{c}_a"].mean()) / q[f"{c}_a"].std(),
                           "age": q.case_age_diagnosis_years})
        if sub == "all":
            for s_ in ["IDHmut-noncodel", "IDHmut-codel"]: dd[s_] = (q.subtype == s_).astype(int).values
        dd = dd.dropna()
        try:
            f = CoxPHFitter().fit(dd, "T", "E"); r = f.summary.loc["x"]
            cox.append(dict(subset=sub, score=c, n=len(dd), events=int(dd.E.sum()), HR=r["exp(coef)"], lo=r["exp(coef) lower 95%"], hi=r["exp(coef) upper 95%"], p=r["p"]))
        except Exception as e: print("cox fail", sub, c, e)
CX = pd.DataFrame(cox); CX.to_csv(f"{OUT}/glass_cox_primary_scores.csv", index=False)
json.dump(dict(coverage=cov, n_rna_samples=int(S.shape[0]), n_pairs=int(len(P)), pairs_by_subtype=P.subtype.value_counts(dropna=False).to_dict(),
               astro_grade_change=A.grade_change.value_counts().to_dict(), pairs_primary_a=int(P.is_primary_a.sum())),
          open(f"{OUT}/glass_summary.json", "w"), indent=1)
pd.set_option("display.width", 250)
print(json.load(open(f"{OUT}/glass_summary.json")))
print(R.round(4).to_string()); print(GI.round(4).to_string()); print(CX.round(4).to_string())
