"""22_renumber_supp_tables.py -- reorder/rename Supplementary_Tables_S1-S16.xlsx to the v2 manuscript numbering."""
import sys, openpyxl
src, dst = sys.argv[1], sys.argv[2]
NEW = [  # (old sheet, new sheet, description)
 ("S1_sections_QC","S1_sections_QC","Visium sections, QC metrics, WHO 2021 type and 1p/19q status (Greenwald 2024, Ravi 2022)"),
 ("S2_niche_frequency","S2_niche_frequency","Fraction of spots in each niche per section"),
 ("S3_LR_niche_enriched","S3_LR_niche_enriched","43 niche-enriched CellChatDB pairs: meta-analysis statistics, likely sender/receiver, GBM vs IDH-mut AUC"),
 ("S4_NCDEG_niche6_malignant","S4_NCDEG_niche6_malignant","Niche-6 malignant NCDEG effects before and after malignant-state adjustment"),
 ("S5_NCDEG_all","S5_NCDEG_all","All significant NCDEGs (malignant and myeloid focal groups)"),
 ("S13_NCDEG_hypoxia_adj","S6_NCDEG_hypoxia_adj","Niche-6 malignant effects for all tested genes: unadjusted, hypoxia-adjusted, hypoxia + state-adjusted"),
 ("S12_IvyGAP","S7_IvyGAP","IvyGAP region means and within-tumour region contrasts (LE, IT, CT, MVP, PAN)"),
 ("S12b_IvyGAP_samples","S7b_IvyGAP_samples","IvyGAP per-sample signature and gene scores"),
 ("S14_GBmap_expression","S8_GBmap_expression","GBmap single-cell expression of ANGPTL4, receptors and tier-1 targets by state; hypoxic vs non-hypoxic MES-like"),
 ("S15_CPTAC_protein","S9_CPTAC_protein","CPTAC GBM protein-level analysis (99 tumours): protein-mRNA agreement, ANGPTL4 correlations, models"),
 ("S6_niche_signatures","S10_niche_signatures","Niche gene signatures used for bulk scoring (N1, N3, N6, N10)"),
 ("S7_survival_meta","S11_survival_meta","Fixed-effect meta-analysis of per-SD Cox HRs (TCGA, CGGA-693, CGGA-325) with multiple-testing correction"),
 ("S8_cox_TCGA","S12_cox_TCGA","Per-cohort Cox models, TCGA (OS and PFS)"),
 ("S9_cox_CGGA","S13_cox_CGGA","Per-cohort Cox models, CGGA (primary and including recurrent)"),
 ("S16_GLASS_longitudinal","S14_GLASS_longitudinal","GLASS paired analysis: changes by subtype and grade change, grade-increase comparison, Cox"),
 ("S16b_GLASS_pairs","S14b_GLASS_pairs","GLASS per-pair scores for earlier (a) and later (b) tumours"),
 ("S11_ICB_GSE121810","S15_ICB_GSE121810","Niche and control scores by arm in GSE121810 (exploratory)"),
 ("S10_targets","S16_targets","48 niche-6 candidate genes with spatial evidence, Open Targets and DepMap annotations, tiers and drug notes"),
]
wb = openpyxl.load_workbook(src)
assert set(o for o,_,_ in NEW) | {"README"} == set(wb.sheetnames), wb.sheetnames
for o, n, _ in NEW: wb[o].title = "tmp_" + n
for o, n, _ in NEW: wb["tmp_" + n].title = n
wb._sheets = [wb["README"]] + [wb[n] for _, n, _ in NEW]
rd = wb["README"]
rd.delete_rows(2, rd.max_row)
for i, (_, n, d) in enumerate(NEW, start=2):
    rd.cell(i, 1, n); rd.cell(i, 2, d)
wb.active = 0
wb.save(dst); print("saved", dst, wb.sheetnames)
