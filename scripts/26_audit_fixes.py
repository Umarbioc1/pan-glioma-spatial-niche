"""26_audit_fixes.py -- corrections from the number audit (23 Nov 2026), applied to v2 (UK) and Neuro-Oncology (US) sources."""
import sys
FIX = [
 # (old fragment, new fragment) -- written for the US source; UK variants handled by a second list
 ("Ravi et al. contributed 16 GBM and 2 IDH-mutant sections plus three non-neoplastic cortex sections, which served as the copy-number reference",
  "Ravi et al. contributed 16 GBM and 2 IDH-mutant sections plus three non-neoplastic cortex sections, two of which served as the copy-number reference"),
 ("and estimated fractions agreed with independent marker scores for neurons (r = 0.77), oligodendrocytes (0.72), myeloid cells (0.53) and MES-like cells (0.50) (Figure 1B)",
  "and estimated fractions agreed with independent marker-program scores across all spots (Spearman ρ) for neurons (0.77), oligodendrocytes (0.72) and MES-like cells (0.50), and more weakly for myeloid cells (0.34) (Figure 1B)"),
 ("It occurred in 13 sections and in both cohorts.",
  "It contributed at least 100 spots in 13 sections from both cohorts (12 GBM and 1 IDH-mutant)."),
 ("(median AUC 0.61–0.77, 13–14 sections)", "(median AUC 0.60–0.77, 13–14 sections)"),
 ("together with VMP1, NAMPT, PDPN, CD44, CHI3L1 and TNC (14 genes in total;",
  "together with VMP1, NAMPT, PDPN, CD44, CHI3L1 and TNC (14 of the 32 up-regulated NCDEGs;"),
 ("followed by VMP1 and CAV1, although no gene passed transcriptome-wide FDR",
  "followed by CLIC4 and ABCA1, with VMP1 and CAV1 the next-ranked NCDEGs, although no gene passed transcriptome-wide FDR"),
 ("but rarely in AC-, NPC- or OPC-like cells (2–13%)", "but rarely in AC-, NPC- or OPC-like cells (0–13%)"),
 ("(HR 1.65, FDR 1 × 10⁻⁴; I² = 0.09)", "(HR 1.65, FDR 7 × 10⁻⁵; I² = 0.09)"),
 ("(HR 1.65, FDR 1×10⁻⁴; I² = 0.09)", "(HR 1.65, FDR 7×10⁻⁵; I² = 0.09)"),
 ("and both CGGA cohorts reproduced the ordering.",
  "and both CGGA cohorts reproduced the lowest scores in oligodendroglioma and the highest in GBM, with IDH-mutant astrocytoma and IDH-wildtype lower-grade glioma similar in CGGA-693."),
 ("We scored them in TCGA lower-grade glioma and GBM (674 tumors)",
  "We scored them in TCGA lower-grade glioma and GBM (674 tumors, 618 with WHO 2021 type and survival data)"),
 ("We scored them in TCGA lower-grade glioma and GBM (674 tumours)",
  "We scored them in TCGA lower-grade glioma and GBM (674 tumours, 618 with WHO 2021 type and survival data)"),
 ("grade increase (P = .22–.26;", "grade increase (all P ≥ .22;"),
 ("grade increase (p = 0.22–0.26;", "grade increase (all p ≥ 0.22;"),
 ("36 tumor sections (97,871 spots) entered niche analysis and three cortex sections served as the copy-number reference.",
  "36 tumor sections (97,871 spots) entered niche analysis; three cortex sections passed quality control, two of which served as the copy-number reference."),
 ("36 tumour sections (97,871 spots) entered niche analysis and three cortex sections served as the copy-number reference.",
  "36 tumour sections (97,871 spots) entered niche analysis; three cortex sections passed quality control, two of which served as the copy-number reference."),
 ("(B) Agreement (Pearson r) between deconvolved cell-type fractions and independent marker-gene scores",
  "(B) Agreement (Spearman ρ across all spots) between deconvolved cell-type fractions and independent marker-program scores"),
]
for path in sys.argv[1:]:
    s = open(path).read(); n = 0
    for a, b in FIX:
        for A, B in ((a, b), (a.replace("tumor", "tumour"), b.replace("tumor", "tumour")), (a.replace("Figure 1B", "Fig 1B"), b.replace("Figure 1B", "Fig 1B"))):
            if A in s: s = s.replace(A, B); n += 1
    open(path, "w").write(s); print(path, n, "fixes")
