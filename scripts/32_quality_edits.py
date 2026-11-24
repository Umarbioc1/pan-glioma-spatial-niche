"""32_quality_edits.py -- title, authors (as in the CREPT ASO manuscript) and language/structure edits for the Neuro-Oncology version.
Usage: python3 32_quality_edits.py ms_no_src.md ms_no_tail.md"""
import sys
src, tail = sys.argv[1:3]
TITLE = "Separating spatial niche from cell state and hypoxia nominates ANGPTL4 as a peri-necrotic ligand in glioblastoma"
E = [
 # ---- title page
 ("| Title | Pan-glioma spatial transcriptomics separates niche from cell state and hypoxia and nominates ANGPTL4 in the glioblastoma peri-necrotic niche |", f"| Title | {TITLE} |"),
 ("| Authors | Umar Farooq¹, **[co-authors]** |", "| Authors | Umar Farooq¹, Wei Deng¹, Jing Zhu¹, Jun Li¹, Xiaotong Yang¹,\\* and Zhijie Chang²,\\* |"),
 ("| Affiliations | ¹ **[Department, institution, Chongqing, China]** |",
  "| Affiliations | ¹ Jinfeng Laboratory, No. 313 Jinyue Road, High-Tech Zone, Chongqing 401329, China; ² State Key Laboratory of Membrane Biology, School of Basic Medical Sciences, School of Medicine, Tsinghua University, Beijing 100084, China |"),
 ("| Corresponding author | Umar Farooq, **[postal address, email, telephone]** |",
  "| Corresponding authors | \\*Xiaotong Yang, Jinfeng Laboratory, No. 313 Jinyue Road, High-Tech Zone, Chongqing 401329, China (yangxiaotong@jflab.ac.cn); \\*Zhijie Chang, School of Basic Medical Sciences, Tsinghua University, Beijing 100084, China (zhijiec@tsinghua.edu.cn) **[add telephone]** |"),
 ("| Funding | **[grant names and numbers, or \"None\"]** |", "| Funding | Chongqing Municipal Science and Technology Project (CSTB2025NSCQBSX0026, to U.F.); Jinfeng Laboratory (JFLKYXM202303AZ-203, to Z.C.) |"),
 ("| Conflict of interest | The authors declare no conflicts of interest. **[confirm]** |", "| Conflict of interest | The authors declare no conflicts of interest. |"),
 ("| Authorship | **[CRediT roles per author]** |", "| Authorship | See Authorship statement |"),
 # ---- abstract
 ("**Background.** Diffuse gliomas differ by WHO type, but it is unclear how their cellular neighborhoods differ and which niche signals are independent of the malignant cell state and hypoxia that define them.",
  "**Background.** Spatial transcriptomics increasingly nominates niche-specific targets in glioma, but in glioblastoma (GBM) niche, mesenchymal-like (MES-like) cell state and hypoxia are entangled. We asked how cellular neighborhoods differ across WHO 2021 diffuse glioma types and which niche signals persist once state and hypoxia are accounted for."),
 ("We analyzed 40 public Visium sections of IDH-wildtype glioblastoma (GBM), IDH-mutant astrocytoma and oligodendroglioma (109,415 spots)",
  "We analyzed 40 public Visium sections of IDH-wildtype GBM, IDH-mutant astrocytoma and oligodendroglioma (109,415 spots)"),
 ("Ten niches emerged. A hypoxic mesenchymal (MES) niche of MES-like tumor cells and macrophages was largely restricted to GBM (10.4% vs 1.2% of spots) and was the densest signaling hub.",
  "Ten niches emerged. A hypoxic MES niche of MES-like tumor cells and macrophages was largely restricted to GBM (10.4% vs 1.2% of spots) and was the densest signaling hub."),
 ("**Conclusions.** Separating niche from cell state and hypoxia identifies ANGPTL4 as a candidate peri-necrotic ligand acting on the adjacent vasculature in GBM.",
  "**Conclusions.** Most apparent niche genes in GBM mark cell state or hypoxia. Testing for this separates them from niche-specific signals and nominates ANGPTL4 as a peri-necrotic ligand that may act on the adjacent vasculature."),
 # ---- importance
 ("Using only public data, we built a 40-section spatial atlas across WHO 2021 diffuse glioma types and explicitly separated niche effects from malignant-state and hypoxia effects.",
  "Using only public data, we built a 40-section spatial atlas across WHO 2021 diffuse glioma types and tested every niche gene against malignant state and hypoxia."),
 # ---- introduction
 ("Two problems limit the use of these findings. Most analyses focus on one tumor type or use data that are not openly available, which hampers comparison across WHO types. And because hypoxia and the MES-like state co-occur, it is rarely tested whether a gene enriched in a niche reflects the niche itself or only the state or hypoxia that defines it, a distinction that matters for target nomination.",
  "Two gaps limit the use of these findings. First, most analyses cover one tumor type or rely on data that are not openly available, which hampers comparison across WHO types. Second, because hypoxia and the MES-like state co-occur, it is rarely tested whether a gene enriched in a niche reflects the niche itself or only the state or hypoxia that defines it. This distinction matters for target nomination: a gene that merely marks MES-like cells or hypoxic tissue is unlikely to be a niche-specific vulnerability."),
 ("We applied it to diffuse glioma using only public data, identified a GBM-restricted hypoxic MES niche, separated niche effects from malignant-state and hypoxia effects, and tested the leading niche ligand, ANGPTL4, in independent resources covering tumor anatomy, single cells, protein, outcome and recurrence.",
  "Here we applied this framework to diffuse glioma using only public data. We asked three questions: which niches distinguish GBM from IDH-mutant gliomas, which niche genes remain after accounting for malignant state and hypoxia, and whether the leading candidate is supported by tumor anatomy, single-cell, protein, outcome and longitudinal data."),
 # ---- results
 ("Niche 6 combined MES-like malignant cells with blood-derived macrophages and was almost restricted to GBM:",
  "Niche 6 combined MES-like malignant cells with blood-derived macrophages and was largely restricted to GBM:"),
 ("To find genes whose expression in a cell type depends on its niche (niche-dependent DEGs, NCDEGs), we used a Python implementation in the spirit of C-SIDE[[17]]. In spots rich in the focal cell type, expression was regressed on niche membership with adjustment for other cell-type fractions and depth; each section contributed one effect, and significance was tested across sections because neighboring spots are not independent.",
  "We next asked which genes change with niche within a given cell type (niche-dependent genes, NCDEGs). Following C-SIDE[[17]], we compared niche and non-niche spots rich in that cell type within each section, adjusting for admixture and depth, and tested effects across sections because neighboring spots are not independent."),
 # ---- discussion
 ("This niche is almost confined to GBM and is the densest signaling hub in the atlas. Separating niche from cell state and hypoxia singles out one ligand, ANGPTL4, whose source, location and receptors we could trace across five independent resources. Four points deserve comment.",
  "This niche is largely confined to GBM and is the densest signaling hub in the atlas. Most of its genes, however, mark the MES-like state or hypoxia rather than the niche itself. Separating these effects singles out one ligand, ANGPTL4, whose source, location and receptors we could trace across five independent resources. Four points deserve comment."),
 ("after accounting for state mix, most niche-6 genes lose enrichment, as do established targets such as VEGFA, CD44 and TNC, and canonical hypoxia genes such as VEGFA, ADM, PGK1 and LDHA lose it after adjustment for hypoxia.",
  "after accounting for state mix, most niche-6 genes lose enrichment, including established targets such as VEGFA, CD44 and TNC, and canonical hypoxia genes such as ADM, PGK1 and LDHA lose it after adjustment for hypoxia."),
 ("We found no earlier report placing this axis in a specific glioma niche. Because the evidence is computational,",
  "To our knowledge, this axis has not previously been placed in a specific glioma niche. The model makes testable predictions: ANGPTL4 protein should localize to pseudopalisading tumor cells next to integrin α5-positive vessels, and blocking ANGPTL4 or integrin α5β1 should reduce peri-necrotic vascular leakiness. Because the evidence is computational,"),
 ("In conclusion, an open, multi-cohort spatial analysis across diffuse glioma types identifies a GBM-restricted hypoxic MES niche, distinguishes niche-specific from state- and hypoxia-driven signals, and nominates ANGPTL4 as a candidate niche ligand produced around necrosis next to receptor-bearing vessels. The niche signature rises with WHO type and grade and stratifies survival at diagnosis.",
  "In conclusion, an open, multi-cohort spatial analysis across diffuse glioma types identifies a GBM-restricted hypoxic MES niche and shows that most of its genes reflect cell state or hypoxia. Testing for this distinguishes niche-specific signals and nominates ANGPTL4 as a candidate ligand produced around necrosis next to receptor-bearing vessels. The niche signature rises with WHO type and grade and stratifies survival at diagnosis."),
]
T = [
 ("## Funding\n\n**[Grant names and numbers, or \"This work received no specific funding.\"]**",
  "## Funding\n\nThis work was supported by the Chongqing Municipal Science and Technology Project (Grant No. CSTB2025NSCQBSX0026) awarded to U.F. and by Jinfeng Laboratory (Grant No. JFLKYXM202303AZ-203) awarded to Z.C."),
 ("The authors declare no conflicts of interest. **[Confirm for every author.]**", "The authors declare no conflicts of interest."),
 ("**[CRediT roles, e.g.: U.F.: conceptualization, data curation, formal analysis, visualization, writing (original draft). Add co-authors.]**",
  "Conceptualization and project design, U.F. and Z.C.; data collection, formal analysis and visualization, U.F.; writing (original draft preparation), U.F.; investigation and data curation, W.D. and J.Z.; writing (review and editing), J.L., X.Y. and Z.C.; supervision, X.Y. and Z.C.; project administration, X.Y.; funding acquisition, U.F. and Z.C. All authors have read and agreed to the submitted version of the manuscript."),
 ("We thank the authors of the original studies, the Ivy Glioblastoma Atlas Project, CPTAC and the GLASS Consortium for making their data publicly available.",
  "We thank Jinfeng Laboratory for providing research facilities and technical support, and the authors of the original studies, the Ivy Glioblastoma Atlas Project, CPTAC and the GLASS Consortium for making their data publicly available."),
]
for path, edits in ((src, E), (tail, T)):
    s = open(path).read(); miss = []
    for a, b in edits:
        if a in s: s = s.replace(a, b)
        else: miss.append(a[:70])
    open(path, "w").write(s); print(path, len(edits) - len(miss), "applied;", "missing:", miss)
