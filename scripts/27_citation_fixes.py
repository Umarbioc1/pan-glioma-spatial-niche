"""27_citation_fixes.py -- reword sentences whose cited sources support a narrower claim (citation audit, 23 Nov 2026)."""
import sys, re
for path in sys.argv[1:]:
    s = open(path).read(); n0 = s
    uk = 'organise' in s
    # Intro: Ravi = spatial segregation of states; Greenwald = hypoxia-organised layers
    s = re.sub(r"In GBM, hypoxia and necrosis organi[sz]e malignant states into layered structures, (?:with MES-like cells close to|and MES-like cells lie close to) hypoxic regions and macrophages",
               ("In GBM, malignant states are spatially segregated and organised by hypoxia into layered structures, with MES-like cells close to hypoxic regions and macrophages" if uk else
                "In GBM, malignant states are spatially segregated and organized by hypoxia into layered structures, with MES-like cells close to hypoxic regions and macrophages"), s)
    # Intro: Karimi = gliomas and brain metastases
    s = s.replace("large spatial studies now span glioma types, treatment and recurrence", "large spatial studies now span glioma types, brain metastases, treatment and recurrence")
    # Discussion: Hara = macrophages drive MES-like; Neftel = MES-like states include a hypoxia-associated subprogram
    s = re.sub(r"Hypoxia and macrophage signals push GBM cells towards? MES-like states ?(\[\[?)3, ?28(\]\]?)",
               lambda m: ("Macrophage signals push GBM cells toward MES-like states" if not uk else "Macrophage signals push GBM cells towards MES-like states") +
               (" " if " [" in m.group(0) else "") + m.group(1) + "28" + m.group(2) + ", a state that includes a hypoxia-associated subprogram" + (" " if " [" in m.group(0) else "") + m.group(1) + "3" + m.group(2), s)
    # Discussion: Takano/Chedeville = HIF/hypoxia induction in glioma cells
    s = s.replace("ANGPTL4 is a known hypoxia-inducible gene in glioma cells", "ANGPTL4 is induced by hypoxia and HIF signalling in glioma cells" if uk else "ANGPTL4 is induced by hypoxia and HIF signaling in glioma cells")
    # Discussion: Garner = STAT3 target, mesenchymal xenografts, STAT3/ANGPTL4 and survival
    s = s.replace("and it is co-regulated with STAT3 in mesenchymal-subtype GBM, where high expression is associated with poor survival",
                  "and it is a STAT3 target that marks mesenchymal-subtype GBM xenografts, with high STAT3 and ANGPTL4 expression associated with poor survival")
    # Discussion: Varn did not report our ANGPTL4 result -> cite GLASS as the data source
    s = re.sub(r"Within 214 GLASS patients, niche-6 scores changed little at recurrence; ANGPTL4 rose nominally in astrocytomas that recurred at a higher grade ?(\[\[?\d+\]\]?|\[\d+\])",
               lambda m: "Within 214 GLASS patients" + (" " if m.group(1).startswith("[") and not m.group(1).startswith("[[") else "") + m.group(1) + ", niche-6 scores changed little at recurrence; ANGPTL4 rose nominally in astrocytomas that recurred at a higher grade", s)
    open(path, "w").write(s); print(path, "changed" if s != n0 else "UNCHANGED")
