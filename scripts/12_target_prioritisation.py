"""
12_target_prioritisation.py -- Pillar 4: niche-resolved targets
Candidates: niche-6 malignant NCDEGs (up, FDR<0.05) + niche-6 significant ligand/receptor genes (CellChatDB).
Open Targets Platform GraphQL (queried via the Open Targets connector, 2026-11-23; raw output ot_raw.json):
  target(ensemblId){approvedSymbol isEssential tractability{label modality value}
    drugAndClinicalCandidates{rows{drug{name maximumClinicalStage drugType}}}
    associatedDiseases(Bs:["MONDO_0018177"]){rows{score}}            # glioblastoma
    depMapEssentiality{screens{geneEffect diseaseFromSource}}}       # median over glioma lines
Tiers: 1 = spatial (NCDEG or LR) & state-adjusted beta>0 with p<0.05 & tractable (SM: approved..high-quality ligand;
          AB: approved..GO CC high conf); 1b = as 1 without tractability; 2 = spatial & tractable; 3 = spatial only.
Drug notes flag irrelevant listings (ROS1-fusion drugs on CD74/SDC4; matrix enzymes on collagens/FN1/LAMA5).
The table-building code is the inline block run in session; outputs: 04_results/12_targets/niche6_target_prioritisation.csv
"""
