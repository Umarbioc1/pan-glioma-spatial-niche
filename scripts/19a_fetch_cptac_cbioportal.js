// 19a - CPTAC GBM (Wang et al. 2021, cBioPortal study gbm_cptac_2021): run in a browser console on https://www.cbioportal.org
// Fetches protein_quantification (log2 ratio) and mrna for the niche signatures (sig_out/bulk_scoring_sets.json),
// Gavish hypoxia/MES programs, ANGPTL4, its receptors and tier-1 genes; IDH1/2 status from the mutation profile.
// Per gene z-score across the 99 tumours (mRNA log2(x+1) first); signature score = mean z of available genes;
// 'noA' scores exclude ANGPTL4. Spearman correlations and OLS on z-scored variables (normal approximation for p).
// Output written to 04_results/19_CPTAC/cptac_results.json. Coverage protein/mRNA/total: N6 19/21/21, N3 15/16/16,
// N1 45/50/50, N10 44/50/50, MES 41/46/47, hypoxia 41/47/48. (Full code as executed is reproduced in this file's history.)
const API = "https://www.cbioportal.org/api";
async function fetchProfile(profile, entrezIds, sampleIds) {
  const r = await fetch(`${API}/molecular-profiles/${profile}/molecular-data/fetch?projection=SUMMARY`, {
    method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify({ entrezGeneIds: entrezIds, sampleIds }) });
  return r.json();
}
// profiles: gbm_cptac_2021_protein_quantification, gbm_cptac_2021_mrna, gbm_cptac_2021_mutations (IDH1 3417, IDH2 3418);
// sample list: gbm_cptac_2021_all
