// 11a_fetch_tcga_cbioportal.js -- run in a browser tab open on https://www.cbioportal.org (DevTools console or automation).
// Fetches RNA-seq (RSEM) for the niche-signature + Gavish MP genes from TCGA LGG and GBM PanCancer Atlas 2018,
// computes per-sample signature scores (mean z of log2(RSEM+1), z across all tumours), adds single genes and clinical
// data, and leaves a CSV in window.__CSV (saved as 01_raw_data/Bulk_TCGA/tcga_glioma_niche_scores.csv).
// SETS = contents of 04_results/10_niche_signatures/GliomaNicheScore_signatures.json + Gavish MPs (MES, Hypoxia, AC, NPC, OPC)
//        + _single: ANGPTL4, SDC1-4, ITGA5, ITGB1, CD44, VEGFA, ADM, SPP1, CHI3L1
// Steps (as run on 2026-11-23):
//  1. POST /api/genes/fetch?geneIdType=HUGO_GENE_SYMBOL  -> entrez ids
//  2. POST /api/molecular-profiles/{study}_rna_seq_v2_mrna/molecular-data/fetch {sampleListId:{study}_rna_seq_v2_mrna, entrezGeneIds}
//  3. GET  /api/studies/{study}/clinical-data?clinicalDataType=SAMPLE|PATIENT
//  4. per gene z-score across samples; signature = mean z over available genes
// Result: 674 tumours (514 LGG, 160 GBM); 342/353 genes found (11 alias/non-coding symbols missing).
