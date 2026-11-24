# Pan-glioma spatial transcriptomics: separating niche, cell state and hypoxia

Code and processed results for:

> Farooq U, Deng W, Zhu J, Li J, Yang X, Chang Z. *Disentangling spatial niche, cell state and hypoxia across diffuse gliomas identifies ANGPTL4 as a peri-necrotic ligand in glioblastoma.* (submitted)

The study builds a 40-section public 10x Visium atlas spanning IDH-wildtype glioblastoma,
IDH-mutant astrocytoma and oligodendroglioma (109,415 spots), defines ten multicellular
niches, and tests every niche-enriched gene against malignant cell state and hypoxia.
ANGPTL4 is the gene least explained by either, and is validated in IvyGAP, GBmap, CPTAC,
TCGA, CGGA and GLASS.

**No new samples were generated. Every input dataset is public.**

## Input data (not redistributed here)

| Dataset | Source |
| --- | --- |
| Visium, 19 sections | Greenwald et al. 2024, Zenodo 12624860 |
| Visium, 28 sections | Ravi et al. 2022, Dryad |
| GBmap single-cell reference | CELLxGENE |
| Malignant meta-programs | Gavish et al. 2023 |
| Ligand–receptor pairs | CellChatDB via LIANA |
| IvyGAP RNA-seq | glioblastoma.alleninstitute.org |
| CPTAC GBM proteomics | cBioPortal, study `gbm_cptac_2021` |
| TCGA PanCancer Atlas | cBioPortal |
| CGGA mRNAseq-693 and -325 | cgga.org.cn, release 2020-05-06 |
| GLASS paired RNA | Synapse `syn17038081` |
| Longitudinal case; immunotherapy cohort | GEO `GSE318238`; `GSE121810` |
| Target annotations | Open Targets Platform, DepMap |

GLASS data are governed by the GLASS data-use terms and are **not** redistributed in this
repository. Only derived summary statistics appear in the supplementary tables. Download
the raw tables from Synapse yourself before running script 20.

## Contents

```
scripts/     numbered analysis scripts, run in order
results/     processed outputs released with the paper (Supplementary Tables S1–S16)
```

Scripts are numbered by stage: 00–06 build the atlas (QC, copy number, deconvolution,
lineage composition patterns, niches); 07 scores ligand–receptor pairs; 09 fits the
niche-dependent gene models with the state and hypoxia adjustments; 10–15 run the bulk,
survival, longitudinal and immunotherapy validations; 16–21 and 29–34 build the figures.

## Running

Python 3.11+ with `scanpy`, `anndata`, `infercnvpy`, `harmonypy`, `scikit-learn`, `scipy`,
`statsmodels`, `lifelines`, `matplotlib`, `adjustText`, `pandas`, `numpy`.

```bash
pip install -r requirements.txt
python3 scripts/01_load_and_qc.py <raw_dir> <out_dir>
```

Each script takes its input and output directories as command-line arguments and prints the
numbers reported in the manuscript. Random seeds are fixed for all clustering.

Figures are restyled for publication by running a figure script through the style hook:

```bash
python3 scripts/33_pubstyle.py all scripts/16_article_figures.py
```

## Licence

Code: MIT (see `LICENSE`). Processed result tables: CC BY 4.0.

## Citation

See `CITATION.cff`. Please also cite the original data sources listed above.
