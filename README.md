# Code Ocean Capsule: CXCL-CXCR2 Signaling in SCLC Metastatic Seeding

**Manuscript:** CXCL-CXCR2 signaling drives cancer-endothelium interactions in
SCLC metastatic seeding  
**Manuscript ID:** CD-26-0546  
**Journal:** Cancer Discovery

## Overview

This capsule contains the code to reproduce the scRNA-seq and bulk RNA-seq
analyses and figures presented in Figure 1 of the manuscript. Human liver
sinusoidal endothelial cells (LSECs) were profiled by 10X Chromium scRNA-seq
and bulk RNA-seq under monoculture and coculture conditions with SCLC cells.

## Pipeline Structure

| Script | Description | Figures |
|--------|-------------|---------|
| `01_preprocessing_qc.R` | Load 10X data, compute QC metrics, filter | QC (supplementary) |
| `02_clustering_annotation.R` | SCTransform, PCA/UMAP, SingleR, endothelial subset | — |
| `03_de_visualization.R` | DE markers, UMAP, composition, heatmap, violin, feature plots | 1b, 1c, 1e, 1f, 1g |
| `04_gsea_enrichment.R` | Hallmark pathway enrichment with escape/UCell | 1h |
| `05_bulk_rnaseq_volcano.R` | Bulk RNA-seq DESeq2 DE and volcano plot | 1i |
| `06_twoway_volcano.R` | Two-way volcano: bulk vs scRNA-seq cluster DE | 1j |
| `07_mobaseq_cxcr2.R` | MOBA-seq CXCR2 KO metastatic dynamics (RP48) | 3b, 3c, 3d, 3e, 3f, 3g |
| `08_mobaseq_h82.R` | MOBA-seq H82 cell line analysis + cross-dataset radar | 3i, 3j, 3k, 3l |
| `09_mobaseq_rac1.R` | MOBA-seq Rac1 KO metastatic dynamics | 5n, 5o, 5p, 5q |
| `10_mobaseq_cxcr_inhibitor.R` | MOBA-seq CXCR2 inhibitor (AZD) treatment | 6l, 6m, 6n |

## Data

Place the following in the `data/` directory:

**scRNA-seq (10X Chromium):**
- `LSEC-B1/sample_filtered_feature_bc_matrix/` — Monoculture
- `LSEC-B2/sample_filtered_feature_bc_matrix/` — Coculture-1 replicate A
- `LSEC-B3/sample_filtered_feature_bc_matrix/` — Coculture-1 replicate B
- `LSEC-B4/sample_filtered_feature_bc_matrix/` — Coculture-2 replicate A
- `LSEC-B5/sample_filtered_feature_bc_matrix/` — Coculture-2 replicate B

Each directory should contain `barcodes.tsv.gz`, `features.tsv.gz`, and
`matrix.mtx.gz` from Cell Ranger.

**Bulk RNA-seq:**
- `20200123_H82-LSEC_Normalized_Counts.txt` — normalized count matrix
- `LSEC_meta.txt` — sample metadata (tab-separated, with `sample` column)

**MOBA-seq:**
- `Moba1st_FilteredSampleInfo.csv` — RP48 barcode count data with sample metadata
- `mobaH82_FilteredSampleInfo.csv` — H82 barcode count data with sample metadata
- `MobaV_FilteredSampleInfo.csv` — Rac1 (MobaV) barcode count data with sample metadata
- `CXCR_inhibitor_FilteredSampleInfo.txt` — CXCR2 inhibitor (AZD) treatment barcode data

## Environment Setup

The `mobaseq` R package is bundled locally in `environment/mobaseq/`.
To prepare this directory, download the package source from the lab's
GitHub repository and place the unzipped package folder there.

## Running

```bash
bash code/run.sh
```

All outputs are written to `results/`.

## Environment

- R 4.3.2 with Seurat v5, SingleR, escape, ComplexHeatmap, presto
- See `environment/Dockerfile` for the full dependency list

## Contact

Andy Xu, xu.a@wustl.edu
