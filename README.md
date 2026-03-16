# Code Ocean Capsule: CXCL-CXCR2 Signaling in SCLC Metastatic Seeding

**Manuscript:** CXCL-CXCR2 signaling drives cancer-endothelium interactions in
SCLC metastatic seeding  
**Manuscript ID:** CD-26-0546  
**Journal:** Cancer Discovery

## Overview

Code to reproduce the figures from scRNA-seq and bulk RNA-seq of human liver sinusoidal endothelial cells (LSECs) under monoculture and coculture conditions with SCLC cells (Figure 1), MOBA-seq barcode tracing of CXCR2 knockout (Figure 3), Rac1 knockout metastatic analysis (Figure 5), and CXCR2 pharmacological inhibition (Figure 6).

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
| `11a_export_for_scvelo.R` | Export Seurat metadata for scVelo | — |
| `11b_scvelo.py` | RNA velocity stream plots (scVelo) | 1d |

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

**RNA velocity (velocyto loom files):**
- `Loomfiles/sample_alignments_B1.loom` — Monoculture
- `Loomfiles/sample_alignments_B2.loom` — Coculture-1 replicate A
- `Loomfiles/sample_alignments_B3.loom` — Coculture-1 replicate B
- `Loomfiles/sample_alignments_B4.loom` — Coculture-2 replicate A
- `Loomfiles/sample_alignments_B5.loom` — Coculture-2 replicate B

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

All outputs are written to `/results/`.

## Environment

- R 4.3.2 with Seurat v5, SingleR, escape, ComplexHeatmap, presto
- See `environment/Dockerfile` for the full dependency list

## Contact

[Your name and email]
