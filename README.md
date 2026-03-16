# Code: CXCL-CXCR2 Signaling in SCLC Metastatic Seeding

**Manuscript:** CXCL-CXCR2 signaling drives cancer-endothelium interactions in
SCLC metastatic seeding  
**Manuscript ID:** CD-26-0546  
**Journal:** Cancer Discovery

## Overview

This capsule contains the code to reproduce the scRNA-seq analyses and
figures presented in Figure 1 of the manuscript. Human liver sinusoidal
endothelial cells (LSECs) were profiled by 10X Chromium scRNA-seq under
monoculture and coculture conditions with SCLC cells.

## Pipeline Structure

| Script | Description | Figures |
|--------|-------------|---------|
| `01_preprocessing_qc.R` | Load 10X data, compute QC metrics, filter | QC (supplementary) |
| `02_clustering_annotation.R` | SCTransform, PCA/UMAP, SingleR, endothelial subset | — |
| `03_de_visualization.R` | DE markers, UMAP, composition, heatmap, violin, feature plots | 1b, 1c, 1e, 1f, 1g |
| `04_gsea_enrichment.R` | Hallmark pathway enrichment with escape/UCell | 1h |

## Data

Place the following in the `data/` directory:

- `LSEC-B1/sample_filtered_feature_bc_matrix/` — Monoculture
- `LSEC-B2/sample_filtered_feature_bc_matrix/` — Coculture-1 replicate A
- `LSEC-B3/sample_filtered_feature_bc_matrix/` — Coculture-1 replicate B
- `LSEC-B4/sample_filtered_feature_bc_matrix/` — Coculture-2 replicate A
- `LSEC-B5/sample_filtered_feature_bc_matrix/` — Coculture-2 replicate B

Each directory should contain `barcodes.tsv.gz`, `features.tsv.gz`, and
`matrix.mtx.gz` from Cell Ranger.

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
