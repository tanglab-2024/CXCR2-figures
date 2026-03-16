#!/usr/bin/env bash
# ==============================================================================
# Master run script for Code Ocean capsule
# Manuscript: CXCL-CXCR2 signaling drives cancer-endothelium interactions
#             in SCLC metastatic seeding
# Manuscript ID: CD-26-0546
# ==============================================================================

set -euo pipefail

echo "============================================="
echo "Starting analysis pipeline"
echo "Date: $(date)"
echo "============================================="

# Step 1: QC, filtering, and preprocessing
echo "[1/4] Preprocessing and QC..."
Rscript code/01_preprocessing_qc.R

# Step 2: Clustering and cell type annotation
echo "[2/4] Clustering and annotation..."
Rscript code/02_clustering_annotation.R

# Step 3: Differential expression and visualization (Figures 1b-g)
echo "[3/4] DE analysis and figure generation..."
Rscript code/03_de_visualization.R

# Step 4: GSEA / pathway enrichment (Figure 1h)
echo "[4/4] Gene set enrichment analysis..."
Rscript code/04_gsea_enrichment.R

echo "============================================="
echo "Pipeline complete. Results in ../results/"
echo "============================================="
