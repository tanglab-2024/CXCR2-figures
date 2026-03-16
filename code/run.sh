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
echo "[4/6] Gene set enrichment analysis..."
Rscript code/04_gsea_enrichment.R

# Step 5: Bulk RNA-seq DE and volcano plot (Figure 1i)
echo "[5/6] Bulk RNA-seq analysis..."
Rscript code/05_bulk_rnaseq_volcano.R

# Step 6: Two-way volcano — bulk vs scRNA-seq (Figure 1j)
echo "[6/7] Two-way volcano plot..."
Rscript code/06_twoway_volcano.R

# Step 7: MOBA-seq CXCR2 metastatic analysis — RP48 (Figure 3b-g)
echo "[7/8] MOBA-seq RP48 analysis..."
Rscript code/07_mobaseq_cxcr2.R

# Step 8: MOBA-seq H82 analysis (Figure 3i-l)
echo "[8/9] MOBA-seq H82 analysis..."
Rscript code/08_mobaseq_h82.R

# Step 9: MOBA-seq Rac1 analysis (Figure 5n-q)
echo "[9/10] MOBA-seq Rac1 analysis..."
Rscript code/09_mobaseq_rac1.R

# Step 10: MOBA-seq CXCR2 inhibitor analysis (Figure 6l-n)
echo "[10/10] MOBA-seq CXCR2 inhibitor analysis..."
Rscript code/10_mobaseq_cxcr_inhibitor.R

echo "============================================="
echo "Pipeline complete. Results in /results/"
echo "============================================="
