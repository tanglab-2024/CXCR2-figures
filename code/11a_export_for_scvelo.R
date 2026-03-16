# ==============================================================================
# 11a_export_for_scvelo.R
# Export Seurat metadata for scVelo RNA velocity analysis
#
# Input:  /results/subLSEC.rds (from 02)
# Output: /results/scvelo/cellID_obs.csv
#         /results/scvelo/cell_embeddings.csv
#         /results/scvelo/clusters.csv
#
# Run this BEFORE 11b_scvelo.py
# ==============================================================================

library(Seurat)

out_dir <- "/results/scvelo"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

subLSEC <- readRDS("/results/subLSEC.rds")

# Export cell IDs, UMAP embeddings, and cluster assignments
write.csv(Cells(subLSEC), file = file.path(out_dir, "cellID_obs.csv"), row.names = FALSE)
write.csv(Embeddings(subLSEC, reduction = "umap"), file = file.path(out_dir, "cell_embeddings.csv"))
write.csv(subLSEC@meta.data$seurat_clusters, file = file.path(out_dir, "clusters.csv"))

message("11a_export_for_scvelo.R complete.")
