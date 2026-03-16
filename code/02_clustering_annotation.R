# ==============================================================================
# 02_clustering_annotation.R
# SCTransform normalization, PCA/UMAP, SingleR annotation, endothelial subset
#
# Input:  LSECfiltered.RData (from 01)
# Output: subLSEC.rds — endothelial cell subset, clustered and annotated
#
# Note: This script includes the intermediate step of annotating the full
#       dataset with SingleR, subsetting endothelial cells, and re-clustering.
# ==============================================================================

library(Seurat)
library(SingleR)
library(SingleCellExperiment)
library(celldex)
library(ggplot2)

out_dir <- "/results"
load(file.path(out_dir, "LSECfiltered.RData"))

# --- 1. Normalize and cluster full dataset ------------------------------------
message("Running SCTransform on full dataset...")
LSECfiltered <- SCTransform(LSECfiltered, vars.to.regress = "mitoRatio")
LSECfiltered <- RunPCA(LSECfiltered)
LSECfiltered <- FindNeighbors(LSECfiltered, dims = 1:40)
LSECfiltered <- FindClusters(LSECfiltered, algorithm = 3, resolution = 0.39)
LSECfiltered <- RunUMAP(LSECfiltered, dims = 1:40)

# --- 2. SingleR annotation (full dataset) ------------------------------------
message("Running SingleR annotation...")
ref <- celldex::HumanPrimaryCellAtlasData()
results <- SingleR(
  test   = as.SingleCellExperiment(LSECfiltered, assay = "RNA"),
  ref    = ref,
  labels = ref$label.fine
)
LSECfiltered$singler_labels <- results$labels

# --- 3. Subset endothelial cells ----------------------------------------------
message("Subsetting endothelial cells...")
endo_labels <- unique(LSECfiltered$singler_labels[
  grepl("^Endothelial", LSECfiltered$singler_labels)
])
LSECendo <- subset(LSECfiltered, subset = singler_labels %in% endo_labels)

# Re-normalize and re-cluster the endothelial subset
LSECendo <- SCTransform(LSECendo, vars.to.regress = "mitoRatio")
LSECendo <- RunPCA(LSECendo)
LSECendo <- FindNeighbors(LSECendo, dims = 1:40)
LSECendo <- FindClusters(LSECendo, algorithm = 3, resolution = 0.39)
LSECendo <- RunUMAP(LSECendo, dims = 1:40)

# --- 4. Remove small / low-quality clusters -----------------------------------
# Clusters 9-14 are small outlier populations removed for downstream analysis
message("Removing small clusters (9-14)...")
subLSEC <- subset(LSECendo, idents = c(9, 10, 11, 12, 13, 14), invert = TRUE)
subLSEC <- RunPCA(subLSEC)
subLSEC <- FindNeighbors(subLSEC, dims = 1:11)
subLSEC <- FindClusters(subLSEC, resolution = 0.17, algorithm = 3)
subLSEC <- RunUMAP(subLSEC, dims = 1:11)

# --- 5. Re-annotate subset with SingleR ---------------------------------------
message("Annotating endothelial subset with SingleR...")
results_sub <- SingleR(
  test   = as.SingleCellExperiment(subLSEC, assay = "RNA"),
  ref    = ref,
  labels = ref$label.fine
)
subLSEC$singler_labels     <- results_sub$labels
subLSEC$singler_labels_abb <- sub("^Endothelial_cells:", "", subLSEC$singler_labels)

# Set sample factor levels for consistent plotting
subLSEC$sample <- factor(
  subLSEC$sample,
  levels = c("Monoculture", "Coculture-1", "Coculture-2")
)

message(sprintf("Final endothelial subset: %d cells, %d genes",
                ncol(subLSEC), nrow(subLSEC)))

# --- 6. Save ------------------------------------------------------------------
saveRDS(subLSEC, file.path(out_dir, "subLSEC.rds"))
message("02_clustering_annotation.R complete.")
