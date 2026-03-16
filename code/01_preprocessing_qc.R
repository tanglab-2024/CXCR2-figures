# ==============================================================================
# 01_preprocessing_qc.R
# Load raw 10X data, perform QC filtering, and save merged/filtered object
#
# Input:  10X Cell Ranger filtered matrices (LSEC-B1 through LSEC-B5)
# Output: LSECfiltered.RData — merged, QC-filtered Seurat object
#
# Figures: QC plots (supplementary)
# ==============================================================================

library(Seurat)
library(dplyr)
library(stringr)
library(ggplot2)

# --- Paths -------------------------------------------------------------------
data_dir <- "../data"
out_dir  <- "../results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- 1. Load 10X matrices ----------------------------------------------------
message("Loading 10X matrices...")
sample_ids <- paste0("LSEC-B", 1:5)
seurat_list <- lapply(sample_ids, function(sid) {
  mat <- Read10X(data.dir = file.path(data_dir, sid, "sample_filtered_feature_bc_matrix"))
  CreateSeuratObject(counts = mat, project = paste0(sid, "_baton"), min.features = 100)
})
names(seurat_list) <- sample_ids

# --- 2. Merge ----------------------------------------------------------------
message("Merging samples...")
LSECmerge <- merge(
  seurat_list[[1]],
  y          = seurat_list[2:5],
  add.cell.ids = paste0("B", 1:5),
  project    = "LSEC-baton_merge"
)

# --- 3. Compute QC metrics ----------------------------------------------------
LSECmerge[["percent.mt"]]       <- PercentageFeatureSet(LSECmerge, pattern = "^MT-")
LSECmerge$log10GenesPerUMI      <- log10(LSECmerge$nFeature_RNA) / log10(LSECmerge$nCount_RNA)
LSECmerge$mitoRatio             <- LSECmerge$percent.mt / 100

# Assign sample labels (monoculture vs coculture conditions)
LSECmerge$sample <- case_when(
  grepl("^B1_", colnames(LSECmerge)) ~ "Monoculture",
  grepl("^B[23]_", colnames(LSECmerge)) ~ "Coculture-1",
  grepl("^B[45]_", colnames(LSECmerge)) ~ "Coculture-2"
)

# Rename metadata columns for clarity
LSECmerge@meta.data <- LSECmerge@meta.data %>%
  dplyr::rename(nUMI = nCount_RNA, nGene = nFeature_RNA)

# --- 4. QC plots (supplementary) ---------------------------------------------
message("Generating QC plots...")
meta <- LSECmerge@meta.data

qc_plots <- list(
  # Cell counts per sample
  ggplot(meta, aes(x = sample, fill = sample)) +
    geom_bar() + theme_classic() + ggtitle("Cells per sample"),
  # UMI distribution
  ggplot(meta, aes(x = nUMI, color = sample, fill = sample)) +
    geom_density(alpha = 0.2) + scale_x_log10() + theme_classic() +
    geom_vline(xintercept = 2000, linetype = "dashed") + ylab("Cell density"),
  # Gene distribution
  ggplot(meta, aes(x = nGene, color = sample, fill = sample)) +
    geom_density(alpha = 0.2) + scale_x_log10() + theme_classic() +
    geom_vline(xintercept = 1000, linetype = "dashed"),
  # Mito ratio
  ggplot(meta, aes(x = mitoRatio, color = sample, fill = sample)) +
    geom_density(alpha = 0.2) + scale_x_log10() + theme_classic() +
    geom_vline(xintercept = 0.2, linetype = "dashed"),
  # Complexity
  ggplot(meta, aes(x = log10GenesPerUMI, color = sample, fill = sample)) +
    geom_density(alpha = 0.2) + theme_classic() +
    geom_vline(xintercept = 0.8, linetype = "dashed")
)

pdf(file.path(out_dir, "QC_plots.pdf"), width = 10, height = 7)
for (p in qc_plots) print(p)
dev.off()

# --- 5. Filter ----------------------------------------------------------------
message("Filtering low-quality cells...")
LSECfiltered <- subset(
  LSECmerge,
  subset = (nUMI >= 2000) & (nGene >= 1000) &
           (log10GenesPerUMI > 0.8) & (mitoRatio < 0.2)
)
LSECfiltered <- JoinLayers(LSECfiltered)

# Keep genes expressed in >= 10 cells
counts  <- GetAssayData(LSECfiltered, layer = "counts")
nonzero <- counts > 0
keep    <- Matrix::rowSums(nonzero) >= 10
LSECfiltered <- CreateSeuratObject(
  counts[keep, ],
  meta.data = LSECfiltered@meta.data
)

message(sprintf("After filtering: %d cells, %d genes", ncol(LSECfiltered), nrow(LSECfiltered)))

# --- 6. Save ------------------------------------------------------------------
save(LSECfiltered, file = file.path(out_dir, "LSECfiltered.RData"))
message("01_preprocessing_qc.R complete.")
