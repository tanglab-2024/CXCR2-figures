# ==============================================================================
# 03_de_visualization.R
# Differential expression, UMAP, composition plots, heatmaps, violin/feature plots
#
# Input:  subLSEC.rds (from 02)
# Output: Figures 1b-g
#   - Fig 1b: UMAP and split UMAP
#   - Fig 1c: Sample composition bar plot
#   - Fig 1e: ComplexHeatmap of top cluster markers
#   - Fig 1f: Violin plots of CXCL genes
#   - Fig 1g: Split feature plots of CXCL genes
# ==============================================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(ggalluvial)
library(presto)
library(ComplexHeatmap)
library(circlize)
library(patchwork)

out_dir <- "../results"
subLSEC <- readRDS(file.path(out_dir, "subLSEC.rds"))

# ==============================================================================
# Figure 1b: UMAP plots
# ==============================================================================
message("Generating Fig 1b: UMAP plots...")

total_cells <- ncol(subLSEC)
avg_umi     <- round(mean(subLSEC$nCount_RNA), 2)

# Full UMAP
p_umap <- DimPlot(subLSEC, reduction = "umap", label = TRUE) +
  annotate("text", x = Inf, y = -Inf, hjust = 1.2, vjust = -3, size = 4,
           label = paste("Total Cells:", total_cells)) +
  annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -1.5, size = 4,
           label = paste("Avg UMI/Cell:", avg_umi))
ggsave(file.path(out_dir, "Fig1b_UMAP.pdf"), p_umap, width = 10, height = 7)

# Split UMAP by sample
p_umap_split <- DimPlot(subLSEC, reduction = "umap", label = TRUE, split.by = "sample")
ggsave(file.path(out_dir, "Fig1b_UMAP_split.pdf"), p_umap_split, width = 14, height = 7)

# SingleR annotation UMAP
p_singler <- DimPlot(subLSEC, reduction = "umap",
                     group.by = "singler_labels_abb", label = TRUE)
ggsave(file.path(out_dir, "Fig1b_SingleR_UMAP.pdf"), p_singler, width = 14, height = 7)

# ==============================================================================
# Figure 1c: Sample composition bar plot (alluvial)
# ==============================================================================
message("Generating Fig 1c: Composition plot...")

cluster_counts <- as.data.frame(table(subLSEC$seurat_clusters, subLSEC$sample))
colnames(cluster_counts) <- c("Cluster", "Sample", "Count")

cluster_pct <- cluster_counts %>%
  group_by(Sample) %>%
  mutate(Total = sum(Count), Percent = (Count / Total) * 100) %>%
  ungroup()

p_comp <- ggplot(cluster_pct, aes(x = Sample, y = Percent, fill = Cluster)) +
  geom_flow(aes(alluvium = Cluster), alpha = 0.5, color = "black",
            curve_type = "linear", width = 0.5) +
  geom_col(width = 0.5, color = "black") +
  geom_text(aes(label = sprintf("%.1f%%", Percent)),
            position = position_stack(vjust = 0.5), size = 3) +
  cowplot::theme_minimal_hgrid() +
  labs(title = "Sample Composition by Cluster",
       x = "Sample", y = "Percent Composition", fill = "Cluster") +
  theme(panel.grid.major = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

ggsave(file.path(out_dir, "Fig1c_composition.pdf"), p_comp, width = 14, height = 7)

# ==============================================================================
# Figure 1e: Heatmap of top cluster markers (ComplexHeatmap)
# ==============================================================================
message("Generating Fig 1e: Marker heatmap...")

# Find markers using presto (fast Wilcoxon)
markers <- RunPrestoAll(subLSEC)

top5 <- markers %>%
  group_by(cluster) %>%
  filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup()

# Build expression matrix (downsampled for visualization)
sub_down <- subset(subLSEC, downsample = 2000)
mat <- GetAssayData(sub_down, assay = "RNA", layer = "data")[top5$gene, ] %>% as.matrix()
mat <- t(scale(t(mat)))

cluster_anno <- sub_down$seurat_clusters
n_clusters   <- length(unique(cluster_anno))

col_fun <- colorRamp2(c(-1, 0, 3), c("#FF00FF", "black", "#FFFF00"))

pdf(file.path(out_dir, "Fig1e_heatmap.pdf"), width = 10, height = 7)
Heatmap(
  mat,
  name              = "Expression",
  col               = col_fun,
  column_split      = factor(cluster_anno[match(colnames(mat), names(cluster_anno))]),
  cluster_columns   = FALSE,
  show_column_dend  = FALSE,
  cluster_column_slices = TRUE,
  column_title_gp   = gpar(fontsize = 8),
  column_gap        = unit(0.5, "mm"),
  cluster_rows      = TRUE,
  show_row_dend     = FALSE,
  row_names_gp      = gpar(fontsize = 4),
  column_title_rot  = 90,
  top_annotation    = HeatmapAnnotation(
    foo = anno_block(gp = gpar(fill = scales::hue_pal()(n_clusters)))
  ),
  show_column_names = FALSE,
  use_raster        = FALSE
)
dev.off()

# ==============================================================================
# Figure 1f: Violin plots of CXCL ligands
# ==============================================================================
message("Generating Fig 1f: Violin plots...")

cxcl_genes <- c("CXCL1", "CXCL2", "CXCL3", "CXCL8")

p_vln <- VlnPlot(subLSEC, features = cxcl_genes, pt.size = 0,
                  ncol = 2, layer = "counts", log = TRUE)
ggsave(file.path(out_dir, "Fig1f_violin.pdf"), p_vln, width = 14, height = 7)

# ==============================================================================
# Figure 1g: Split feature plots of CXCL ligands
# ==============================================================================
message("Generating Fig 1g: Feature plots...")

sample_subsets <- list(
  Monoculture  = subset(subLSEC, subset = sample == "Monoculture"),
  `Coculture-1` = subset(subLSEC, subset = sample == "Coculture-1"),
  `Coculture-2` = subset(subLSEC, subset = sample == "Coculture-2")
)

for (gene in cxcl_genes) {
  plots <- lapply(sample_subsets, function(ss) {
    FeaturePlot(ss, features = gene, min.cutoff = "q10")
  })
  combined <- plot_grid(plotlist = plots, ncol = 3,
                        labels = names(sample_subsets))
  ggsave(file.path(out_dir, paste0("Fig1g_", gene, "_feature.pdf")),
         combined, width = 16, height = 6)
}

message("03_de_visualization.R complete.")
