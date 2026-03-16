# ==============================================================================
# 06_twoway_volcano.R
# Two-way volcano plot comparing bulk RNA-seq DE with scRNA-seq cluster DE
#
# Input:  results/LSEC_RNAseq_result.csv (from 05)
#         results/subLSEC.rds (from 02)
# Output: Figure 1j — Two-way volcano (bulk vs scRNA-seq EC4 vs EC1)
#
# Note: Seurat clusters 0-4 correspond to EC1-EC5 in the manuscript.
# ==============================================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(scales)

out_dir <- "/results"

# --- 1. Load bulk RNA-seq DE results ------------------------------------------
message("Loading bulk DE results...")
bulk_data <- read.csv(file.path(out_dir, "LSEC_RNAseq_result.csv")) %>%
  dplyr::select(gene, log2FoldChange, padj) %>%
  dplyr::rename(bulk_log2FC = log2FoldChange, bulk_padj = padj)

# --- 2. Compute scRNA-seq DE (EC4 vs EC1 = Seurat cluster 3 vs 0) -------------
message("Computing scRNA-seq DE: EC4 vs EC1 (Seurat clusters 3 vs 0)...")
subLSEC <- readRDS(file.path(out_dir, "subLSEC.rds"))
Idents(subLSEC) <- "seurat_clusters"

# Also compute EC4 vs all others (EC1-3, EC5 = Seurat clusters 0,1,2,4)
de_genes_ec4_vs_others <- RunPresto(
  subLSEC,
  ident.1 = 3,
  ident.2 = c(0, 1, 2, 4),
  logfc.threshold = 0,
  test.use = "wilcox"
)
de_genes_ec4_vs_others$gene <- rownames(de_genes_ec4_vs_others)
write.csv(de_genes_ec4_vs_others,
          file.path(out_dir, "EC4_vs_others_DEGs.csv"), row.names = FALSE)

# EC4 vs EC1 (used for two-way volcano, Fig 1j)
de_genes_ec4_vs_ec1 <- RunPresto(
  subLSEC,
  ident.1 = 3,
  ident.2 = 0,
  logfc.threshold = 0,
  test.use = "wilcox"
)
de_genes_ec4_vs_ec1$gene <- rownames(de_genes_ec4_vs_ec1)
write.csv(de_genes_ec4_vs_ec1,
          file.path(out_dir, "EC4_vs_EC1_DEGs.csv"), row.names = FALSE)

sc_data <- de_genes_ec4_vs_ec1 %>%
  dplyr::select(gene, avg_log2FC, p_val_adj) %>%
  dplyr::rename(sc_log2FC = avg_log2FC, sc_padj = p_val_adj)

# --- 3. Merge and annotate ----------------------------------------------------
message("Merging bulk and single-cell DE results...")
merged_data <- inner_join(sc_data, bulk_data, by = "gene") %>%
  filter(!is.na(sc_padj) & !is.na(bulk_padj)) %>%
  mutate(
    # Use the less significant p-value (more conservative)
    max_pval = pmax(sc_padj, bulk_padj, na.rm = TRUE),
    max_pval = case_when(
      is.na(max_pval) ~ 1,
      max_pval == 0   ~ 1e-300,
      max_pval > 1    ~ 1,
      TRUE            ~ max_pval
    ),
    pval_exponent = pmin(-log10(max_pval), 50),

    # Concordant DE categories
    outline_category = case_when(
      sc_padj < 1e-5 & bulk_padj < 1e-5 &
        sc_log2FC > 1 & bulk_log2FC > 1   ~ "Both Up",
      sc_padj < 1e-5 & bulk_padj < 1e-5 &
        sc_log2FC < -1 & bulk_log2FC < -1 ~ "Both Down",
      TRUE ~ "Other"
    ),
    label_gene = outline_category %in% c("Both Up", "Both Down")
  )

genes_to_label <- merged_data %>% filter(label_gene)

# --- 4. Figure 1j: Two-way volcano plot --------------------------------------
message("Generating Fig 1j: Two-way volcano...")

p_twoway <- ggplot(merged_data, aes(x = bulk_log2FC, y = sc_log2FC)) +
  # Points colored by max p-value
  geom_point(aes(color = pval_exponent), alpha = 0.7, size = 1.5) +
  # Outline concordant DE genes
  geom_point(
    data = merged_data %>% filter(outline_category %in% c("Both Up", "Both Down")),
    aes(shape = outline_category),
    color = "black", size = 2, stroke = 1, fill = NA
  ) +
  # Gene labels
  geom_text_repel(
    data = genes_to_label,
    aes(label = gene),
    size = 3, max.overlaps = 20, box.padding = 0.5, point.padding = 0.2
  ) +
  # Color scale
  scale_color_gradient(
    low = "grey80", high = "red3",
    name = expression(-log[10](max~p[adj])),
    guide = guide_colorbar(title.position = "top", title.hjust = 0.5),
    limits = c(0, 5), na.value = "grey90", oob = squish
  ) +
  # Shape legend
  scale_shape_manual(
    values = c("Both Up" = 21, "Both Down" = 21),
    name = "Concordant DE",
    guide = guide_legend(override.aes = list(color = "black", fill = NA, stroke = 1))
  ) +
  # Reference lines
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", alpha = 0.7) +
  geom_hline(yintercept = c(-1, 1), linetype = "dashed", alpha = 0.7) +
  labs(
    title = "Two-way Volcano: Bulk RNA-seq vs scRNA-seq",
    subtitle = "Outlines: concordant DE (padj < 1e-5 & |LFC| > 1 in both)",
    x = "Bulk RNA-seq Log2FC (Coculture vs Monoculture)",
    y = "scRNA-seq Log2FC (EC4 vs EC1)"
  ) +
  theme_minimal() +
  theme(
    plot.title    = element_text(hjust = 0.5, size = 12),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    legend.position = "bottom",
    legend.box      = "horizontal"
  ) +
  coord_cartesian(xlim = c(-2.5, 2.5), ylim = c(-5, 5))

ggsave(file.path(out_dir, "Fig1j_twoway_volcano.pdf"), p_twoway, width = 10, height = 10)

message("06_twoway_volcano.R complete.")
