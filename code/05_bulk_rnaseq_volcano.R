# ==============================================================================
# 05_bulk_rnaseq_volcano.R
# Bulk RNA-seq: DESeq2 differential expression and volcano plot
#
# Input:  data/20200123_H82-LSEC_Normalized_Counts.txt — normalized count matrix
#         data/LSEC_meta.txt — sample metadata
# Output: Figure 1i — Volcano plot (coculture vs monoculture, LFC > 1 cutoff)
#         LSEC_RNAseq_result.csv — full DESeq2 results
# ==============================================================================

library(DESeq2)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(grid)
library(gtable)

out_dir <- "../results"
data_dir <- "../data"

# --- Helper: add flagged labels to pheatmap -----------------------------------
# Displays only selected row labels with connector lines (for dense heatmaps)
add.flag <- function(pheatmap, kept.labels, repel.degree) {
  heatmap <- pheatmap$gtable
  new.label <- heatmap$grobs[[which(heatmap$layout$name == "row_names")]]
  new.label$label <- ifelse(new.label$label %in% kept.labels,
                            new.label$label, "")

  repelled.y <- function(d, d.select, k = repel.degree) {
    strip.npc <- function(dd) {
      if (!"unit.arithmetic" %in% class(dd)) return(as.numeric(dd))
      d1 <- strip.npc(dd$arg1)
      d2 <- strip.npc(dd$arg2)
      fn <- dd$fname
      return(lazyeval::lazy_eval(paste(d1, fn, d2)))
    }
    full.range     <- sapply(seq_along(d), function(i) strip.npc(d[i]))
    selected.range <- sapply(seq_along(d[d.select]), function(i) strip.npc(d[d.select][i]))
    return(unit(seq(
      from = max(selected.range) + k * (max(full.range) - max(selected.range)),
      to   = min(selected.range) - k * (min(selected.range) - min(full.range)),
      length.out = sum(d.select)
    ), "npc"))
  }

  new.y.positions <- repelled.y(new.label$y, d.select = new.label$label != "")
  new.flag <- segmentsGrob(
    x0 = new.label$x,
    x1 = new.label$x + unit(0.15, "npc"),
    y0 = new.label$y[new.label$label != ""],
    y1 = new.y.positions
  )
  new.label$x <- new.label$x + unit(0.15, "npc")
  new.label$y[new.label$label != ""] <- new.y.positions
  heatmap <- gtable::gtable_add_grob(x = heatmap, grobs = new.flag, t = 4, l = 4)
  heatmap$grobs[[which(heatmap$layout$name == "row_names")]] <- new.label
  grid.newpage()
  grid.draw(heatmap)
  invisible(heatmap)
}

# --- 1. Load and clean count data ---------------------------------------------
message("Loading count data...")
data <- read.table(
  file.path(data_dir, "20200123_H82-LSEC_Normalized_Counts.txt"),
  sep = "\t", header = TRUE
)

# Keep only genes with unique identifiers
unique_values <- unique(data[, 1])
filtered_data <- data.frame()
for (value in unique_values) {
  subset_data <- data[data[, 1] == value, ]
  if (nrow(subset_data) == 1) {
    filtered_data <- rbind(filtered_data, subset_data)
  }
}
rownames(filtered_data) <- filtered_data[, 1]
filtered_data <- filtered_data[, -1]

# Select LSEC columns (columns 10 onward)
countdata <- filtered_data[, 10:ncol(filtered_data)]

# --- 2. Heatmap with CXCL labels (supplementary) -----------------------------
message("Generating bulk RNA-seq heatmap...")
df <- countdata[which(rowSums(countdata) > 10), ]
df <- na.omit(df)

pdf(file.path(out_dir, "bulk_rnaseq_heatmap.pdf"), width = 10, height = 12)
heatmap_obj <- pheatmap(
  df,
  color = colorRampPalette(c("magenta", "black", "yellow"))(50),
  scale = "row", show_rownames = FALSE, clustering_method = "ward.D"
)

# Extract top up/down regulated genes by clustering order
order <- heatmap_obj$tree_row$order
reordered_df <- df[order, ]
rows_to_keep <- 1227
rows_total   <- nrow(reordered_df)
mat_heatmap  <- as.matrix(reordered_df)
mat_heatmap  <- mat_heatmap[c(1:rows_to_keep, (rows_total - 1379 + 1):rows_total), ]

z <- t(scale(t(mat_heatmap)))
heat_lab <- pheatmap(
  z,
  color = colorRampPalette(c("magenta", "black", "yellow"))(50),
  legend_breaks = c(-2, -1, 0, 1, 2, max(z)),
  legend_labels = c("-2", "-1", "0", "1", "2", "z-score\n"),
  show_rownames = TRUE, clustering_method = "ward.D"
)
add.flag(heat_lab, kept.labels = paste0("CXCL", 1:17), repel.degree = 0.2)
dev.off()

# --- 3. DESeq2 differential expression ---------------------------------------
message("Running DESeq2...")
data_clean <- na.omit(countdata)
Metadata   <- read.csv(file.path(data_dir, "LSEC_meta.txt"), sep = "\t", header = TRUE)

dds <- DESeqDataSetFromMatrix(
  countData = data_clean,
  colData   = Metadata,
  design    = ~sample
)
dds <- dds[rowSums(counts(dds)) > 1]
ds2 <- DESeq(dds, betaPrior = FALSE)

# Export results
ds2_results <- results(ds2, contrast = c("sample", "Coculture", "Monoculture"))
voldata <- as.data.frame(ds2_results)
voldata$gene <- rownames(voldata)

write.csv(voldata, file.path(out_dir, "LSEC_RNAseq_result.csv"), row.names = FALSE)

# --- 4. Figure 1i: Volcano plot (LFC > 1 cutoff) -----------------------------
message("Generating Fig 1i: Volcano plot...")

voldata$label <- ifelse(
  (voldata$padj < 0.05) & (abs(voldata$log2FoldChange) > 1), "Yes", "No"
)

n_up   <- voldata %>% filter(padj < 0.05, log2FoldChange > 0) %>% nrow()
n_down <- voldata %>% filter(padj < 0.05, log2FoldChange < 0) %>% nrow()

cxcl_genes <- voldata %>%
  filter(gene %in% c("CXCL1", "CXCL2", "CXCL3")) %>%
  mutate(neg_log_padj = -log10(padj))

p_volcano <- ggplot(voldata, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = label)) +
  scale_color_manual(values = c("#d2dae2", "#b50719")) +
  geom_text_repel(
    data = cxcl_genes,
    aes(x = log2FoldChange, y = neg_log_padj, label = gene),
    size = 3, fontface = "bold", nudge_y = -1.5
  ) +
  labs(
    title = "Volcano Plot: Coculture vs Monoculture",
    x = expression(log[2](FC)),
    y = expression(-log[10](padj))
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = 4, color = "gray50") +
  geom_vline(xintercept = c(-1, 1), linetype = 4, color = "gray50") +
  annotate("text", x = -2.65, y = -log10(0.05) + 2,
           label = "padj = 0.05", size = 3, hjust = 0, color = "gray50") +
  annotate("text", x = -1, y = 145, label = "LFC = -1",
           size = 3, hjust = 1.1, color = "gray50") +
  annotate("text", x = 1, y = 145, label = "LFC = 1",
           size = 3, hjust = -0.1, color = "gray50") +
  annotate("text", x = 2.4, y = 140,
           label = paste0("Up: ", n_up, "\n(p < 0.05, LFC > 0)"),
           size = 3, hjust = 1, color = "black") +
  annotate("text", x = -2.4, y = 140,
           label = paste0("Down: ", n_down, "\n(p < 0.05, LFC < 0)"),
           size = 3, hjust = 0, color = "black") +
  coord_cartesian(ylim = c(0, 150), xlim = c(-2.5, 2.5)) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "none")

ggsave(file.path(out_dir, "Fig1i_volcano.pdf"), p_volcano, width = 8, height = 8)

message("05_bulk_rnaseq_volcano.R complete.")
