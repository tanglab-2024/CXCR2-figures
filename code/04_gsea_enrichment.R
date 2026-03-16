# ==============================================================================
# 04_gsea_enrichment.R
# Gene set enrichment analysis using escape (UCell) on Hallmark gene sets
#
# Input:  subLSEC.rds (from 02)
# Output: Figure 1h — Hallmark pathway enrichment heatmap across clusters
# ==============================================================================

library(Seurat)
library(SingleCellExperiment)
library(escape)
library(msigdbr)
library(ggplot2)

out_dir <- "../results"
subLSEC <- readRDS(file.path(out_dir, "subLSEC.rds"))

# --- 1. Convert to SingleCellExperiment ---------------------------------------
message("Converting to SCE and fetching Hallmark gene sets...")
sce <- as.SingleCellExperiment(subLSEC, assay = "RNA")
GS.hallmark <- getGeneSets(library = "H")

# --- 2. Run escape enrichment scoring ----------------------------------------
message("Running UCell enrichment scoring (this may take several minutes)...")
sce <- runEscape(
  sce,
  method         = "UCell",
  gene.sets      = GS.hallmark,
  groups         = 5000,
  min.size       = 5,
  new.assay.name = "escape.UCell"
)

# --- 3. Generate enrichment heatmap (Figure 1h) ------------------------------
message("Generating Fig 1h: Enrichment heatmap...")

pdf(file.path(out_dir, "Fig1h_escape_heatmap.pdf"), width = 12, height = 8)
heatmapEnrichment(
  sce,
  group.by     = "ident",
  assay        = "escape.UCell",
  scale        = TRUE,
  cluster.rows = TRUE,
  palette      = "Blue-Red 2"
)
dev.off()

# --- 4. Optional: rank and plot top 10 pathways ------------------------------
message("Ranking top enriched pathways...")
enrich_mat <- escape.matrix(
  sce,
  method    = "UCell",
  gene.sets = GS.hallmark,
  groups    = 5000,
  min.size  = 5
)

mean_scores  <- apply(enrich_mat, 1, mean)
top10_sets   <- names(sort(mean_scores, decreasing = TRUE))[1:10]

# Save ranked pathway scores
write.csv(
  data.frame(pathway = names(sort(mean_scores, decreasing = TRUE)),
             mean_score = sort(mean_scores, decreasing = TRUE)),
  file.path(out_dir, "hallmark_pathway_rankings.csv"),
  row.names = FALSE
)

message("04_gsea_enrichment.R complete.")
