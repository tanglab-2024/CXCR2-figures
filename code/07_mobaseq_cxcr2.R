# ==============================================================================
# 07_mobaseq_cxcr2.R
# MOBA-seq analysis of CXCR2 knockout metastatic dynamics
#
# Input:  data/Moba1st_FilteredSampleInfo.csv — MOBA-seq barcode count data
# Output: Figure 3b — Seeding comparison (Ctrl vs CXCR2, jitter by tissue/time)
#         Figure 3c — Bootstrap metastatic burden (gene-level)
#         Figure 3d — Bootstrap metastatic burden (sgid-level)
#         Figure 3e — Dormancy cutoff tissue plot
#         Figure 3f — Dormancy analysis
#         Figure 3g — Supermet jitter plot
#
# Requires: mobaseq package (bundled in environment/mobaseq/)
#
# Note: Seurat clusters 0-4 correspond to EC1-EC5 in the manuscript.
# ==============================================================================

library(mobaseq)
library(dplyr)
library(ggplot2)

out_dir  <- "../results"
data_dir <- "../data"

# ==============================================================================
# 1. Load and filter MOBA-seq data
# ==============================================================================
message("Loading MOBA-seq count data...")
counts_df <- load_counts(
  file.path(data_dir, "Moba1st_FilteredSampleInfo.csv"),
  remove.sgids = "none",
  keep.distance = 0
)

# Pre-transplant histogram and filtering
result <- plot_pretran_histogram(
  counts_df,
  percentile.threshold = 0.99,
  filter.data = TRUE,
  return.list = TRUE
)
ggsave(file.path(out_dir, "pretran_histogram.pdf"),
       plot = result$plot, width = 8, height = 8)
counts_filtered <- result$filtered_df

# Filter to CXCR2 sgRNAs and controls, assign gene labels
counts_filtered_CXCR <- counts_filtered %>%
  filter(sgid %in% c("sgDummy", "sgCXCR2-1", "sgCXCR2-2", "sgCXCR2-5")) %>%
  mutate(gene = case_when(
    grepl("Dummy", sgid, ignore.case = TRUE) ~ "Ctrl",
    grepl("CXCR2", sgid, ignore.case = TRUE) ~ "CXCR2",
    TRUE ~ NA_character_
  ))

# Apply tissue-specific cell number thresholds
counts_filtered_CXCR <- counts_filtered_CXCR %>%
  filter(
    (tissue == "preTran" & cell_num >= 2) |
      (tissue == "Blood" & cell_num >= 1) |
      (tissue == "Brain" & cell_num >= 10) |
      (!tissue %in% c("preTran", "Blood", "Brain") & cell_num >= 10)
  )

# ==============================================================================
# 2. Figure 3b: Seeding comparison (jitter by tissue and time point)
# ==============================================================================
message("Generating Fig 3b: Seeding comparison...")
plot_cell_num_jitter_faceted(
  counts_filtered_CXCR,
  facet.by = "tissue",
  color.by = "gene",
  time.points = c("2days", "1week", "2weeks", "3weeks"),
  filter.facets = c("Liver", "Lung", "Brain"),
  filter.colors = c("Ctrl", "CXCR2"),
  per.mouse.annotation = TRUE,
  title = "Cell Number Seeding Comparison: Ctrl vs CXCR2",
  facet.nrow = 1
)
ggsave(file.path(out_dir, "Fig3b_seeding_comp.pdf"),
       width = 14, height = 8, dpi = 300)

# ==============================================================================
# 3. Compute metastatic statistics (bootstrap)
# ==============================================================================
message("Computing pre-transplant data...")
pretran_data <- generate_pretran_data(counts_filtered_CXCR)

message("Computing dormancy cutoffs...")
dormancy_results <- dormancy_cutoff_by_genotype_tissue_time_v4(
  counts_filtered_CXCR,
  tissues = c("Liver", "Lung", "Brain"),
  time.points = c("3weeks"),
  min.colonies = 1,
  bw.adjust = 3,
  split.by.tissue = FALSE
)

geno_cutoffs <- dormancy_results$cutoffs_df %>%
  dplyr::select(mouse_genotype, dormant_peak, dormancy_cutoff)

message("Computing bootstrap metastatic statistics (gene-level)...")
results_gene <- calculate_metastatic_statistics(
  counts_filtered_CXCR,
  tissues = c("Liver", "Lung", "Brain"),
  statistics = c("seeding", "burden", "size_percentile",
                 "supermet", "dormancy_cutoff", "peak_mode"),
  pretran.cell.num = pretran_data$pretran_cell_num,
  pretran.cell.num.gene = pretran_data$pretran_cell_num_gene,
  min.colonies.seeding = 1,
  min.colonies.burden = 1,
  R = 1000,
  dormancy.cutoff = geno_cutoffs
)
write.csv(results_gene,
          file.path(out_dir, "moba_cxcr2_gene_statistics.csv"),
          row.names = FALSE)

message("Computing bootstrap metastatic statistics (sgid-level)...")
results_sgid <- calculate_metastatic_statistics(
  counts_filtered_CXCR,
  tag = "sgid",
  tissues = c("Liver", "Lung", "Brain"),
  statistics = c("seeding", "burden", "size_percentile",
                 "supermet", "dormancy_cutoff", "peak_mode"),
  pretran.cell.num = pretran_data$pretran_cell_num,
  pretran.cell.num.gene = pretran_data$pretran_cell_num_gene,
  R = 1000,
  dormancy.cutoff = geno_cutoffs
)
write.csv(results_sgid,
          file.path(out_dir, "moba_cxcr2_sgid_statistics.csv"),
          row.names = FALSE)

# ==============================================================================
# 4. Figures 3c-d: Bootstrap metastatic burden plots
# ==============================================================================
message("Generating Figs 3c-d: Bootstrap burden plots...")

plots <- plot_all_bootstrap_combinations(
  results.df.sgid = results_sgid,
  results.df.gene = results_gene,
  plot.type = "point"
)

# Fig 3c: Gene-level burden
ggsave(file.path(out_dir, "Fig3c_burden_gene.pdf"),
       plot = plots$met_burden_gene, width = 10, height = 8, dpi = 300)

# Fig 3d: sgid-level burden
ggsave(file.path(out_dir, "Fig3d_burden_sgid.pdf"),
       plot = plots$met_burden_sgid, width = 10, height = 8, dpi = 300)

# ==============================================================================
# 5. Figure 3e: Dormancy cutoff tissue plot
# ==============================================================================
message("Generating Fig 3e: Dormancy cutoff plot...")
plot_dormancy_cutoffs_v2(
  counts.df = counts_filtered_CXCR,
  cutoffs_df = dormancy_results$cutoffs_df,
  value.col = "log2_cell_num",
  time.points = c("3weeks"),
  tissues = c("Liver", "Lung", "Brain"),
  genotypes = unique(counts_filtered_CXCR$mouse_genotype),
  split.by.tissue = FALSE,
  min.colonies = 1
)
ggsave(file.path(out_dir, "Fig3e_dormant_tissues.pdf"),
       width = 8, height = 8, dpi = 300)

# ==============================================================================
# 6. Figure 3f: Dormancy analysis
# ==============================================================================
# TODO: Confirm time points and parameters for manuscript figure
message("Generating Fig 3f: Dormancy analysis...")
results_dormancy <- plot_dormancy_analysis_v2(
  data = counts_filtered_CXCR,
  tissues = c("Liver", "Lung", "Brain"),
  genotypes = c("NSG"),
  time.points = c("3weeks"),
  include.ctrl = TRUE,
  use.constraints = TRUE,
  tag.type = "gene",
  constraints = geno_cutoffs,
  facet.by = c("tissue", "gene"),
  show.mode.lines = TRUE,
  min.colonies = 100
)
ggsave(file.path(out_dir, "Fig3f_dormancy.pdf"),
       width = 10, height = 8, dpi = 300)

# ==============================================================================
# 7. Figure 3g: Supermet jitter plot
# ==============================================================================
# TODO: Confirm target tissue for manuscript figure
message("Generating Fig 3g: Supermet jitter...")
plot_supermet_jitter(
  counts_filtered_CXCR,
  genes = c("Ctrl", "CXCR2"),
  time.points = c("2days", "1week", "2weeks", "3weeks"),
  target.tissue = "Liver"
)
ggsave(file.path(out_dir, "Fig3g_supermet.pdf"),
       width = 10, height = 8, dpi = 300)

message("07_mobaseq_cxcr2.R complete.")
