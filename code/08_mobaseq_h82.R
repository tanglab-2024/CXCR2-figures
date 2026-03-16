# ==============================================================================
# 08_mobaseq_h82.R
# MOBA-seq analysis of H82 cell line with CXCR2/CXCR4/NF1B knockouts
#
# Input:  data/mobaH82_FilteredSampleInfo.csv — H82 MOBA-seq barcode count data
# Output: Figure 3i  — Colony jitter plot (Liver, 3 weeks)
#         Figure 3j  — Burden dot plot (Liver, 3 weeks)
#         Figure 3k  — Seeding dot plot (Liver, 3 weeks)
#         Figure 3l  — Cross-dataset radar comparison (H82 vs RP48, Liver)
#
# Requires: mobaseq package (bundled in environment/mobaseq/)
#
# Note: This script also requires radar_ready_dataC from 07_mobaseq_cxcr2.R
#       for the cross-dataset radar plot (Fig 3l). Run 07 before this script.
# ==============================================================================

library(mobaseq)
library(dplyr)
library(ggplot2)

out_dir  <- "../results"
data_dir <- "../data"

# ==============================================================================
# 1. Load and filter H82 MOBA-seq data
# ==============================================================================
message("Loading H82 MOBA-seq count data...")
counts_df_H <- load_counts(
  file.path(data_dir, "mobaH82_FilteredSampleInfo.csv"),
  remove.sgids = c("None", "sgDummy"),
  keep.distance = 0
)

# Adjust preTran cell numbers (10x count for this dataset)
counts_df_H <- counts_df_H %>%
  mutate(cell_num = ifelse(sample_id == "preTran", 10 * count, cell_num))

# Pre-transplant histogram
result <- plot_pretran_histogram(
  counts_df_H,
  percentile.threshold = 0.99,
  filter.data = TRUE,
  return.list = TRUE
)
ggsave(file.path(out_dir, "H82_pretran_histogram.pdf"),
       plot = result$plot, width = 8, height = 8)

# Filter samples: remove Spi sgids, keep only cas9 samples + preTran
counts_H_filt <- counts_df_H %>%
  filter(
    sgid != "None",
    !grepl("^Spi", sgid, ignore.case = TRUE),
    grepl("cas9", sample_id, ignore.case = TRUE) | sample_id == "preTran"
  ) %>%
  mutate(
    gene = case_when(
      grepl("safe", sgid, ignore.case = TRUE) ~ "Ctrl",
      TRUE ~ sub("^sg", "", sgid, ignore.case = TRUE)
    ),
    time_point = case_when(
      is.na(time_point) ~ NA_character_,
      grepl("1wk", time_point, ignore.case = TRUE)  ~ "1week",
      grepl("3wk", time_point, ignore.case = TRUE)  ~ "3weeks",
      grepl("SubQ", time_point, ignore.case = TRUE)  ~ "3weeks",
      TRUE ~ time_point
    ),
    mouse_genotype = case_when(
      sample_id != "preTran" ~ "TEST",
      TRUE ~ mouse_genotype
    )
  )

# Apply tissue-specific cell number thresholds
counts_H_filt_2 <- counts_H_filt %>%
  filter(
    (tissue == "Cells" & cell_num >= 2) |
      (tissue == "Liver" & cell_num >= 100) |
      (!tissue %in% c("preTran", "Liver"))
  )

# Add log2 cell number
counts_H_filt_2 <- counts_H_filt_2 %>%
  mutate(log2_cell_num = log2(cell_num))

# Subset to key sgRNAs for focused analyses
counts_H_filt_f <- counts_H_filt_2 %>%
  filter(sgid %in% c("sgSAFE4", "sgNF1B", "sgCXCR2", "sgCXCR4"))

# ==============================================================================
# 2. Figure 3i: Colony jitter plot (Liver, 3 weeks)
# ==============================================================================
message("Generating Fig 3i: Colony jitter plot...")
plot_cell_num_jitter(
  counts_H_filt_f,
  group.by = "sgid",
  color.by = "gene",
  time.points = "3weeks",
  tissues = "Liver",
  gene.order = c("Ctrl", "CXCR2", "CXCR4", "NF1B"),
  colors = c("#196aaf", "#F28285", "#A48f7d", "#78A891")
)
ggsave(file.path(out_dir, "Fig3i_colony_jitter.pdf"),
       width = 10, height = 8, dpi = 300)

# ==============================================================================
# 3. Compute dormancy cutoffs and metastatic statistics
# ==============================================================================
message("Computing dormancy cutoffs for H82...")
resultsH <- dormancy_cutoff_by_genotype_tissue_time_v4(
  counts_H_filt_2,
  tissues = c("Liver"),
  time.points = c("1week", "3weeks"),
  min.colonies = 1,
  bw.adjust = 1,
  split.by.tissue = FALSE
)

geno_cutoffsH <- resultsH$cutoffs_df %>%
  dplyr::select(mouse_genotype, dormant_peak, dormancy_cutoff)

message("Computing pre-transplant data...")
pretran_dataH <- generate_pretran_data(counts_H_filt_2)

message("Computing bootstrap metastatic statistics (sgid-level)...")
results_H82 <- calculate_metastatic_statistics(
  counts_H_filt_2,
  tag = "sgid",
  tissues = c("Liver", "Subq"),
  statistics = c("seeding", "burden", "size_percentile",
                 "peak_mode", "dormancy_cutoff"),
  pretran.cell.num = pretran_dataH$pretran_cell_num,
  pretran.cell.num.gene = pretran_dataH$pretran_cell_num_gene,
  min.colonies.seeding = 10,
  min.colonies.burden = 10,
  min.colonies.peakmode = 25,
  min.colonies.dormancy = 25,
  R = 1000,
  dormancy.cutoff = geno_cutoffsH
)
write.csv(results_H82,
          file.path(out_dir, "h82_stats_sgid.csv"), row.names = FALSE)

# ==============================================================================
# 4. Figures 3j-k: Burden and seeding dot plots (Liver, 3 weeks)
# ==============================================================================
genes <- c("NF1B", "CXCR2", "CXCR4", "Ctrl")

message("Generating Fig 3j: Burden dot plot...")
plot_comparison_dots(
  results_H82, genes,
  metric = "burden",
  value.type = "fc",
  plot.type = "bar",
  filter.genotype = "TEST",
  filter.tissue = "Liver",
  filter.timepoint = "3weeks"
)
ggsave(file.path(out_dir, "Fig3j_burden_dot.pdf"),
       width = 8, height = 8, dpi = 300)

message("Generating Fig 3k: Seeding dot plot...")
plot_comparison_dots(
  results_H82, genes,
  metric = "seeding",
  value.type = "fc",
  plot.type = "bar",
  filter.genotype = "TEST",
  filter.tissue = "Liver",
  filter.timepoint = "3weeks"
)
ggsave(file.path(out_dir, "Fig3k_seeding_dot.pdf"),
       width = 8, height = 8, dpi = 300)

# ==============================================================================
# 5. Figure 3l: Cross-dataset radar comparison (H82 vs RP48)
# ==============================================================================
message("Generating Fig 3l: Cross-dataset radar...")

# Normalize H82 data for radar
radar_ready_dataH <- normalize_for_radar(
  results_H82,
  genes = c("CXCR2"),
  tissues = c("Liver"),
  time.points = c("1week", "3weeks"),
  genotypes = c("TEST"),
  dormancy.metric = "relative_dormancy_escape",
  check.ctrl.mean = TRUE
)

# Load RP48 radar data from script 07 results
# (07_mobaseq_cxcr2.R must have run first)
results_gene_rp48 <- read.csv(file.path(out_dir, "moba_cxcr2_gene_statistics.csv"))
radar_ready_dataC <- normalize_for_radar(
  results_gene_rp48,
  genes = c("CXCR2"),
  tissues = c("Liver"),
  time.points = c("1week", "3weeks"),
  genotypes = c("NSG"),
  dormancy.metric = "relative_dormancy_escape",
  check.ctrl.mean = TRUE
)

# Cross-dataset radar
plot_radar_cross_dataset(
  data1 = radar_ready_dataH,
  data2 = radar_ready_dataC,
  gene = "CXCR2",
  label1 = "H82",
  label2 = "RP48",
  tissues = c("Liver"),
  time.points = c("1week", "3weeks"),
  dormancy.metric = "relative_dormancy_escape",
  exclude.metrics = c("SuperMet%")
)
ggsave(file.path(out_dir, "Fig3l_radar_cross_dataset.pdf"),
       width = 12, height = 8, dpi = 300)

message("08_mobaseq_h82.R complete.")
