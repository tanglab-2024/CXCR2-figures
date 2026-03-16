# ==============================================================================
# 10_mobaseq_cxcr_inhibitor.R
# MOBA-seq analysis of CXCR2 inhibitor (AZD) treatment effect on metastasis
#
# Input:  data/CXCR_inhibitor_FilteredSampleInfo.txt — inhibitor MOBA-seq data
# Output: Figure 6l — Colony jitter + tumor histogram (Control vs AZD)
#         Figure 6m — Boxplot of metastatic burden by treatment
#         Figure 6n — Boxplot of metastatic seeding by treatment
#
# Requires: mobaseq package (bundled in environment/mobaseq/)
# ==============================================================================

library(mobaseq)
library(dplyr)
library(stringr)
library(ggplot2)

out_dir  <- "../results"
data_dir <- "../data"

# ==============================================================================
# 1. Load and annotate CXCR inhibitor MOBA-seq data
# ==============================================================================
message("Loading CXCR inhibitor MOBA-seq data...")
counts_df <- load_counts(file.path(data_dir, "CXCR_inhibitor_FilteredSampleInfo.txt"))

# Assign treatment groups based on mouse ID
counts_df <- counts_df %>%
  mutate(
    treatment = case_when(
      mouse_id %in% c(73149, 73150, 73152, 73153, 73155) ~ "Control",
      mouse_id %in% c(73148, 73143, 73145, 73142, 73141) ~ "AZD(2mg/ml)",
      TRUE ~ NA_character_
    )
  )

# Filter to sgDummy barcodes with tissue-specific thresholds
subset_df <- counts_df %>%
  filter(
    sgid == "sgDummy",
    (str_detect(tolower(sample_id), "pretran") & cell_num >= 2) |
      (tissue == "Liver" & cell_num >= 10000)
  ) %>%
  mutate(
    gene = case_when(
      sgid == "sgDummy" ~ "Dummy",
      grepl("sgNT|sgsafe", sgid, ignore.case = TRUE) ~ "Ctrl",
      grepl("mmu-sg[0-9]+", sgid) ~ sub("mmu-sg[0-9]+", "", sgid),
      TRUE ~ NA_character_
    )
  )

# ==============================================================================
# 2. Figure 6l: Colony jitter + tumor histogram
# ==============================================================================
message("Generating Fig 6l: Colony jitter and histogram...")

# Jitter plot
p_jitter <- plot_cell_num_jitter(
  subset_df,
  exclude.samples = c("preTran"),
  group.by = "treatment",
  color.by = "treatment",
  alpha = 0.2
)
ggsave(file.path(out_dir, "Fig6l_seeding_jitter.pdf"),
       plot = p_jitter, width = 9, height = 8, dpi = 300)

# Tumor histogram
p_hist <- plot_tumor_histogram(subset_df, group.by = "treatment")
ggsave(file.path(out_dir, "Fig6l_seeding_hist.pdf"),
       plot = p_hist, width = 6, height = 8, dpi = 300)

# ==============================================================================
# 3. Figure 6m: Metastatic burden boxplot
# ==============================================================================
message("Generating Fig 6m: Burden boxplot...")
p_burden <- plot_boxplot_stats(
  subset_df,
  group.by = "treatment",
  metric = "total",
  label.points = TRUE
)
ggsave(file.path(out_dir, "Fig6m_burden_boxplot.pdf"),
       plot = p_burden, width = 8, height = 8, dpi = 300)

# ==============================================================================
# 4. Figure 6n: Metastatic seeding boxplot
# ==============================================================================
message("Generating Fig 6n: Seeding boxplot...")
p_seeding <- plot_boxplot_stats(
  subset_df,
  group.by = "treatment",
  metric = "count",
  label.points = TRUE
)
ggsave(file.path(out_dir, "Fig6n_seeding_boxplot.pdf"),
       plot = p_seeding, width = 8, height = 8, dpi = 300)

message("10_mobaseq_cxcr_inhibitor.R complete.")
