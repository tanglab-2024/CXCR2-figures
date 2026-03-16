# ==============================================================================
# 09_mobaseq_rac1.R
# MOBA-seq analysis of Rac1 knockout metastatic dynamics
#
# Input:  data/MobaV_FilteredSampleInfo.csv — Rac1 MOBA-seq barcode count data
# Output: Figure 5n — Composition bar plot (NSG, Liver, 3 weeks)
#         Figure 5o — Colony jitter + histogram (NSG, Liver, 3 weeks)
#         Figure 5p — Bootstrap seeding bar plot (NSG, Liver)
#         Figure 5q — Density mode plot (NSG, Liver)
#
# Requires: mobaseq package (bundled in environment/mobaseq/)
# ==============================================================================

library(mobaseq)
library(dplyr)
library(stringr)
library(ggplot2)

out_dir  <- "/results"
data_dir <- "/data"

# ==============================================================================
# 1. Load and filter Rac1 MOBA-seq data
# ==============================================================================
message("Loading Rac1 MOBA-seq count data...")
counts_df_rac1 <- load_counts(
  file.path(data_dir, "MobaV_FilteredSampleInfo.csv"),
  remove.sgids = "sgDummy",
  keep.bc.end = c("CTGA", "GGAA"),
  keep.distance = 0
)

# Recode tissue names
counts_df_rac1 <- counts_df_rac1 %>%
  mutate(tissue = recode(tissue, "Cells" = "preTran", "WholeBlood" = "Blood"))

# Pre-transplant histogram and filtering
result_rac1 <- plot_pretran_histogram(
  counts_df_rac1,
  percentile.threshold = 0.99,
  filter.data = TRUE,
  return.list = TRUE
)
ggsave(file.path(out_dir, "Rac1_pretran_histogram.pdf"),
       plot = result_rac1$plot, width = 8, height = 8)
counts_df_rac1_filt <- result_rac1$filtered_df

# Filter to Rac1/control sgRNAs and apply tissue-specific thresholds
counts_df_rac1_filt2 <- counts_df_rac1_filt %>%
  filter(
    sgid != "sgDummy",
    distance == 0,
    str_detect(sgid, regex("Rac1|NT2|NT3|Safe4", ignore_case = TRUE)),
    (tissue == "preTran" & cell_num >= 2) |
      (tissue == "WholeBlood" & cell_num >= 1) |
      (tissue == "Brain" & cell_num >= 10) |
      (tissue == "BoneMarrow" & cell_num >= 5) |
      (!tissue %in% c("preTran", "Blood", "Brain", "BoneMarrow") & cell_num >= 100)
  ) %>%
  mutate(
    gene = case_when(
      str_detect(sgid, regex("^sg(nt|safe)", ignore_case = TRUE)) ~ "Ctrl",
      str_detect(sgid, regex("mmu-sg[0-9]+", ignore_case = TRUE)) ~
        str_replace(sgid, regex("mmu-sg[0-9]+", ignore_case = TRUE), ""),
      TRUE ~ NA_character_
    ),
    gene = str_replace(gene, "-$", "")
  )

# ==============================================================================
# 2. Figure 5n: Composition bar plot (NSG, Liver, 3 weeks)
# ==============================================================================
message("Generating Fig 5n: Composition bar plot...")
p_comp <- plot_composition(
  counts_df_rac1_filt2 %>%
    filter(tissue %in% c("Liver", "preTran")) %>%
    filter(time_point == "3weeks" | is.na(time_point)) %>%
    filter(mouse_genotype == "NSG" | is.na(mouse_genotype)),
  stack.by = "sgid",
  color.by = "sgid",
  main.gene = "Rac1",
  show.labels = TRUE,
  label.threshold = 1,
  pretran.right = TRUE
)
ggsave(file.path(out_dir, "Fig5n_composition.pdf"),
       plot = p_comp, width = 9, height = 9, dpi = 300)

# ==============================================================================
# 3. Figure 5o: Colony jitter + tumor histogram (NSG, Liver, 3 weeks)
# ==============================================================================
message("Generating Fig 5o: Colony jitter and histogram...")

# Jitter plot
p_jitter <- plot_cell_num_jitter(
  counts_df_rac1_filt2,
  group.by = "sgid",
  color.by = "gene",
  tissues = "Liver",
  mouse.genotypes = "NSG",
  exclude.samples = c("preTran"),
  time.points = "3weeks"
)
ggsave(file.path(out_dir, "Fig5o_colony_jitter.pdf"),
       plot = p_jitter, width = 9, height = 9, dpi = 300)

# Tumor histogram
p_hist <- plot_tumor_histogram(
  counts_df_rac1_filt2 %>%
    filter(tissue == "Liver", mouse_genotype == "NSG"),
  bins = 50,
  group.by = "gene",
  group.values = c("Ctrl", "Rac1"),
  colors = c("Ctrl" = "#4472C4", "Rac1" = "#E8897E"),
  time.points = "3weeks"
)
ggsave(file.path(out_dir, "Fig5o_tumor_histogram.pdf"),
       plot = p_hist, width = 4, height = 8, dpi = 300)

# ==============================================================================
# 4. Compute dormancy cutoffs and metastatic statistics
# ==============================================================================
message("Computing pre-transplant data...")
pretran_data_rac1 <- generate_pretran_data(counts_df_rac1_filt2)

message("Computing dormancy cutoffs...")
resultsRac1 <- dormancy_cutoff_by_genotype_tissue_time_v4(
  counts_df_rac1_filt2,
  tissues = c("Brain", "Lung", "Liver", "BoneMarrow"),
  time.points = c("3weeks"),
  min.colonies = 1,
  bw.adjust = 1,
  split.by.tissue = FALSE
)

geno_cutoffsRac <- resultsRac1$cutoffs_df %>%
  dplyr::select(mouse_genotype, dormant_peak, dormancy_cutoff)

message("Computing bootstrap metastatic statistics (sgid-level)...")
results_rac1_sgid <- calculate_metastatic_statistics(
  counts_df_rac1_filt2,
  tag = "sgid",
  tissues = c("Brain", "Lung", "Liver", "BoneMarrow"),
  statistics = c("seeding", "burden", "size_percentile",
                 "supermet", "dormancy_cutoff", "peak_mode"),
  pretran.cell.num = pretran_data_rac1$pretran_cell_num,
  pretran.cell.num.gene = pretran_data_rac1$pretran_cell_num_gene,
  R = 1000,
  dormancy.cutoff = geno_cutoffsRac
)
write.csv(results_rac1_sgid,
          file.path(out_dir, "rac1_stats_sgid.csv"), row.names = FALSE)

message("Computing bootstrap metastatic statistics (gene-level)...")
results_rac1_gene <- calculate_metastatic_statistics(
  counts_df_rac1_filt2,
  tag = "gene",
  tissues = c("Brain", "Lung", "Liver", "BoneMarrow"),
  statistics = c("seeding", "burden", "size_percentile",
                 "supermet", "dormancy_cutoff", "peak_mode"),
  pretran.cell.num = pretran_data_rac1$pretran_cell_num,
  pretran.cell.num.gene = pretran_data_rac1$pretran_cell_num_gene,
  R = 1000,
  dormancy.cutoff = geno_cutoffsRac
)
write.csv(results_rac1_gene,
          file.path(out_dir, "rac1_stats_gene.csv"), row.names = FALSE)

# ==============================================================================
# 5. Figure 5p: Bootstrap seeding bar plot (NSG, Liver)
# ==============================================================================
message("Generating Fig 5p: Bootstrap seeding bar plot...")
plot_boot_met_seed(
  results_rac1_sgid,
  tissues = "Liver",
  mouse.genotypes = "NSG",
  time.points = c("1week", "3weeks"),
  facet.by = "time_point",
  treatment.name = "Rac1"
)
ggsave(file.path(out_dir, "Fig5p_boot_seeding.pdf"),
       width = 10, height = 8, dpi = 300)

# ==============================================================================
# 6. Figure 5q: Density mode plot (NSG, Liver)
# ==============================================================================
message("Generating Fig 5q: Density mode plot...")
p_density <- plot_density_modes(
  counts_df_rac1_filt2 %>% filter(mouse_genotype == "NSG"),
  genes = c("Ctrl", "Rac1"),
  time.points = c("1week", "3weeks"),
  tissues = "Liver",
  xlim = c(1e2, 1e6),
  color.by = "gene"
)
ggsave(file.path(out_dir, "Fig5q_density_modes.pdf"),
       plot = p_density, width = 9, height = 6, dpi = 300)

message("09_mobaseq_rac1.R complete.")
