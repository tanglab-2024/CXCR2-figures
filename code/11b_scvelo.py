# ==============================================================================
# 11b_scvelo.py
# RNA velocity analysis using scVelo
#
# Input:  /data/Loomfiles/sample_alignments_B*.loom — velocyto loom files
#         /results/scvelo/cellID_obs.csv (from 11a)
#         /results/scvelo/cell_embeddings.csv (from 11a)
#         /results/scvelo/clusters.csv (from 11a)
# Output: Figure 1d — RNA velocity stream plots (monoculture, coculture-1, coculture-2)
#         as both PDF and SVG
#
# Run 11a_export_for_scvelo.R first to generate the Seurat exports.
# ==============================================================================

import matplotlib
matplotlib.use('Agg')
import numpy as np
import anndata
import scanpy
import scvelo as scv
import pandas as pd
import os

out_dir = "/results"
scvelo_dir = "/results/scvelo"
data_dir = "/data/Loomfiles"

# Cluster color palette (matches Seurat UMAP colors for EC1-EC5)
palette = ["#F8766D", "#A3A500", "#00BF7D", "#00B0F6", "#E76BF3"]

# Set scVelo figure parameters
scv.set_figure_params(
    style='scvelo', dpi=100, dpi_save=600, frameon=None,
    vector_friendly=True, transparent=True, fontsize=12,
    figsize=(10, 7), format='pdf'
)

# --- 1. Import Seurat metadata ------------------------------------------------
print("Loading Seurat metadata...")
sample_obs = pd.read_csv(os.path.join(scvelo_dir, "cellID_obs.csv"))
umap_cord = pd.read_csv(os.path.join(scvelo_dir, "cell_embeddings.csv"))
cell_clusters = pd.read_csv(os.path.join(scvelo_dir, "clusters.csv"))

# --- 2. Import and rename loom files ------------------------------------------
print("Loading loom files...")
B1 = anndata.read_loom(os.path.join(data_dir, "sample_alignments_B1.loom"))
B2 = anndata.read_loom(os.path.join(data_dir, "sample_alignments_B2.loom"))
B3 = anndata.read_loom(os.path.join(data_dir, "sample_alignments_B3.loom"))
B4 = anndata.read_loom(os.path.join(data_dir, "sample_alignments_B4.loom"))
B5 = anndata.read_loom(os.path.join(data_dir, "sample_alignments_B5.loom"))

for b in [B1, B2, B3, B4, B5]:
    b.var_names_make_unique()

# Rename loom cell IDs to match Seurat barcode format (B1_, B2_, etc.)
loom_prefixes = {
    B1: ('sample_alignments_2O1JF:', 'B1_'),
    B2: ('sample_alignments_6RIZA:', 'B2_'),
    B3: ('sample_alignments_OFUKD:', 'B3_'),
    B4: ('sample_alignments_6PF4O:', 'B4_'),
    B5: ('sample_alignments_QZ0JK:', 'B5_'),
}
for bdata, (old_prefix, new_prefix) in loom_prefixes.items():
    bdata.obs = bdata.obs.rename(index=lambda x: x.replace(old_prefix, new_prefix).replace('x', ''))

# Clean Seurat cell IDs to match
sample_obs.x = sample_obs.x.replace({"-1": ""}, regex=True)

# --- 3. Filter loom cells to match Seurat subset -----------------------------
print("Filtering loom cells to Seurat subset...")
B1_filt = B1[np.isin(B1.obs.index, sample_obs["x"])]
B2_filt = B2[np.isin(B2.obs.index, sample_obs["x"])]
B3_filt = B3[np.isin(B3.obs.index, sample_obs["x"])]
B4_filt = B4[np.isin(B4.obs.index, sample_obs["x"])]
B5_filt = B5[np.isin(B5.obs.index, sample_obs["x"])]

total_loom = sum(len(b.obs.index) for b in [B1_filt, B2_filt, B3_filt, B4_filt, B5_filt])
print(f"Seurat cells: {len(sample_obs)}, Loom cells: {total_loom}")


# --- Helper function: add UMAP + clusters to an AnnData object ----------------
def add_umap_and_clusters(adata, umap_cord, sample_obs, cell_clusters):
    """Add Seurat UMAP coordinates and cluster labels to an AnnData object."""
    idx = pd.DataFrame(adata.obs.index, columns=['CellID'])
    umap = umap_cord.rename(columns={'Unnamed: 0': 'CellID'})
    umap.CellID = umap.CellID.replace({"-1": ""}, regex=True)
    umap_ordered = idx.merge(umap, on="CellID").iloc[:, 1:]
    adata.obsm['X_umap'] = umap_ordered.values

    cluster_mapping = pd.Series(cell_clusters['x'].values, index=sample_obs['x'].values)
    adata.obs["clusters"] = adata.obs.index.map(cluster_mapping).astype("category")
    return adata


# --- Helper function: run scVelo and save stream plot -------------------------
def run_scvelo(adata, name, umap_cord, sample_obs, cell_clusters, palette, out_dir):
    """Run scVelo pipeline and save velocity stream plot as PDF and SVG."""
    adata = add_umap_and_clusters(adata, umap_cord, sample_obs, cell_clusters)

    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
    scv.tl.velocity(adata, mode="stochastic")
    scv.tl.velocity_graph(adata)

    # Set scVelo figure directory to /results
    scv.settings.figdir = out_dir

    # Save as PDF
    scv.pl.velocity_embedding_stream(
        adata, basis='umap', color="clusters",
        save=f'Fig1d_stream_{name}.pdf', palette=palette
    )
    # Save as SVG
    scv.settings.set_figure_params(format='svg')
    scv.pl.velocity_embedding_stream(
        adata, basis='umap', color="clusters",
        save=f'Fig1d_stream_{name}.svg', palette=palette
    )
    scv.settings.set_figure_params(format='pdf')

    return adata


# --- 4. Run scVelo on each sample group ---------------------------------------

# All samples merged
print("Running scVelo on all samples...")
B_All = B1_filt.concatenate(B2_filt, B3_filt, B4_filt, B5_filt, index_unique=None)
B_All = run_scvelo(B_All, "all", umap_cord, sample_obs, cell_clusters, palette, out_dir)

# Monoculture (B1 only) — Figure 1d
print("Running scVelo on monoculture (B1)...")
B1_filt = run_scvelo(B1_filt, "monoculture", umap_cord, sample_obs, cell_clusters, palette, out_dir)

# Coculture-1 (B2 + B3) — Figure 1d
print("Running scVelo on coculture-1 (B2+B3)...")
B_Co1 = B2_filt.concatenate(B3_filt, index_unique=None)
B_Co1 = run_scvelo(B_Co1, "coculture1", umap_cord, sample_obs, cell_clusters, palette, out_dir)

# Coculture-2 (B4 + B5) — Figure 1d
print("Running scVelo on coculture-2 (B4+B5)...")
B_Co2 = B4_filt.concatenate(B5_filt, index_unique=None)
B_Co2 = run_scvelo(B_Co2, "coculture2", umap_cord, sample_obs, cell_clusters, palette, out_dir)

print("11b_scvelo.py complete.")
