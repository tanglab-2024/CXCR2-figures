# ==============================================================================
# Dockerfile for Code Ocean compute capsule
# Manuscript: CXCL-CXCR2 signaling drives cancer-endothelium interactions
#             in SCLC metastatic seeding (CD-26-0546)
# ==============================================================================

FROM registry.codeocean.com/codeocean/r-base:4.3.2-ubuntu22.04

# System dependencies for R packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libhdf5-dev \
    libgeos-dev \
    libglpk-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    python3 \
    python3-pip \
    && rm -rf /var/lib/apt/lists/*

# Install BiocManager and CRAN packages
RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org')"

RUN R -e " \
  install.packages(c( \
    'Seurat', 'SeuratObject', 'remotes', 'dplyr', 'tibble', 'stringr', \
    'ggplot2', 'cowplot', 'patchwork', 'ggrepel', 'ggalluvial', \
    'circlize', 'scales', 'reticulate', 'pheatmap', 'ggpubr', \
    'lazyeval', 'gtable' \
  ), repos='https://cloud.r-project.org') \
"

# Bioconductor packages
RUN R -e " \
  BiocManager::install(c( \
    'SingleR', 'SingleCellExperiment', 'celldex', \
    'escape', 'msigdbr', 'fgsea', \
    'presto', 'ComplexHeatmap', \
    'clusterProfiler', 'org.Hs.eg.db', \
    'glmGamPoi', 'MAST', 'GSVA', \
    'pcaMethods', 'DESeq2', 'apeglm', 'tximport' \
  ), ask=FALSE, update=FALSE) \
"

# SeuratWrappers from GitHub
RUN R -e "remotes::install_github('satijalab/seurat-wrappers')"

# mobaseq package (bundled locally in environment/mobaseq/)
# To prepare: download the package source and place in environment/mobaseq/
COPY mobaseq/ /tmp/mobaseq/
RUN R -e "install.packages('/tmp/mobaseq', repos=NULL, type='source')"
