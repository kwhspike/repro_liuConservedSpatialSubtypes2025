# R package dependencies for this pipeline.
# Tested on R 4.4.1 (Ubuntu 20.04), see sessionInfo_output.txt.
#
# Seurat v5 is required: the SCTransform / assay handling in P*_seurat/S1.3
# uses v5 syntax.
#
# Usage:
#   Rscript install.R

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

install.packages(c(
  "Seurat",      # v5: single-sample processing + integration + Leiden (S4.4)
  "data.table",
  "qs",
  "dbscan",      # frNN fixed-radius neighborhood search (NMF/S4.3*)
  "pheatmap",
  "ggplot2",
  "future",
  "doParallel",  # parallel nmfEstimateRank (optional, falls back to serial)
  "plyr",
  "dplyr",
  "reticulate"   # Seurat Leiden (algorithm = 4) calls python leidenalg
))

# S4.3 was run with NMF 0.27 (NMF >= 0.27 + R >= 4.2 has a verbose-output bug,
# see the note in NMF/S4.3_CAF_NMF.R). Install this exact version:
remotes::install_version("NMF", version = "0.27")

# The Seurat Leiden clustering in NMF/S4.4_CAF_Leiden_clustering.R needs
# python package leidenalg inside the environment reticulate binds to:
#   export RETICULATE_PYTHON=<conda env>/bin/python   # e.g. the squidpy_env
# (leidenalg is already included in environment.yml)
