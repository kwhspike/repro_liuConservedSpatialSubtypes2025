# Export epithelial SCT_snn_res.0.3 cluster labels (with cell_id) to csv,
# for writing back into the Xenium zarr.
# Run from the repository root.
suppressPackageStartupMessages(library(qs))

obj <- qread("P1_seurat/P1_epithelial.rds")

cluster_col <- "SCT_snn_res.0.3"
stopifnot(cluster_col %in% colnames(obj@meta.data))

df <- data.frame(
  cell_id = obj@meta.data$cell_id,
  SCT_snn_res_0.3 = obj@meta.data[[cluster_col]],
  row.names = NULL,
  stringsAsFactors = FALSE
)
# drop any cells missing cell_id
df <- df[!is.na(df$cell_id) & df$cell_id != "", ]

out <- "P1_seurat/P1_epi_SCT_snn_res_0.3.csv"
write.csv(df, out, row.names = FALSE, quote = TRUE)
cat("wrote", nrow(df), "cells to", out, "\n")
print(table(df$SCT_snn_res_0.3))
