# Export non-epithelial integration celltype annotations (with cell_id) to
# per-sample csv files, for writing back into the Xenium zarrs.
suppressPackageStartupMessages(library(qs))

obj <- qread("seurat_integration/integra_non_epi_process.qs")

stopifnot("celltype" %in% colnames(obj@meta.data))
stopifnot("cell_id" %in% colnames(obj@meta.data))

md <- obj@meta.data
# drop any cells missing cell_id
md <- md[!is.na(md$cell_id) & md$cell_id != "", ]

for (sample in c("Xenium_P1", "Xenium_P2", "Xenium_P5")) {
  sub <- md[md$orig.ident == sample, c("cell_id", "celltype")]
  short <- sub("Xenium_", "", sample)  # P1 / P2 / P5
  out <- paste0("seurat_integration/", short, "_non_epi_celltype.csv")
  write.csv(sub, out, row.names = FALSE, quote = TRUE)
  cat("wrote", nrow(sub), "cells to", out, "\n")
  print(table(sub$celltype))
}
