# Select epithelial clusters and split each sample into epithelial /
# non-epithelial objects, based on the marker dot plot.
# Run from the repository root (R >= 4.4, Seurat v5).
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

source("R/crc_marker_list.R")

P1=qs::qread("P1_seurat/P1_process.qs")

# 只保留 panel 中实际存在的 genes
crc_marker_list <- lapply(
  crc_marker_list,
  function(x) x[x %in% rownames(P1)]
)

DefaultAssay(P1) <- "SCT"

Idents(P1) <- P1$SCT_snn_res.0.3

markers_P1 <- FindAllMarkers(
  P1,
  only.pos = TRUE,
  min.pct = 0.05,
  logfc.threshold = 0.1
)

top10_P1 <- markers_P1 %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 10
  ) %>%
  ungroup()


DotPlot(
  P1,
  features = crc_marker_list,
  assay = "SCT",
  dot.scale = 2.5
) +
  scale_color_gradientn(
    colors = c(
      "#2166AC",
      "#1A9850",
      "#FEE08B",
      "#D73027"
    )
  ) +
  RotatedAxis()

epi_clusters <- c("0", "1", "4", "11", "12", "13", "14", "15", "16")
epi <- subset(
  P1,
  idents = epi_clusters
)
qs::qsave(epi,"P1_seurat/P1_epithelial.rds")

non_epi <- subset(
  P1,
  idents = setdiff(levels(Idents(P1)), epi_clusters)
)

qs::qsave(non_epi,"P1_seurat/P1_non_epithelial.rds")

pdf("P1_seurat/P1_umap_all.pdf")
DimPlot(P1,label = T)
dev.off()
