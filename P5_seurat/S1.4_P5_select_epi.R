# Select epithelial clusters and split each sample into epithelial /
# non-epithelial objects, based on the marker dot plot.
# Run from the repository root (R >= 4.4, Seurat v5).
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

source("R/crc_marker_list.R")

P5=qs::qread("P5_seurat/P5_process.qs")


# 只保留 panel 中实际存在的 genes
crc_marker_list <- lapply(
  crc_marker_list,
  function(x) x[x %in% rownames(P5)]
)

DefaultAssay(P5) <- "SCT"

Idents(P5) <- P5$SCT_snn_res.0.3

markers_P5 <- FindAllMarkers(
  P5,
  only.pos = TRUE,
  min.pct = 0.05,
  logfc.threshold = 0.1
)

toP20_P5 <- markers_P5 %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 10
  ) %>%
  ungroup()


DotPlot(
  P5,
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

epi_clusters <- c("0", "2", "5", "8", "9", "16", "17")
epi <- subset(
  P5,
  idents = epi_clusters
)
qs::qsave(epi,"P5_seurat/P5_epithelial.rds")

non_epi <- subset(
  P5,
  idents = setdiff(levels(Idents(P5)), epi_clusters)
)

qs::qsave(non_epi,"P5_seurat/P5_non_epithelial.rds")

pdf("P5_seurat/P5_umap_all.pdf")
DimPlot(P5,label = T)
dev.off()

