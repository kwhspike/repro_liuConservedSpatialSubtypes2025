# Select epithelial clusters and split each sample into epithelial /
# non-epithelial objects, based on the marker dot plot.
# Run from the repository root (R >= 4.4, Seurat v5).
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

source("R/crc_marker_list.R")

P2=qs::qread("P2_seurat/P2_process.qs")

P2 <- FindClusters(
  P2,
  resolution = 0.4
)


# 只保留 panel 中实际存在的 genes
crc_marker_list <- lapply(
  crc_marker_list,
  function(x) x[x %in% rownames(P2)]
)

DefaultAssay(P2) <- "SCT"

Idents(P2) <- P2$SCT_snn_res.0.4

markers_P2 <- FindAllMarkers(
  P2,
  only.pos = TRUE,
  min.pct = 0.05,
  logfc.threshold = 0.1
)

toP20_P2 <- markers_P2 %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 10
  ) %>%
  ungroup()


DotPlot(
  P2,
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

epi_clusters <- c("0", "1", "2", "5", "6", "13", "14")
epi <- subset(
  P2,
  idents = epi_clusters
)
qs::qsave(epi,"P2_seurat/P2_epithelial.rds")

non_epi <- subset(
  P2,
  idents = setdiff(levels(Idents(P2)), epi_clusters)
)

qs::qsave(non_epi,"P2_seurat/P2_non_epithelial.rds")

pdf("P2_seurat/P2_umap_all.pdf")
DimPlot(P2,label = T)
dev.off()
