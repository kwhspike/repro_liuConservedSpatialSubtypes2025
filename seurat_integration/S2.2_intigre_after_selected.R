# Integrate the non-epithelial cells of all three samples (Seurat rpca),
# cluster, and annotate major cell types.
# Run from the repository root (R >= 4.4, Seurat v5).
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

source("R/crc_marker_list.R")

P1=qs::qread("P1_seurat/P1_non_epithelial.rds")
P2=qs::qread("P2_seurat/P2_non_epithelial.rds")
P5=qs::qread("P5_seurat/P5_non_epithelial.rds")

obj.list <- list(
  P1 = P1,
  P2 = P2,
  P5 = P5
)

rm(P1, P2, P5)
gc()

# ------------------------------ normalize -----------------------------------
obj.list <- lapply(obj.list, function(x) {
  DefaultAssay(x) <- "Xenium"
  x <- NormalizeData(x, verbose = FALSE)
  return(x)
})
# ------------------------------ integration features -------------------------
features <- SelectIntegrationFeatures(
  object.list = obj.list,
  nfeatures = 1000
)
# ------------------------------ PCA ------------------------------------------
obj.list <- lapply(obj.list, function(x) {
  x <- ScaleData(
    x,
    features = features,
    verbose = FALSE
  )
  x <- RunPCA(
    x,
    features = features,
    npcs = 30,
    verbose = FALSE
  )
  return(x)
})
# ------------------------------ anchors --------------------------------------
anchors <- FindIntegrationAnchors(
  object.list = obj.list,
  anchor.features = features,
  reduction = "rpca",
  dims = 1:30
)
# ------------------------------ integrate ------------------------------------
obj <- IntegrateData(
  anchorset = anchors,
  dims = 1:30
)

qs::qsave(obj,"seurat_integration/integra_non_epi.qs")

####################
DefaultAssay(obj) <- "integrated"
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(
  obj,
  npcs = 30,
  verbose = FALSE
)
obj <- RunUMAP(
  obj,
  reduction = "pca",
  dims = 1:30
)

obj <- FindNeighbors(
  obj,
  reduction = "pca",
  dims = 1:30
)

obj <- FindClusters(
  obj,
  resolution = 0.5
)
##################现在开始分群
pdf("seurat_integration/all_sample_dimplot.pdf")
DimPlot(obj,group.by = "orig.ident")
dev.off()
pdf("seurat_integration/all_cluster_dimplot.pdf")
DimPlot(obj,label = T)
dev.off()

Idents(obj) <- obj$integrated_snn_res.0.5 

markers_obj <- FindAllMarkers(
  obj,
  only.pos = TRUE,
  min.pct = 0.05,
  logfc.threshold = 0.1
)
qs::qsave(markers_obj,"seurat_integration/markers_nonepi.qs")
top10_obj <- markers_obj %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 10
  ) %>%
  ungroup()
print(top10_obj,n=300)

#################发现有上皮混杂，继续细分
DotPlot(
  obj,
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
##################注释细胞
# 根据 cluster 进行大类细胞类型注释
celltype <- c(
  "0"  = "T cell",
  "1"  = "Fibroblast",
  "2"  = "CAF cell",
  "3"  = "Macrophage",
  "4"  = "Plasma cell",
  "5"  = "T cell",
  "6"  = "Smooth muscle cell",
  "7"  = "B cell",
  "8"  = "Tumor cell",
  "9"  = "Monocyte",
  "10" = "Epithelial cell",
  "11" = "Pericyte",
  "12" = "Tumor cell",
  "13" = "Mast cell",
  "14" = "Fibroblast",
  "15" = "Neural/Schwann cell",
  "16" = "Dendritic cell",
  "17" = "Cycling cell",
  "18" = "Tumor cell",
  "19" = "CAF cell"
)
# 添加到 Seurat 对象
obj$celltype <- plyr::mapvalues(
  x = obj$seurat_clusters,
  from = names(celltype),
  to = celltype
)

# 查看每个 cluster 对应的细胞类型
table(obj$celltype)
table(obj$seurat_clusters, obj$celltype)
cells_keep <- colnames(obj)[
  !as.character(obj$celltype) %in% c("20", "21", "22")
]

obj <- subset(obj, cells = cells_keep)
obj$celltype <- droplevels(obj$celltype)
qs::qsave(obj,"seurat_integration/integra_non_epi_process.qs")
####################
celltype <- c(
  "8"  = "8",
  "10" = "10",
  "12" = "12",
  "18" = "18"
)
obj$celltype <- plyr::mapvalues(
  x = obj$seurat_clusters,
  from = names(celltype),
  to = celltype
)
table(obj$celltype)
table(obj$seurat_clusters, obj$celltype)
cells_keep <- colnames(obj)[
  as.character(obj$celltype) %in% c("12", "18", "8","10")
]

obj <- subset(obj, cells = cells_keep)
obj$celltype <- droplevels(obj$celltype)
qs::qsave(obj,"seurat_integration/integra_non_epi_epi.qs")
