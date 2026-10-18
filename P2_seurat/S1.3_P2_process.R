P2=readRDS("P2_seurat/P2_seurat.rds")
library(Seurat)
##########P2降维聚类
P2 <- subset(P2, subset = nFeature_Xenium >= 10)
library(future)

options(future.globals.maxSize = 2 * 1024^3)
P2 <- SCTransform(
  P2,
  assay = "Xenium",
  verbose = FALSE,
  variable.features.n = 422
)

P2 <- RunPCA(
  P2,
  assay = "SCT",
  features = rownames(P2),
  npcs = 30,
  verbose = FALSE
)

P2 <- FindNeighbors(
  P2,
  dims = 1:30
)

P2 <- FindClusters(
  P2,
  resolution = 0.4   # downstream steps (S1.5 / S1.6 / S1.7 / S4.1) use SCT_snn_res.0.4
)

P2 <- RunUMAP(
  P2,
  dims = 1:30
)

qs::qsave(P2,"P2_seurat/P2_process.qs")