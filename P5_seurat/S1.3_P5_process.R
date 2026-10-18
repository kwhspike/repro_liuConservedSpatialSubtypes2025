P5=readRDS("P5_seurat/P5_seurat.rds")
library(Seurat)
##########P5降维聚类
P5 <- subset(P5, subset = nFeature_Xenium >= 10)
P5 <- SCTransform(
  P5,
  assay = "Xenium",
  verbose = FALSE,
  variable.features.n = 422
)

P5 <- RunPCA(
  P5,
  assay = "SCT",
  features = rownames(P5),
  npcs = 30,
  verbose = FALSE
)

P5 <- FindNeighbors(
  P5,
  dims = 1:30
)

P5 <- FindClusters(
  P5,
  resolution = 0.3
)

P5 <- RunUMAP(
  P5,
  dims = 1:30
)

qs::qsave(P5,"P5_seurat/P5_process.qs")