P1=readRDS("P1_seurat/P1_seurat.rds")
library(Seurat)
##########P1降维聚类
P1 <- subset(P1, subset = nFeature_Xenium >= 10)
P1 <- SCTransform(
  P1,
  assay = "Xenium",
  verbose = FALSE,
  variable.features.n = 422
)

P1 <- RunPCA(
  P1,
  assay = "SCT",
  features = rownames(P1),
  npcs = 30,
  verbose = FALSE
)

P1 <- FindNeighbors(
  P1,
  dims = 1:30
)

P1 <- FindClusters(
  P1,
  resolution = 0.3
)

P1 <- RunUMAP(
  P1,
  dims = 1:30
)

qs::qsave(P1,"P1_seurat/P1_process.qs")