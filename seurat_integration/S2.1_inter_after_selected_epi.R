library(Seurat)
P1=qs::qread("P1_seurat/P1_epithelial.rds")
P2=qs::qread("P2_seurat/P2_epithelial.rds")
P5=qs::qread("P5_seurat/P5_epithelial.rds")

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
DefaultAssay(obj) <- "integrated"

qs::qsave(obj,"seurat_integration/integra_epi.qs")