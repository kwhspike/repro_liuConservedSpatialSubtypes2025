# Build Seurat v5 object from exported Xenium P5 10x-style files
# and run basic QC -> normalization -> PCA -> UMAP -> clustering.

library(Seurat)
library(data.table)

# Run from the repository root. Intermediate files live inside the repository;
# only the raw Xenium zarr stores are read from XENIUM_BASE (see S1.2_export).

tenx_dir <- "P5_seurat/tenx"
out_rds  <- "P5_seurat/P5_seurat.rds"
out_meta <- "P5_seurat/P5_seurat_meta.csv"

# ---- read matrix -----------------------------------------------------------
counts <- Read10X(data.dir = tenx_dir, gene.column = 2, unique.features = TRUE)
# counts: genes x cells (dgCMatrix)
stopifnot(all(!duplicated(rownames(counts))))

obj <- CreateSeuratObject(
  counts = counts,
  project = "Xenium_P5",
  assay = "Xenium",
  min.cells = 0,
  min.features = 0
)

# ---- add metadata ----------------------------------------------------------
meta <- fread(file.path(tenx_dir, "metadata.csv.gz"), data.table = FALSE)
rownames(meta) <- meta$cell_id
meta <- meta[colnames(obj), , drop = FALSE]
stopifnot(all(rownames(meta) == colnames(obj)))

obj <- AddMetaData(obj, metadata = meta)
print(head(obj@meta.data))

# ---- basic processing ------------------------------------------------------
obj <- NormalizeData(obj, normalization.method = "LogNormalize", verbose = FALSE)
obj <- FindVariableFeatures(obj, nfeatures = 2000, verbose = FALSE)
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
obj <- FindNeighbors(obj, dims = 1:30, verbose = FALSE)
obj <- FindClusters(obj, resolution = 0.5, verbose = FALSE)
obj <- RunUMAP(obj, dims = 1:30, verbose = FALSE)

print(obj)
print(table(Idents(obj)))

# ---- save ------------------------------------------------------------------
saveRDS(obj, out_rds)

# export metadata (cell_id + cluster + umap + coordinates) for python mapping
umap_emb <- Embeddings(obj, "umap")
export <- data.frame(
  cell_id = colnames(obj),
  cluster = as.character(Idents(obj)),
  umap1 = umap_emb[, 1],
  umap2 = umap_emb[, 2],
  x_centroid = obj$x_centroid,
  y_centroid = obj$y_centroid,
  stringsAsFactors = FALSE
)
fwrite(export, out_meta)

cat("saved:", out_rds, "\n")
cat("saved:", out_meta, "\n")
