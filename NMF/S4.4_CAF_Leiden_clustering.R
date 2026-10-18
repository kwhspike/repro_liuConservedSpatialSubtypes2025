#!/usr/bin/env Rscript
# =============================================================================
# S4.4_CAF_Leiden_clustering.R
#
# S4.3b_CAF_NMF_k5.R 的后续步骤: 选定 k 后对 W 矩阵做 Leiden 聚类 -> 内皮空间 niche,
# 统计每个 niche 的平均【原始】邻域组成 (与原文 groupby('subtype').mean() 一致)。
# 本版本不做自动注释: 聚类结果 (05 表 + Fig3/Fig4) 供人工命名亚群。
# 注意: 原文方法不含逐样本归一化 (邻域组成 / 全切片组成), 故本脚本不做该归一化。
#
# 输入 (均由 S4.3b_CAF_NMF_k5.R 生成, 运行前请先跑 S4.3b 并选定 radius):
#   - k5_r{radius}/01_CAF_neighbor_comp_per_sample_r{radius}.qs  (每样本邻域组成 + 坐标)
#   - k5_r{radius}/03_NMF_W_cell_by_factor.qs  (W 矩阵)
#   - k5_r{radius}/04_NMF_H_factor_by_celltype.csv
#   - NMF/input/cell_metadata.csv          (仅用于补 sample_id / cell_type 列)
#
# 输出 (选定 k 的半径子目录内):
#   - 05_cluster_mean_neighbor_comp.csv          各 niche 平均邻域组成 (人工命名依据)
#   - 08_CAF_spatial_clusters_per_cell.csv      每个内皮细胞的坐标 + cluster
#   - figures/Fig1~Fig4, Fig5                    空间分布 / H 热图 / 组成热图 / 组成条形图 / UMAP
# 人工命名后, 如需在 08 表加 annotation 列, 自建一个 cluster -> 名称的 CSV 再 merge 即可。
#
# 运行:  Rscript NMF/S4.4_CAF_Leiden_clustering.R
#
# 注意: Seurat 的 Leiden (algorithm = 4) 依赖 reticulate + python 包 leidenalg。
#       运行前把 reticulate 指向带 leidenalg 的环境 (如 squidpy_env):
#         export RETICULATE_PYTHON=~/anaconda3/envs/squidpy_env/bin/python
#       Leiden resolution 在下方 config 中手动调整, 每改一次重跑本脚本即可。
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(qs)          # qread
  library(Seurat)      # FindClusters(algorithm = 4) = Leiden, 需 v4.3+
  library(ggplot2)
  library(pheatmap)
})

## ============================ 0. 配置 ======================================
## !! 选定邻域半径后, 只改下面 radius 一项, 路径自动生成; 每改一次重跑本脚本 !!
radius <- 80   # 邻域半径 (um): 20 / 40 / 60 / 80

config <- list(
  radius = radius,

  input_path  = "NMF/input/cell_metadata.csv",
  sample_col   = "sample_id",     # 与 S4.3b config 一致
  celltype_col = "cell_type",
  target_label = "CAF cell",  # 分析对象, 与 S4.3b target_labels 一致
  qs_per_sample = sprintf("NMF/output_CAF/k5_r%d/01_CAF_neighbor_comp_per_sample_r%d.qs",
                          radius, radius),
  qs_W        = sprintf("NMF/output_CAF/k5_r%d/03_NMF_W_cell_by_factor.qs", radius),
  csv_H       = sprintf("NMF/output_CAF/k5_r%d/04_NMF_H_factor_by_celltype.csv", radius),

  ## ---- Leiden 聚类参数 (手动调整这里) ----
  leiden_resolution = 0.05,  # Leiden resolution, 聚类数随其增大而增多

  ## ---- 输出: 与该 radius 的 NMF 结果同目录, 便于不同半径的结果对比 ----
  outdir = sprintf("NMF/output_CAF/k5_r%d", radius)
)

dir.create(config$outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(config$outdir, "figures"), showWarnings = FALSE)

## ======================= 1. 读取 S4.3 的结果 ==============================
per_sample <- qread(config$qs_per_sample)
W          <- qread(config$qs_W)
k_nmf      <- ncol(W)

## ---- 合并所有样本的【原始】邻域组成矩阵 (与原文一致, 不归一化) ----
## 注意: Xenium 的 cell_id 仅样本内唯一, 需构造全局唯一 id (sample|cell_id),
##       否则跨样本同 id 细胞在后续 merge 时会错配/重复
comp_all <- do.call(rbind, lapply(names(per_sample), function(s) {
  m <- per_sample[[s]]$raw
  rownames(m) <- paste(s, rownames(m), sep = "|")
  m
}))
stopifnot(nrow(W) == nrow(comp_all),
          all(rownames(W) == sub("^[^|]+\\|", "", rownames(comp_all))))
rownames(W) <- rownames(comp_all)   # 统一换成全局唯一 id
cat(sprintf("[INFO] W: %d 内皮细胞 x %d factor\n", nrow(W), k_nmf))

## ======================= 2. Leiden 聚类 (Seurat) ==========================
# 需要 reticulate 指向含 leidenalg 的 python 环境 (如 squidpy_env):
#   export RETICULATE_PYTHON=~/anaconda3/envs/squidpy_env/bin/python
# W 是 NMF 得分矩阵, 不做 NormalizeData (Seurat v5 需显式填充 data layer,
# 否则 ScaleData 报 "No layer matching pattern 'data' found")
W_t <- as(t(W), "dgCMatrix")
rownames(W_t) <- gsub("_", "-", rownames(W_t))   # Seurat 会把 feature 名中的 _ 替换为 -
CAF_seu <- CreateSeuratObject(counts = W_t, assay = "NMF")
CAF_seu <- SetAssayData(CAF_seu, layer = "data", new.data = W_t)
CAF_seu <- ScaleData(CAF_seu, features = rownames(CAF_seu), verbose = FALSE)
CAF_seu <- RunPCA(CAF_seu, features = rownames(CAF_seu),
                   npcs = k_nmf, approx = FALSE, verbose = FALSE)  # 仅 5 个特征, 用精确 SVD
CAF_seu <- FindNeighbors(CAF_seu, dims = 1:k_nmf, verbose = FALSE)
CAF_seu <- FindClusters(CAF_seu, resolution = config$leiden_resolution,
                         algorithm = 4, verbose = FALSE)   # 4 = Leiden
clusters <- as.character(Idents(CAF_seu))
names(clusters) <- colnames(CAF_seu)
cat(sprintf("[INFO] Leiden resolution = %g, 内皮空间 niche (clusters):\n",
            config$leiden_resolution))
print(table(clusters))

## ============== 3. 各 niche 的平均邻域组成 (供手动注释) ====================
# 每个 cluster: 对其成员内皮细胞的【原始】邻域组成取均值 (同原文 groupby().mean()),
# 按原始细胞类型直接输出, 不做大类合并
cl_names <- sort(unique(clusters))
cluster_profile <- t(vapply(cl_names, function(cl)
  colMeans(comp_all[names(clusters)[clusters == cl], , drop = FALSE]),
  numeric(ncol(comp_all))))
rownames(cluster_profile) <- cl_names

fwrite(as.data.table(cluster_profile, keep.rownames = "cluster"),
       file.path(config$outdir, "05_cluster_mean_neighbor_comp.csv"))
cat("[INFO] 各 niche 平均邻域组成已保存至 05_cluster_mean_neighbor_comp.csv,",
    "请结合 Fig3/Fig4 人工命名亚群。\n", sep = "\n")

## ======================= 4. 输出每个内皮细胞的结果表 ========================
# 不做自动注释; 人工命名后如需加 annotation 列, 见 S4.5_CAF_manual_annotation.R
cell_anno <- data.table(uid = names(clusters), cluster = clusters)
coords_all <- rbindlist(lapply(names(per_sample), function(s) {
  cd <- as.data.table(per_sample[[s]]$coords)
  cd[, `:=`(uid = paste(s, cell_id, sep = "|"), sample_id = s)]
  cd
}))
cell_anno <- merge(cell_anno, coords_all, by = "uid")
cell_anno[, cell_type := paste(config$target_label, collapse = ", ")]  # 全部为分析对象细胞
cell_anno <- cell_anno[, .(cell_id, cluster, x, y, sample_id, cell_type)]
fwrite(cell_anno, file.path(config$outdir, "08_CAF_spatial_clusters_per_cell.csv"))

## ============================ 5. 绘图 ======================================
# Fig1: 各样本空间分布, 按内皮空间 niche (Leiden cluster) 着色
p1 <- ggplot(cell_anno, aes(x = x, y = y, color = cluster)) +
  geom_point(size = 0.4) + coord_fixed() +
  facet_wrap(~ sample_id) +
  theme_bw() + labs(title = "CAF spatial niches")
ggsave(file.path(config$outdir, "figures", "Fig1_CAF_spatial_clusters.pdf"),
       p1, width = 12, height = 5)

# Fig2: NMF factor-by-celltype 热图 (H 按行归一化展示相对组成)
H <- as.matrix(fread(config$csv_H), rownames = 1)
H_show <- sweep(H, 1, rowSums(H), "/")
pdf(file.path(config$outdir, "figures", "Fig2_NMF_H_factor_by_celltype.pdf"),
    width = 8, height = 4)
pheatmap(H_show, cluster_rows = TRUE, cluster_cols = TRUE,
         main = "NMF factor-by-celltype (row-normalized)", fontsize = 8)
dev.off()

# Fig3: niche 平均邻域组成热图 (细胞类型分辨率)
pdf(file.path(config$outdir, "figures", "Fig3_cluster_neighbor_comp_heatmap.pdf"),
    width = 9, height = 5)
pheatmap(cluster_profile, scale = "column",
         main = "Mean normalized neighborhood composition per CAF niche",
         fontsize = 8)
dev.off()

# Fig4: niche 平均邻域组成堆叠条形图 (细胞类型分辨率, 供手动注释参考)
dfp <- as.data.table(cluster_profile, keep.rownames = "cluster")
dfp <- melt(dfp, id.vars = "cluster", variable.name = "celltype", value.name = "score")
p4 <- ggplot(dfp, aes(x = cluster, y = score, fill = celltype)) +
  geom_col(position = "fill") + theme_bw() +
  labs(title = "Neighborhood composition per CAF niche (cell type level)",
       y = "relative score")
ggsave(file.path(config$outdir, "figures", "Fig4_cluster_neighbor_comp_bar.pdf"),
       p4, width = 7, height = 4)

# Fig5 (可选): W 矩阵的 UMAP 可视化, 检查聚类结构
CAF_seu <- RunUMAP(CAF_seu, dims = 1:k_nmf, verbose = FALSE)
p5 <- DimPlot(CAF_seu, reduction = "umap", pt.size = 0.3) +
  ggtitle("CAF neighbor-factor space (Leiden clusters)")
ggsave(file.path(config$outdir, "figures", "Fig5_W_umap.pdf"), p5,
       width = 6, height = 5)

cat("[DONE] 所有结果输出至: ", config$outdir, "\n")
