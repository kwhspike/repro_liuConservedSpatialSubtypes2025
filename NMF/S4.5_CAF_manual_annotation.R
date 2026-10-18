#!/usr/bin/env Rscript
# =============================================================================
# S4.5_CAF_manual_annotation.R
#
# S4.4 聚类后的人工命名: 按 05 表各 cluster 的主导邻域成分逐个命名,
# 并归入最终 5 大分组 (与 NMF 的 5 个 factor 对应)。
#
# 输入:  NMF/output_CAF/k5_r{radius}/08_CAF_spatial_clusters_per_cell.csv
# 输出:  - 07_cluster_annotation.csv          (cluster -> annotation + group)
#        - 08_CAF_spatial_clusters_per_cell.csv  (原表增加 annotation / group 两列)
#
# 运行:  Rscript NMF/S4.5_CAF_manual_annotation.R
# 改命名: 只改下方 anno_map, 重跑即可。
# =============================================================================

suppressPackageStartupMessages(library(data.table))

radius <- 80   # 与 S4.4 config 一致
outdir <- sprintf("NMF/output_CAF/k5_r%d", radius)

## ---- cluster -> 名称 + 大分组 (依据 05 表主导邻域成分, 2026-10-08 人工命名) ----
anno_map <- data.table(
  cluster    = as.character(1:12),
  annotation = c(
    "CAF-Stromal",          # 1  Fibro 0.24 + CAF 0.24
    "CAF-Tumor-Boundary",   # 2  Tumor 0.44 + T 0.12 (肿瘤浸润边界)
    "CAF-Epi",              # 3  Epi 0.53
    "CAF-Immune",           # 4  T 0.30 + B 0.14
    "CAF-Tumor",            # 5  Tumor 0.84
    "CAF-Muscle",           # 6  Smooth muscle 0.27
    "CAF-Tumor-Core",       # 7  Tumor 0.94 (肿瘤核心)
    "CAF-Vascular",         # 8  CAF 0.37 + Fibro 0.31
    "CAF-Muscle-Neural",    # 9  Smooth muscle 0.59 + Neural 0.20
    "CAF-Epi-Adherent",     # 10 Epi 0.88 (贴 epithelium 最紧)
    "CAF-Muscle-Neural",    # 11 Smooth muscle 0.62
    "CAF-Muscle-Neural"     # 12 Smooth muscle 0.63
  ),
  group = c(
    "CAF-Stromal",          # 1
    "CAF-Tumor",            # 2
    "CAF-Epi",              # 3
    "CAF-Immune",           # 4
    "CAF-Tumor",            # 5
    "CAF-Muscle",           # 6
    "CAF-Tumor",            # 7
    "CAF-Stromal",          # 8
    "CAF-Muscle",           # 9
    "CAF-Epi",              # 10
    "CAF-Muscle",           # 11
    "CAF-Muscle"            # 12
  )
)

## ---- 读 per-cell 表, merge 注释 ----
f08 <- file.path(outdir, "08_CAF_spatial_clusters_per_cell.csv")
cell_anno <- fread(f08)

cell_anno[, cluster := as.character(cluster)]   # 08 表 cluster 为整数, 统一为字符
stopifnot(all(cell_anno$cluster %in% anno_map$cluster))
cell_anno <- merge(cell_anno, anno_map, by = "cluster", all.x = TRUE)
setorder(cell_anno, cell_id)

fwrite(anno_map, file.path(outdir, "07_cluster_annotation.csv"))
fwrite(cell_anno, f08)   # 覆盖原表, 增加 annotation / group 两列

cat("[INFO] 亚群命名 (cluster -> annotation / group):\n")
print(anno_map)
cat("\n[INFO] 各命名亚群的细胞数:\n")
print(cell_anno[, .N, by = .(annotation, group)][order(-N)])
cat("\n[DONE] 已更新: ", f08, "\n")
