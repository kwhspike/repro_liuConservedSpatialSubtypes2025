#!/usr/bin/env Rscript
# =============================================================================
# S4.3_CAF_NMF.R  (由原 S4.3_rank_survey.R + S4.4_CAF_spatial_subtype_NMF.R 合并)
#
# 内皮细胞 (CAF) 空间 niche 分析: 邻域向量 + NMF, **全部细胞** (不抽样)
#   1) 以每个内皮细胞为圆心, 80 um 半径圆形邻域, 统计邻域细胞类型组成
#      -> cell-by-celltype 原始组成矩阵 (按原文, NMF 前不做归一化)
#   2) nmfEstimateRank 扫 k = rank_min:rank_max (brunet + nrun 次重复),
#      保存 cophenetic / dispersion / rss 等指标 (Fig0 + 00_rank_summary.csv)
#   3) 把 **每个 k 的最优分解结果全部导出**: W (cell-by-factor) 和
#      H (factor-by-celltype) 保存到 output_CAF/NMF_k{k}/ 子目录,
#      供对比后选择合适的 k 做下游 Leiden 聚类 (S4.4)
#
# 注意: 样本级归一化 (邻域组成 / 全切片组成) 按原文在 NMF 之后的 S4.5 中进行。
#
# 输入: NMF/input/cell_metadata.csv (由 NMF/S4.2_prepare_cell_metadata.py 生成)
#
# 运行 (建议限制 BLAS 线程, 避免并行 worker 互相抢占):
#   OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
#     nohup Rscript NMF/S4.3_CAF_NMF.R > NMF/output_CAF/S4.3_NMF.log 2>&1 &
#
# 完成后: 看 figures/Fig0_rank_survey.pdf + 00_rank_summary.csv 选 k,
#         本流程选定 k = 5: 运行 S4.3b_CAF_NMF_k5.R (固定 k=5, 4 个半径),
#         再按 S4.4 头部说明跑 S4.4b (scanpy Leiden) -> S4.4 (统计/绘图)
#         -> S4.5_CAF_manual_annotation.R (人工命名)。
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dbscan)      # frNN: 固定半径邻域搜索
  library(NMF)         # nmfEstimateRank(), 需 v0.27
  library(qs)          # qsave/qread
})

## ============================ 0. 配置 ======================================
config <- list(
  ## ---- 输入 ----
  input_path   = "NMF/input/cell_metadata.csv",
  sample_col   = "sample_id",
  celltype_col = "cell_type",
  x_col        = "x_centroid",    # 空间坐标, 单位 um
  y_col        = "y_centroid",

  ## ---- 分析对象: 哪些 cell_type 标签算作内皮细胞 ----
  target_labels = c("CAF cell"),

  ## ---- 邻域 / NMF 参数 ----
  radius     = 80,        # 邻域半径 (um)
  seed       = 42,
  rank_min   = 2,
  rank_max   = 10,
  nmf_method = "brunet",
  survey_nrun  = 20,      # rank 扫描每个 k 的重复运行次数
  nmf_ncores   = 16,      # NMF 并行核数 (1 = 串行; 需 R 包 doParallel)

  ## ---- 缓存: 邻域矩阵构建约需数分钟, 重跑本脚本时默认直接复用已保存的 01 结果 ----
  reuse_cache = TRUE,

  outdir = "NMF/output_CAF"
)

dir.create(config$outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(config$outdir, "figures"), showWarnings = FALSE)
set.seed(config$seed)

## ======================= 1. 构建邻域组成矩阵 ==============================
f_cache <- file.path(config$outdir, "01_CAF_neighbor_comp_per_sample.qs")
if (config$reuse_cache && file.exists(f_cache)) {
  per_sample <- qread(f_cache)
  cat("[INFO] 复用缓存的邻域矩阵: ", f_cache, "\n")
} else {
  meta <- fread(config$input_path)
  stopifnot(all(c(config$sample_col, config$celltype_col,
                  config$x_col, config$y_col) %in% colnames(meta)))
  setnames(meta,
           old = c(config$sample_col, config$celltype_col, config$x_col, config$y_col),
           new = c("sample_id", "cell_type", "x", "y"))
  meta <- meta[!is.na(x) & !is.na(y) & !is.na(cell_type) & cell_type != ""]
  cat(sprintf("[INFO] 总细胞数: %d, 样本数: %d\n",
              nrow(meta), length(unique(meta$sample_id))))

  all_celltypes <- sort(unique(meta$cell_type))
  cat(sprintf("[INFO] 细胞类型 (%d): %s\n",
              length(all_celltypes), paste(all_celltypes, collapse = ", ")))

  per_sample <- list()
  for (s in unique(meta$sample_id)) {
    m_s <- meta[sample_id == s]
    coords <- as.matrix(m_s[, .(x, y)])
    # 固定半径邻域 (frNN 返回含自身的邻居, 需剔除自身)
    nn <- dbscan::frNN(coords, eps = config$radius)
    idx_target <- which(m_s$cell_type %in% config$target_labels)

    if (length(idx_target) == 0) {
      warning(sprintf("[WARN] 样本 %s 中没有目标细胞 (请检查 target_labels)", s))
      next
    }

    # 邻域组成矩阵: 行 = 内皮细胞, 列 = 细胞类型, 值为邻域内比例
    comp <- matrix(0, nrow = length(idx_target), ncol = length(all_celltypes),
                   dimnames = list(m_s$cell_id[idx_target], all_celltypes))
    for (i in seq_along(idx_target)) {
      nbr <- nn$id[[idx_target[i]]]
      nbr <- nbr[nbr != idx_target[i]]            # 剔除自身
      if (length(nbr) == 0) next
      tab <- table(m_s$cell_type[nbr])
      comp[i, names(tab)] <- as.numeric(tab) / length(nbr)
    }

    # 剔除邻域内没有任何邻居的内皮细胞
    keep <- rowSums(comp) > 0
    if (any(!keep))
      warning(sprintf("[WARN] 样本 %s: %d 个内皮细胞在 %g um 内无邻居, 已剔除",
                      s, sum(!keep), config$radius))
    comp <- comp[keep, , drop = FALSE]

    per_sample[[s]] <- list(
      raw    = comp,   # 原始邻域组成 (NMF 输入; 样本级归一化在 S4.5 聚类后进行)
      coords = m_s[match(rownames(comp), cell_id), .(cell_id, x, y)]
    )
    cat(sprintf("[INFO] 样本 %s: 内皮细胞 = %d, 平均邻域细胞数 = %.1f\n",
                s, nrow(comp),
                mean(sapply(nn$id[idx_target], length) - 1)))
  }

  qsave(per_sample, f_cache)
}

comp_all <- do.call(rbind, lapply(per_sample, `[[`, "raw"))
cat(sprintf("[INFO] 内皮-邻域矩阵 (全部细胞): %d 行 x %d 列\n",
            nrow(comp_all), ncol(comp_all)))
fwrite(data.table(cell_id = rownames(comp_all), as.data.table(comp_all)),
       file.path(config$outdir, "02_CAF_neighbor_comp.csv"))

## ============================ 2. rank 扫描 =================================
# nrun 次重复分解可并行 (foreach 后端); 未装 doParallel 时自动回退串行
nmf_opts <- list()
if (config$nmf_ncores > 1) {
  if (requireNamespace("doParallel", quietly = TRUE)) {
    doParallel::registerDoParallel(cores = config$nmf_ncores)
    nmf_opts <- "p"   # 勿加 "t" 进度输出: NMF v0.27 在 R>=4.2 下 verbose 有 bug
    cat(sprintf("[INFO] NMF 并行: doParallel 后端, %d 核\n", config$nmf_ncores))
  } else {
    warning("[WARN] 未安装 doParallel, NMF 串行运行。install.packages('doParallel')")
  }
}
cat(sprintf("[INFO] 全细胞 rank 扫描: k = %d:%d, method = %s, nrun = %d (耗时较长, 请耐心等待) ...\n",
            config$rank_min, config$rank_max, config$nmf_method, config$survey_nrun))
est <- nmfEstimateRank(comp_all,
                       range = seq(config$rank_min, config$rank_max),
                       method = config$nmf_method,
                       nrun = config$survey_nrun,
                       seed = config$seed,
                       .options = nmf_opts)

## 指标汇总表 (cophenetic / dispersion / rss / evar ...)
rank_stats <- summary(est)
print(rank_stats)
fwrite(as.data.table(rank_stats, keep.rownames = "rank"),
       file.path(config$outdir, "00_rank_summary.csv"))

pdf(file.path(config$outdir, "figures", "Fig0_rank_survey.pdf"),
    width = 10, height = 6)
plot(est)   # cophenetic / dispersion / rss 随 k 变化; 选稳定性恶化前的最大 k
dev.off()

qsave(est, file.path(config$outdir, "00_rank_survey.qs"))

## ============ 3. 导出每个 k 的最优分解 (W/H) 供选 k 对比 / 下游分析 ============
for (k in seq(config$rank_min, config$rank_max)) {
  fit_k <- fit(est[[as.character(k)]])   # 该 rank 下目标函数最优的一次分解
  W <- basis(fit_k)                      # cell-by-factor
  H <- coef(fit_k)                       # factor-by-celltype
  colnames(W) <- paste0("Factor_", seq_len(k))
  rownames(H) <- paste0("Factor_", seq_len(k))

  kdir <- file.path(config$outdir, sprintf("NMF_k%d", k))
  dir.create(kdir, showWarnings = FALSE, recursive = TRUE)
  fwrite(data.table(cell_id = rownames(W), as.data.table(W)),
         file.path(kdir, "03_NMF_W_cell_by_factor.csv"))
  qsave(W, file.path(kdir, "03_NMF_W_cell_by_factor.qs"))
  fwrite(data.table(factor = rownames(H), as.data.table(H)),
         file.path(kdir, "04_NMF_H_factor_by_celltype.csv"))
  cat(sprintf("[INFO] k = %d 的 W/H 已导出至 %s\n", k, kdir))
}

cat("[DONE] 全部完成。\n")
cat("  1) 用 figures/Fig0_rank_survey.pdf + 00_rank_summary.csv 选 k\n")
cat("  2) 本流程选定 k = 5: Rscript NMF/S4.3b_CAF_NMF_k5.R\n")
cat("  3) 再跑 S4.4b (python Leiden) 和 S4.4 (R 统计/绘图):\n")
cat("     ~/anaconda3/envs/squidpy_env/bin/python NMF/S4.4b_CAF_Leiden_clustering.py\n")
cat("     Rscript NMF/S4.4_CAF_Leiden_clustering.R\n")
cat("  4) 人工命名: Rscript NMF/S4.5_CAF_manual_annotation.R\n")
