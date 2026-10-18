# =============================================================================
# S4.5_radius_stability.R
#
# 比较 NMF(k=5) 在 40 / 60 / 80 um 邻域半径下的输入组成矩阵稳定性
#
# 分析内容:
#   1. 读取三个半径的 cell x celltype 邻域组成矩阵 (02_CAF_neighbor_comp.csv)
#   2. 统一细胞与 cell-type 后, 逐细胞比较 composition 稳定性
#      (Pearson 相关 / MAE / Bray-Curtis 距离)
#   3. 比较各 cell type 平均比例随半径的变化
#   4. 邻居数量比较: 原始脚本未保存 neighbor_n, 本脚本从比例分数反推
#      (每个非零值 = k / n, 用有理数重构恢复分母 n), 并对 r20 与日志中的
#      均值 (~12-15) 交叉验证
#   5. 汇总为结果表 + 图
#
# 输出: NMF/output_CAF/radius_stability/
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(qs)
})

config <- list(
  radii  = c(40, 60, 80),
  outdir = "NMF/output_CAF/radius_stability"
)

dir.create(config$outdir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# 1. 读取三个半径的 composition 矩阵
# -----------------------------------------------------------------------------

read_comp <- function(radius) {
  path <- sprintf("NMF/output_CAF/k5_r%d/02_CAF_neighbor_comp.csv", radius)

  df <- fread(path, data.table = FALSE)

  cell_id <- df[[1]]
  mat <- as.matrix(df[, -1, drop = FALSE])

  storage.mode(mat) <- "numeric"
  rownames(mat) <- cell_id

  mat
}

comp <- lapply(setNames(config$radii, paste0("r", config$radii)), read_comp)

cat("[INFO] 各半径矩阵维度:\n")
print(sapply(comp, dim))

# -----------------------------------------------------------------------------
# 2. 统一细胞与 cell-type (取三者交集)
# -----------------------------------------------------------------------------

common_cells <- Reduce(intersect, lapply(comp, rownames))

# 三个矩阵列完全一致, 这里仍用 union + 补 0 以保留通用性
common_types <- Reduce(union, lapply(comp, colnames))

align_comp <- function(mat, cells, types) {
  out <- matrix(
    0,
    nrow = length(cells),
    ncol = length(types),
    dimnames = list(cells, types)
  )

  out[rownames(mat), colnames(mat)] <- mat
  out[cells, types, drop = FALSE]
}

A <- lapply(comp, align_comp, cells = common_cells, types = common_types)

cat(sprintf("[INFO] 三半径共有的内皮细胞数: %d\n", length(common_cells)))
cat(sprintf("[INFO] cell-type 列数: %d\n", length(common_types)))

# 验证比例矩阵归一化
cat("\n[INFO] 行和(归一化)检查:\n")
for (r in names(A)) {
  rs <- rowSums(A[[r]])
  cat(sprintf("  %s: min = %.4f, max = %.4f\n", r, min(rs), max(rs)))
}

# -----------------------------------------------------------------------------
# 3. 逐细胞 composition 稳定性比较
# -----------------------------------------------------------------------------

compare_comp <- function(A_, B_) {
  bc <- rowSums(abs(A_ - B_)) / 2

  data.frame(
    mean_BC   = mean(bc, na.rm = TRUE),
    median_BC = median(bc, na.rm = TRUE),
    q90_BC    = unname(quantile(bc, 0.90, na.rm = TRUE)),
    mean_MAE  = mean(rowMeans(abs(A_ - B_)), na.rm = TRUE)
  )
}

pairs <- list(
  `40_vs_60` = list(A[["r40"]], A[["r60"]]),
  `60_vs_80` = list(A[["r60"]], A[["r80"]]),
  `40_vs_80` = list(A[["r40"]], A[["r80"]])
)

pair_summary <- do.call(
  rbind,
  lapply(names(pairs), function(nm) {
    compare_comp(pairs[[nm]][[1]], pairs[[nm]][[2]])
  })
)
rownames(pair_summary) <- names(pairs)

cat("\n[INFO] 半径间逐细胞 composition 距离汇总:\n")
print(pair_summary)

fwrite(
  data.frame(comparison = rownames(pair_summary), pair_summary),
  file.path(config$outdir, "01_pairwise_distance_summary.csv")
)

# 逐细胞 40 vs 80: Pearson 相关 + MAE + BC
cell_cor_40_80 <- vapply(seq_len(nrow(A[["r40"]])), function(i) {
  x <- A[["r40"]][i, ]
  y <- A[["r80"]][i, ]

  if (sd(x) == 0 || sd(y) == 0) {
    return(NA_real_)
  }

  cor(x, y, method = "pearson")
}, numeric(1))

cell_mae_40_80 <- rowMeans(abs(A[["r40"]] - A[["r80"]]))
bc_40_80       <- rowSums(abs(A[["r40"]] - A[["r80"]])) / 2
bc_40_60       <- rowSums(abs(A[["r40"]] - A[["r60"]])) / 2
bc_60_80       <- rowSums(abs(A[["r60"]] - A[["r80"]])) / 2

cat("\n[INFO] 40 vs 80 逐细胞 Pearson 相关:\n")
print(summary(cell_cor_40_80))
cat("\n[INFO] 40 vs 80 逐细胞 MAE:\n")
print(summary(cell_mae_40_80))
cat("\n[INFO] 40 vs 80 逐细胞 Bray-Curtis:\n")
print(summary(bc_40_80))

# -----------------------------------------------------------------------------
# 4. 各 cell type 平均比例随半径的变化
# -----------------------------------------------------------------------------

mean_prop <- sapply(A, colMeans)

type_comparison <- data.frame(
  cell_type  = common_types,
  radius40   = mean_prop[, "r40"],
  radius60   = mean_prop[, "r60"],
  radius80   = mean_prop[, "r80"]
)

type_comparison$diff_40_80    <- type_comparison$radius80 - type_comparison$radius40
type_comparison$abs_diff_40_80 <- abs(type_comparison$diff_40_80)

type_comparison <- type_comparison[order(-type_comparison$abs_diff_40_80), ]

cat("\n[INFO] 各 cell type 平均比例 (按 |diff 40-80| 排序):\n")
print(type_comparison, row.names = FALSE)

fwrite(
  type_comparison,
  file.path(config$outdir, "02_celltype_mean_proportions.csv")
)

# -----------------------------------------------------------------------------
# 5. 邻居数量: 从比例分数反推 (v = k / n, 重构分母 n)
# -----------------------------------------------------------------------------

# 每行: 所有非零值共享同一分母 n (邻居数).
# 候选 n 来自每个非零值的 round(k / v) (k = 1..8), 逐一验证全行一致性.
infer_neighbor_n <- function(mat, max_n = 2000, tol = 1e-6) {
  n_cells <- nrow(mat)

  res <- integer(n_cells)

  for (i in seq_len(n_cells)) {
    v <- mat[i, ]
    v <- v[v > 0]

    if (length(v) == 0) {
      res[i] <- 0L
      next
    }

    # 候选分母: 任一值都可能是 k/n -> n ~ k/v
    cand <- unique(as.integer(round(
      as.vector(outer(1:8, v, `/`))
    )))
    cand <- cand[cand >= 1 & cand <= max_n]

    ok <- cand[
      vapply(cand, function(n) {
        max(abs(v * n - round(v * n))) < tol * n
      }, logical(1))
    ]

    if (length(ok) == 0) {
      res[i] <- NA_integer_
    } else {
      # 多个候选时取最小 (首个满足整数约束的 n)
      res[i] <- min(ok)
    }
  }

  names(res) <- rownames(mat)
  res
}

cat("\n[INFO] 从比例矩阵反推邻居数量...\n")

neighbor_n <- lapply(comp, infer_neighbor_n)

# 只保留共有细胞
neighbor_n <- lapply(neighbor_n, `[`, common_cells)

n_summary <- data.frame(
  radius  = paste0(config$radii, " um"),
  mean_n  = sapply(neighbor_n, mean,   na.rm = TRUE),
  median  = sapply(neighbor_n, median, na.rm = TRUE),
  Q1      = sapply(neighbor_n, function(x) unname(quantile(x, 0.25, na.rm = TRUE))),
  Q3      = sapply(neighbor_n, function(x) unname(quantile(x, 0.75, na.rm = TRUE))),
  min     = sapply(neighbor_n, min,    na.rm = TRUE),
  max     = sapply(neighbor_n, max,    na.rm = TRUE),
  n_NA    = sapply(neighbor_n, function(x) sum(is.na(x)))
)

cat("\n[INFO] 邻居数量分布 (反推):\n")
print(n_summary)

fwrite(
  n_summary,
  file.path(config$outdir, "03_neighbor_count_summary.csv")
)

# 交叉验证: r20 反推均值应与报告 (~12-15) 一致
comp20 <- read_comp(20)
n20 <- infer_neighbor_n(comp20)

cat("\n[INFO] [验证] r20 反推邻居数 (报告均值 ~12-15):\n")
print(summary(n20))

# -----------------------------------------------------------------------------
# 6. 保存逐细胞指标 + 图
# -----------------------------------------------------------------------------

cell_metrics <- data.frame(
  cell_id       = common_cells,
  n_40          = neighbor_n[["r40"]],
  n_60          = neighbor_n[["r60"]],
  n_80          = neighbor_n[["r80"]],
  cor_40_80     = cell_cor_40_80,
  mae_40_80     = cell_mae_40_80,
  bc_40_60      = bc_40_60,
  bc_60_80      = bc_60_80,
  bc_40_80      = bc_40_80
)

fwrite(
  cell_metrics,
  file.path(config$outdir, "04_cell_level_metrics.csv")
)

# 图: Bray-Curtis 距离分布
pdf(file.path(config$outdir, "fig1_BC_distribution.pdf"), width = 7, height = 5)

par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))

hist(bc_40_60, breaks = 60, freq = FALSE, col = "steelblue",
     main = "40 vs 60 um", xlab = "Bray-Curtis distance")
abline(v = mean(bc_40_60), col = "red", lwd = 2)
mtext(sprintf("mean = %.3f", mean(bc_40_60)), side = 3, line = 0.3, cex = 0.8)

hist(bc_60_80, breaks = 60, freq = FALSE, col = "steelblue",
     main = "60 vs 80 um", xlab = "Bray-Curtis distance")
abline(v = mean(bc_60_80), col = "red", lwd = 2)
mtext(sprintf("mean = %.3f", mean(bc_60_80)), side = 3, line = 0.3, cex = 0.8)

hist(bc_40_80, breaks = 60, freq = FALSE, col = "steelblue",
     main = "40 vs 80 um", xlab = "Bray-Curtis distance")
abline(v = mean(bc_40_80), col = "red", lwd = 2)
mtext(sprintf("mean = %.3f", mean(bc_40_80)), side = 3, line = 0.3, cex = 0.8)

dev.off()

# 图: 各 cell type 平均比例变化 (40 vs 80)
pdf(file.path(config$outdir, "fig2_celltype_prop_change.pdf"), width = 8, height = 5)

tc <- type_comparison[order(type_comparison$diff_40_80), ]

par(mar = c(5, 10, 3, 1))

barplot(
  tc$diff_40_80,
  names.arg = tc$cell_type,
  horiz = TRUE,
  las = 1,
  col = ifelse(tc$diff_40_80 > 0, "tomato", "steelblue"),
  xlab = "mean proportion difference (80 um - 40 um)",
  main = "Cell-type mean neighborhood proportion change"
)

abline(v = 0, col = "black", lwd = 1)
dev.off()

# 图: 邻居数分布 (反推)
pdf(file.path(config$outdir, "fig3_neighbor_count.pdf"), width = 7, height = 5)

par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))

for (r in names(neighbor_n)) {
  x <- neighbor_n[[r]]

  hist(x, breaks = 80, freq = FALSE, col = "grey70",
       main = paste0("radius = ", gsub("r", "", r), " um"),
       xlab = "inferred neighbor count")

  abline(v = median(x, na.rm = TRUE), col = "red", lwd = 2)
  mtext(sprintf("median = %.0f", median(x, na.rm = TRUE)),
        side = 3, line = 0.3, cex = 0.8)
}

dev.off()

# -----------------------------------------------------------------------------
# 7. 汇总打印
# -----------------------------------------------------------------------------

cat("\n\n")
cat("############################################################\n")
cat("[DONE] 半径稳定性比较完成\n")
cat(sprintf("输出目录: %s/\n", config$outdir))
cat("############################################################\n")
