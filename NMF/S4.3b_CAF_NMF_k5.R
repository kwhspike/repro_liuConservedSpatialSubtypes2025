# =============================================================================
# S4.3b_CAF_NMF_k5.R
#
# RStudio 分步运行版本
# CAF spatial niche NMF
# 同时运行 4 个邻域半径 (20 / 40 / 60 / 80 um), 每个 radius 单独运行 k = 5 NMF
# 每个 radius 的所有分解结果 (W / H / comp / fit 对象) 都保存到独立子目录:
#   NMF/output_CAF/k5_r{radius}/
# =============================================================================


# =============================================================================
# 0. 加载包 + 配置
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dbscan)
  library(NMF)
  library(qs)
})

config <- list(
  
  # 输入
  input_path   = "NMF/input/cell_metadata.csv",
  
  sample_col   = "sample_id",
  celltype_col = "cell_type",
  x_col        = "x_centroid",
  y_col        = "y_centroid",
  
  # CAF
  target_labels = c("CAF cell"),
  
  # 邻域半径 (um): 同时运行这 4 个
  radii = c(20, 40, 60, 80),
  
  # NMF
  rank = 5,
  seed = 42,
  nmf_method = "brunet",
  nmf_nrun = 20,
  
  # 并行
  nmf_ncores = 16,
  
  # 是否复用邻域矩阵缓存 (每个 radius 有独立缓存文件)
  reuse_cache = TRUE,
  
  # 输出
  outdir = "NMF/output_CAF"
)

dir.create(
  config$outdir,
  showWarnings = FALSE,
  recursive = TRUE
)

set.seed(config$seed)

cat("[INFO] Configuration loaded.\n")
cat(sprintf("[INFO] radii = %s um, rank = %d\n",
            paste(config$radii, collapse = ", "), config$rank))

# =============================================================================
# 0b. 一次性读取 metadata + 注册 NMF 并行后端
# =============================================================================

cat("[INFO] 读取 cell metadata...\n")

meta <- fread(config$input_path)

stopifnot(
  all(
    c(
      config$sample_col,
      config$celltype_col,
      config$x_col,
      config$y_col
    ) %in% colnames(meta)
  )
)

setnames(
  meta,
  old = c(
    config$sample_col,
    config$celltype_col,
    config$x_col,
    config$y_col
  ),
  new = c(
    "sample_id",
    "cell_type",
    "x",
    "y"
  )
)

# 如果没有 cell_id，就生成一个
if (!"cell_id" %in% colnames(meta)) {
  meta[, cell_id := paste0("cell_", .I)]
}

# 删除无效数据
meta <- meta[
  !is.na(x) &
    !is.na(y) &
    !is.na(cell_type) &
    cell_type != ""
]

cat(
  sprintf(
    "[INFO] 总细胞数: %d\n",
    nrow(meta)
  )
)

cat(
  sprintf(
    "[INFO] 样本数: %d\n",
    length(unique(meta$sample_id))
  )
)

all_celltypes <- sort(
  unique(meta$cell_type)
)

cat(
  sprintf(
    "[INFO] cell type 数量: %d\n",
    length(all_celltypes)
  )
)

# NMF 并行 (注册一次, 所有 radius 共用)
nmf_opts <- list()

if (config$nmf_ncores > 1) {
  
  if (requireNamespace("doParallel", quietly = TRUE)) {
    
    doParallel::registerDoParallel(
      cores = config$nmf_ncores
    )
    
    nmf_opts <- "p"
    
    cat(
      sprintf(
        "[INFO] 使用 %d 个 CPU 核心运行 NMF\n",
        config$nmf_ncores
      )
    )
    
  } else {
    
    warning(
      "没有安装 doParallel，将使用串行 NMF"
    )
  }
}

# =============================================================================
# 主循环: 对每个 radius 构建邻域矩阵 -> NMF k = 5 -> 保存全部结果
# =============================================================================

for (radius in config$radii) {
  
  cat("\n\n")
  cat("############################################################\n")
  cat(sprintf("## radius = %d um\n", radius))
  cat("############################################################\n\n")
  
  # 该 radius 的输出子目录
  rdir <- file.path(
    config$outdir,
    sprintf("k5_r%d", radius)
  )
  
  dir.create(
    rdir,
    showWarnings = FALSE,
    recursive = TRUE
  )
  
  # ---------------------------------------------------------------------------
  # 1. 读取 / 构建 CAF 邻域组成矩阵 (radius 专属缓存)
  # ---------------------------------------------------------------------------
  
  f_cache <- file.path(
    rdir,
    sprintf("01_CAF_neighbor_comp_per_sample_r%d.qs", radius)
  )
  
  if (config$reuse_cache && file.exists(f_cache)) {
    
    cat("[INFO] 发现已有邻域矩阵缓存。\n")
    cat("[INFO] 直接读取:\n")
    cat(f_cache, "\n\n")
    
    per_sample <- qread(f_cache)
    
  } else {
    
    per_sample <- list()
    
    for (s in unique(meta$sample_id)) {
      
      cat(
        sprintf(
          "\n[INFO] 正在处理样本: %s\n",
          s
        )
      )
      
      m_s <- meta[
        sample_id == s
      ]
      
      coords <- as.matrix(
        m_s[, .(x, y)]
      )
      
      # 固定半径邻域
      nn <- dbscan::frNN(
        coords,
        eps = radius
      )
      
      # 找到 CAF
      idx_target <- which(
        m_s$cell_type %in% config$target_labels
      )
      
      if (length(idx_target) == 0) {
        
        warning(
          sprintf(
            "样本 %s 中没有 CAF cell",
            s
          )
        )
        
        next
      }
      
      # -----------------------------------------------------------------------
      # 行 = CAF cell
      # 列 = cell type
      # 值 = 邻域内该 cell type 所占比例
      # -----------------------------------------------------------------------
      
      comp <- matrix(
        0,
        nrow = length(idx_target),
        ncol = length(all_celltypes),
        dimnames = list(
          m_s$cell_id[idx_target],
          all_celltypes
        )
      )
      
      neighbor_counts <- numeric(
        length(idx_target)
      )
      
      for (i in seq_along(idx_target)) {
        
        cell_idx <- idx_target[i]
        
        nbr <- nn$id[[cell_idx]]
        
        # 删除自己
        nbr <- nbr[
          nbr != cell_idx
        ]
        
        neighbor_counts[i] <- length(nbr)
        
        if (length(nbr) == 0) {
          next
        }
        
        tab <- table(
          m_s$cell_type[nbr]
        )
        
        comp[
          i,
          names(tab)
        ] <- as.numeric(tab) / length(nbr)
      }
      
      # 删除邻域内完全没有邻居的 CAF
      keep <- rowSums(comp) > 0
      
      if (any(!keep)) {
        
        warning(
          sprintf(
            "样本 %s: %d 个 CAF 没有邻居，已删除",
            s,
            sum(!keep)
          )
        )
      }
      
      comp <- comp[
        keep,
        ,
        drop = FALSE
      ]
      
      per_sample[[as.character(s)]] <- list(
        
        raw = comp,
        
        coords = m_s[
          match(
            rownames(comp),
            cell_id
          ),
          .(
            cell_id,
            x,
            y
          )
        ]
      )
      
      cat(
        sprintf(
          "[INFO] %s: CAF = %d, 平均邻域细胞数 = %.1f\n",
          s,
          nrow(comp),
          mean(neighbor_counts)
        )
      )
    }
    
    # 保存缓存
    qsave(
      per_sample,
      f_cache
    )
    
    cat("\n[INFO] 邻域矩阵已经缓存:\n")
    cat(f_cache, "\n")
  }
  
  cat(sprintf("\n[INFO] radius = %d um 邻域矩阵构建/读取完成。\n", radius))
  
  # ---------------------------------------------------------------------------
  # 2. 合并所有样本 + 检查 NMF 输入
  # ---------------------------------------------------------------------------
  
  comp_all <- do.call(
    rbind,
    lapply(
      per_sample,
      `[[`,
      "raw"
    )
  )
  
  cat(
    sprintf(
      "[INFO] NMF 输入矩阵: %d cells x %d cell types\n",
      nrow(comp_all),
      ncol(comp_all)
    )
  )
  
  # 检查
  cat(
    sprintf(
      "[INFO] NA 数量: %d\n",
      sum(is.na(comp_all))
    )
  )
  
  cat(
    sprintf(
      "[INFO] 0 值数量: %d\n",
      sum(comp_all == 0)
    )
  )
  
  cat(
    sprintf(
      "[INFO] 非零值数量: %d\n",
      sum(comp_all > 0)
    )
  )
  
  # 每一行应该是邻域组成比例
  row_sums <- rowSums(comp_all)
  
  cat(
    sprintf(
      "[INFO] 行和范围: %.4f - %.4f\n",
      min(row_sums),
      max(row_sums)
    )
  )
  
  # 保存 NMF 输入矩阵
  f_comp <- file.path(
    rdir,
    "02_CAF_neighbor_comp.csv"
  )
  
  fwrite(
    data.table(
      cell_id = rownames(comp_all),
      as.data.table(comp_all)
    ),
    f_comp
  )
  
  cat("[INFO] NMF 输入矩阵已保存:\n")
  cat(f_comp, "\n")
  
  # ---------------------------------------------------------------------------
  # 3. 运行 NMF (k = 5)
  # ---------------------------------------------------------------------------
  
  cat("\n")
  cat("============================================================\n")
  cat(sprintf("radius = %d um, 开始 NMF\n", radius))
  cat("rank  = 5\n")
  cat("method =", config$nmf_method, "\n")
  cat("nrun  =", config$nmf_nrun, "\n")
  cat("============================================================\n\n")
  
  set.seed(config$seed)
  
  fit_k5 <- nmf(
    comp_all,
    rank = config$rank,
    method = config$nmf_method,
    nrun = config$nmf_nrun,
    seed = config$seed,
    .options = nmf_opts
  )
  
  cat(sprintf("\n[DONE] radius = %d um, NMF k = 5 完成。\n", radius))
  
  # ---------------------------------------------------------------------------
  # 4. 保存分解结果 (所有结果)
  # ---------------------------------------------------------------------------
  
  W <- basis(fit_k5)   # CAF cell x factor
  H <- coef(fit_k5)    # factor x cell type
  
  factor_names <- paste0("Factor_", seq_len(config$rank))
  
  colnames(W) <- factor_names
  rownames(H) <- factor_names
  
  # W: csv + qs
  fwrite(
    data.table(
      cell_id = rownames(W),
      as.data.table(W)
    ),
    file.path(rdir, "03_NMF_W_cell_by_factor.csv")
  )
  
  qsave(
    W,
    file.path(rdir, "03_NMF_W_cell_by_factor.qs")
  )
  
  # H: csv + qs
  fwrite(
    data.table(
      factor = rownames(H),
      as.data.table(H)
    ),
    file.path(rdir, "04_NMF_H_factor_by_celltype.csv")
  )
  
  qsave(
    H,
    file.path(rdir, "04_NMF_H_factor_by_celltype.qs")
  )
  
  # 完整 nmf 结果对象 (basis/coef/residuals 等)
  qsave(
    fit_k5,
    file.path(rdir, "05_NMF_fit_k5.qs")
  )
  
  cat(sprintf("[INFO] radius = %d um 的 W/H/fit 已保存至 %s\n", radius, rdir))
  
}

cat("\n############################################################\n")
cat("[DONE] 所有 radius (20/40/60/80 um) 的 k = 5 NMF 全部完成。\n")
cat(sprintf("结果目录: %s/k5_r{20,40,60,80}/\n", config$outdir))
cat("############################################################\n")
