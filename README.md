# Xenium CRC CAF 空间邻域与 NMF 分析流程

本仓库参考 Liu et al. [1] 提出的 **CAF 空间邻域分析框架**，将邻域细胞类型组成、
NMF 与 Leiden 聚类用于 10x Genomics **Xenium Human Colon Cancer (CRC) Add-on FFPE**
公开数据的 3 个样本 (P1 / P2 / P5)。这是针对该数据集的流程实现，不是对原文多平台、
多癌种研究的完整复现。

---

## 1. 数据

- **下载**:10x Genomics 数据集
  [Human CRC (dataset-human-crc)](https://www.10xgenomics.com/platforms/visium/product-family/dataset-human-crc),
  需要其中 3 个 Xenium Add-on FFPE 样本:
  `Xenium_V1_Human_Colon_Cancer_P1_CRC_Add_on_FFPE_outs`、
  `Xenium_V1_Human_Colon_Cancer_P2_CRC_Add_on_FFPE_outs`、
  `Xenium_V1_Human_Colon_Cancer_P5_CRC_Add_on_FFPE_outs`。
- 数据**不包含**在本仓库中。下载解压后得到一个包含 3 个 `*_outs/` 文件夹的目录。
- 所有脚本通过环境变量 **`XENIUM_BASE`** 定位该目录:

  ```bash
  export XENIUM_BASE=/path/to/your/xenium/download   # 指向 *_outs 文件夹的父目录
  ```

  也可以在每个脚本内直接修改 `BASE` 的默认值。
- 所有脚本**统一在仓库根目录下运行**;中间产物 (`.rds` / `.qs` / `tenx/` /
  `NMF/output_CAF/` 等) 全部生成在仓库内部,已被 `.gitignore` 排除。

## 2. 环境

| 组件 | 版本 | 安装 |
|---|---|---|
| R | 4.4.1 (Ubuntu 20.04 验证,见 `sessionInfo_output.txt`) | `Rscript install.R` |
| Python | 3.11 (conda 环境 `squidpy_env`) | `conda env create -f environment.yml` |

- R 侧关键包:Seurat **v5**、data.table、qs、dbscan、**NMF 0.27**、pheatmap、
  ggplot2、future、doParallel、plyr、dplyr、reticulate(`install.R` 自动安装)。
- R 脚本 `NMF/S4.4_CAF_Leiden_clustering.R` 的 Leiden 聚类
  (Seurat `algorithm = 4`) 通过 reticulate 调用 python 包 `leidenalg`,
  需要先让 reticulate 指向 conda 环境:

  ```bash
  conda activate squidpy_env
  export RETICULATE_PYTHON=$CONDA_PREFIX/bin/python
  ```

## 3. 分析流程

### S1 — 单样本处理 (每个样本结构相同,以 P1 为例)

| 步骤 | 脚本 | 说明 |
|---|---|---|
| S1.1 | `S1.1_convert_xenium_to_zarr.py` | 10x `*_outs/` → SpatialData `.zarr` (三个样本) |
| S1.2 | `P1_seurat/S1.2_export_p1_to_seurat.py` | zarr table → 10x 格式 (`P1_seurat/tenx/`) |
| S1.2 | `P1_seurat/S1.2_build_seurat.R` | 读入 10x 文件,构建 Seurat v5 对象 + QC |
| S1.3 | `P1_seurat/S1.3_P1_process.R` | SCTransform → PCA → 聚类 → UMAP |
| S1.4 | `P1_seurat/S1.4_P1_select_epi.R` | marker DotPlot,人工挑选上皮簇,拆分上皮/非上皮 |
| S1.5 | `P1_seurat/S1.5_export_epi_clusters.R` | 导出上皮簇标签 csv |
| S1.6 | `P1_seurat/S1.6_add_epi_clusters_to_zarr.py` | 簇标签写回 zarr (非上皮细胞标 `other`) |
| S1.7 | `P1_seurat/S1.7_plot_epi_clusters_spatial.py` | 空间图可视化 |

上皮聚类 resolution:**P1 = 0.3,P2 = 0.3,P5 = 0.3**

### S2 — 跨样本整合 (Seurat rpca integration)

| 步骤 | 脚本 | 说明 |
|---|---|---|
| S2.1 | `seurat_integration/S2.1_inter_after_selected_epi.R` | 三样本**上皮**整合 (`integra_epi.qs`) |
| S2.2 | `seurat_integration/S2.2_intigre_after_selected.R` | 三样本**非上皮**整合 → 聚类 → 大类细胞注释 (`integra_non_epi_process.qs`) |

### S3 — 非上皮细胞类型注释写回 zarr

| 步骤 | 脚本 | 说明 |
|---|---|---|
| S3.1 | `seurat_integration/S3.1_export_non_epi_celltype.R` | 导出非上皮 celltype csv |
| S3.2 | `seurat_integration/S3.2_add_non_epi_celltype_to_zarr.py` | 写入 zarr 的 `celltype` 列 |
| S3.3 | `seurat_integration/S3.3_plot_non_epi_celltype_spatial.py` | 空间图可视化 |

### S4 — CAF 空间 niche:邻域组成 + NMF

| 步骤 | 脚本 | 说明 |
|---|---|---|
| S4.1 | `seurat_integration/S4.1_merge_epi_celltype_to_zarr.py` | 依据 S1.7 空间图把上皮簇分为 Tumor / Normal,合并进 zarr 的 `celltype` 列 (Normal 簇:P1 {4,11},P2 {6,14},P5 {2,8,17}) |
| S4.2 | `NMF/S4.2_prepare_cell_metadata.py` | 从三个 zarr 汇总所有已注释细胞的 metadata (`NMF/input/cell_metadata.csv`) |
| S4.3 | `NMF/S4.3_CAF_NMF.R` | rank 扫描 (k = 2–10, brunet, nrun = 20),输出 cophenetic/dispersion/rss 指标 |
| S4.3b | `NMF/S4.3b_CAF_NMF_k5.R` | **本流程 canonical**:固定 k = 5,同时跑 4 个邻域半径 (20/40/60/80 µm) |
| S4.4b | `NMF/S4.4b_CAF_Leiden_clustering.py` | scanpy 对 W 矩阵做 Leiden 聚类 (复现原文社区实现) |
| S4.4 | `NMF/S4.4_CAF_Leiden_clustering.R` | R 侧 Leiden (Seurat algorithm 4) + 统计/绘图 (Fig1–Fig5) |
| S4.5 | `NMF/S4.5_CAF_manual_annotation.R` | 按各 cluster 主导邻域成分人工命名 (改 `anno_map` 后重跑即可) |
| S4.5 | `NMF/S4.5_radius_stability.R` | 40/60/80 µm 半径的邻域组成稳定性比较 |

关键参数:NMF **k = 5**,canonical 邻域半径 **80 µm**,Leiden resolution
0.05 (S4.4)/ 0.1 (S4.4b,scanpy 复现对照)。NMF 运行建议限制 BLAS 线程:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  nohup Rscript NMF/S4.3_CAF_NMF.R > NMF/output_CAF/S4.3_NMF.log 2>&1 &
```

### 运行顺序一览

```text
S1.1 → (P1|P2|P5) × [S1.2 py → S1.2 R → S1.3 → S1.4 → S1.5 → S1.6 → S1.7]
     → S2.1 → S2.2 → S3.1 → S3.2 → S3.3
     → S4.1 → S4.2 → S4.3(选 k) → S4.3b → S4.4b → S4.4 → S4.5(命名)
     → (可选) S4.5_radius_stability
```

Python 脚本用 `squidpy_env` 运行,例如:

```bash
$CONDA_PREFIX/bin/python P1_seurat/S1.2_export_p1_to_seurat.py
```

## 4. 目录结构

```text
├── S1.1_convert_xenium_to_zarr.py     # 原始 outs → zarr
├── environment.yml                    # conda squidpy_env (python 3.11)
├── install.R                          # R 依赖安装 (R 4.4.1)
├── sessionInfo_output.txt             # 运行环境记录
├── R/
│   └── crc_marker_list.R              # 共享的 CRC/TME marker 基因列表
├── P1_seurat/  P2_seurat/  P5_seurat/ # 单样本流程 (脚本同名, P* 区分)
├── seurat_integration/                # S2 整合 + S3 注释写回 + S4.1 上皮注释
└── NMF/                               # S4 CAF niche NMF 全部分析
```

## 5. 引用

[1] Liu Y, Sinjab A, Min J, Han G, Paradiso F, Zhang Y, et al. Conserved spatial
subtypes and cellular neighborhoods of cancer-associated fibroblasts revealed by
single-cell spatial multi-omics. *Cancer Cell*. 2025;43(5):905–924.e6.
doi:10.1016/j.ccell.2025.03.004.

- 数据:10x Genomics Xenium Human Colon Cancer Add-on FFPE,
  <https://www.10xgenomics.com/platforms/visium/product-family/dataset-human-crc>

## License

[MIT](LICENSE)
