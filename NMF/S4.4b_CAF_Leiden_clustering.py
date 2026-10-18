#!/usr/bin/env python
# =============================================================================
# S4.4b_CAF_Leiden_clustering.py
#
# S4.3b_CAF_NMF_k5.R 的后续步骤: 对 W 矩阵 (cell x factor) 做 Leiden 聚类。
# 用 scanpy 实现 (与原文社区复现一致):
#   sc.pp.neighbors(n_neighbors=20, use_rep='X_nmf', metric='cosine')
#   sc.tl.leiden(resolution=0.1)
#
# 输入:  NMF/output_CAF/k5_r{radius}/03_NMF_W_cell_by_factor.csv
# 输出:  NMF/output_CAF/k5_r{radius}/06_leiden_clusters.csv  (cell_id, cluster)
#
# 运行 (squidpy_env 环境):
#   ~/anaconda3/envs/squidpy_env/bin/python NMF/S4.4b_CAF_Leiden_clustering.py
#
# 之后再用 R 跑 S4.4_CAF_Leiden_clustering.R 读取聚类结果做下游统计/绘图。
# =============================================================================

import argparse

import pandas as pd
import scanpy as sc


def main() -> None:
    parser = argparse.ArgumentParser(description="CAF W-matrix Leiden clustering")
    parser.add_argument("--radius", type=int, default=80, choices=[20, 40, 60, 80],
                        help="neighborhood radius (um), must match S4.3b output dir")
    parser.add_argument("--resolution", type=float, default=0.1,
                        help="Leiden resolution (larger -> more clusters)")
    args = parser.parse_args()

    rdir = f"NMF/output_CAF/k5_r{args.radius}"
    w_path = f"{rdir}/03_NMF_W_cell_by_factor.csv"
    out_path = f"{rdir}/06_leiden_clusters.csv"

    W = pd.read_csv(w_path, index_col=0)
    print(f"[INFO] W: {W.shape[0]} cells x {W.shape[1]} factors")

    adata = sc.AnnData(X=W.values)
    adata.obs_names = W.index
    adata.obsm["X_nmf"] = W.values

    sc.pp.neighbors(adata, n_neighbors=20, use_rep="X_nmf", metric="cosine")
    sc.tl.leiden(adata, resolution=args.resolution, random_state=0, key_added="leiden")

    out = pd.DataFrame({
        "cell_id": W.index,
        "cluster": adata.obs["leiden"].astype(str).values,
    })
    out.to_csv(out_path, index=False)
    print(f"[INFO] Leiden resolution = {args.resolution}, clusters:")
    print(out["cluster"].value_counts().sort_index())
    print(f"[DONE] saved to {out_path}")


if __name__ == "__main__":
    main()
