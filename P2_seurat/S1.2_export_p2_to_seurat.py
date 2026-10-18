"""Export Xenium P2 zarr (spatialdata) table to 10x-style files for Seurat.

Outputs (in <out_dir>):
    matrix.mtx.gz      raw counts (cells x genes, scipy sparse)
    barcodes.tsv.gz    cell_id
    features.tsv.gz    gene_ids, gene_symbol, feature_type
    metadata.csv.gz    full obs + x_centroid/y_centroid merged from cell_circles
"""

import warnings

warnings.filterwarnings("ignore")

import gzip
import os
import shutil

import pandas as pd
import scipy.io
import scipy.sparse as sp
import spatialdata as sd

# Raw data directory (containing the Xenium *.zarr stores) is set via the
# XENIUM_BASE environment variable; adjust the default below to your local download location.
BASE = os.environ.get("XENIUM_BASE", "/path/to/Xenium")
ZARR_PATH = f"{BASE}/Xenium_P2_CRC_Add_on_FFPE.zarr"
# 10x-style export + all downstream intermediates live inside the repository
OUT_DIR = "P2_seurat/tenx"


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    sdata = sd.read_zarr(ZARR_PATH)
    table = sdata.tables["table"]
    print(f"table: {table.shape[0]} cells x {table.shape[1]} genes")

    # --- counts matrix -----------------------------------------------------
    X = table.X
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    assert sp.issparse(X), "X is not sparse"
    print("X dtype:", X.dtype, "| max:", X.max(), "| sum:", X.sum())

    # 10x matrix convention: genes x barcodes
    scipy.io.mmwrite(os.path.join(OUT_DIR, "matrix.mtx"), X.T.tocsr())
    with open(os.path.join(OUT_DIR, "matrix.mtx"), "rb") as f_in, gzip.open(
        os.path.join(OUT_DIR, "matrix.mtx.gz"), "wb"
    ) as f_out:
        shutil.copyfileobj(f_in, f_out)
    os.remove(os.path.join(OUT_DIR, "matrix.mtx"))

    # --- barcodes / features ----------------------------------------------
    cell_ids = table.obs["cell_id"].astype(str).tolist()
    with gzip.open(os.path.join(OUT_DIR, "barcodes.tsv.gz"), "wt") as f:
        f.write("\n".join(cell_ids) + "\n")

    features = pd.DataFrame(
        {
            "gene_ids": table.var["gene_ids"].astype(str),
            "gene_symbol": table.var.index.astype(str),
            "feature_type": table.var["feature_types"].astype(str),
        }
    )
    features.to_csv(
        os.path.join(OUT_DIR, "features.tsv.gz"),
        sep="\t",
        header=False,
        index=False,
        compression="gzip",
    )

    # --- metadata: obs + centroids from cell_circles -----------------------
    obs = table.obs.copy()
    cc = sdata.shapes["cell_circles"]
    centroids = pd.DataFrame(
        {"x_centroid": cc.geometry.x, "y_centroid": cc.geometry.y},
        index=cc.index.astype(str),
    )
    obs = obs.set_index("cell_id", drop=False)
    obs = obs.join(centroids, how="left")
    obs = obs.reset_index(drop=True)

    n_missing = obs["x_centroid"].isna().sum()
    print(f"cells without centroid: {n_missing}")
    assert n_missing == 0, "some cells lack coordinates!"

    obs.to_csv(os.path.join(OUT_DIR, "metadata.csv.gz"), index=False, compression="gzip")

    print("done. files in", OUT_DIR)
    for fn in sorted(os.listdir(OUT_DIR)):
        print(" ", fn, os.path.getsize(os.path.join(OUT_DIR, fn)))


if __name__ == "__main__":
    main()
