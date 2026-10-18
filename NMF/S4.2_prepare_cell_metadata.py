# Export the merged celltype annotation (S3.2 non-epi + S4.1 epithelial
# tumor/normal) from all three Xenium zarrs into a single cell metadata table:
# one row per annotated cell, columns exactly as required by
# NMF/S4.3_CAF_NMF.R (和 S4.4_CAF_Leiden_clustering.R)。
# Cells still labeled "other" (not covered by either annotation, mostly
# QC-filtered cells) are excluded.
# 注意: zarr table 的 obs 里**没有** x_centroid/y_centroid 列, 空间坐标在
# table.obsm['spatial'] (与 obs 逐行对齐, 单位 um, 已与 cell_circles 几何全量核对)
import os

import numpy as np
import pandas as pd
import spatialdata as sd

# Raw data directory (containing the Xenium *.zarr stores) is set via the
# XENIUM_BASE environment variable; adjust the default below to your local
# download location.
BASE = os.environ.get("XENIUM_BASE", "/path/to/Xenium")
SAMPLES = {
    "P1": f"{BASE}/Xenium_P1_CRC_Add_on_FFPE.zarr",
    "P2": f"{BASE}/Xenium_P2_CRC_Add_on_FFPE.zarr",
    "P5": f"{BASE}/Xenium_P5_CRC_Add_on_FFPE.zarr",
}

frames = []
for sample, xenium_path in SAMPLES.items():
    sdata = sd.read_zarr(xenium_path)
    table = sdata.tables["table"]
    obs = table.obs
    if "spatial" not in table.obsm:
        table = table.to_memory()
    xy = np.asarray(table.obsm["spatial"], dtype=float)
    assert len(xy) == len(obs), "obsm['spatial'] not aligned with obs"
    df = pd.DataFrame({
        "cell_id": obs["cell_id"].astype(str),
        "sample_id": sample,
        "cell_type": obs["celltype"].astype(str),
        "x_centroid": xy[:, 0],
        "y_centroid": xy[:, 1],
    })
    n_total = len(df)
    df = df[df["cell_type"] != "other"]
    print(f"[{sample}] total cells: {n_total}, annotated: {len(df)}, "
          f"excluded 'other': {n_total - len(df)}")
    frames.append(df)

meta = pd.concat(frames, ignore_index=True)
out = f"{BASE}/NMF/input/cell_metadata.csv"
meta.to_csv(out, index=False)
print(f"wrote {len(meta)} cells to {out}")
print("cell_type counts:")
print(meta["cell_type"].value_counts())
print("per-sample counts:")
print(meta["sample_id"].value_counts())
