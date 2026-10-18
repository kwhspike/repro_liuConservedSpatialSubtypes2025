# Add non-epithelial integration celltype annotations (exported from Seurat)
# as an obs column of each Xenium zarr table, so they can be visualized with
# spatialdata_plot / squidpy.
import os

import pandas as pd
import spatialdata as sd

# Raw data directory (containing the Xenium *.zarr stores) is set via the
# XENIUM_BASE environment variable; adjust the default below to your local
# download location. The annotation csvs live inside the repository
# (written by seurat_integration/S3.1_export_non_epi_celltype.R).
BASE = os.environ.get("XENIUM_BASE", "/path/to/Xenium")
SAMPLES = {
    "P1": f"{BASE}/Xenium_P1_CRC_Add_on_FFPE.zarr",
    "P2": f"{BASE}/Xenium_P2_CRC_Add_on_FFPE.zarr",
    "P5": f"{BASE}/Xenium_P5_CRC_Add_on_FFPE.zarr",
}

OBS_COL = "celltype"

for sample, xenium_path in SAMPLES.items():
    cluster_csv = f"seurat_integration/{sample}_non_epi_celltype.csv"

    sdata = sd.read_zarr(xenium_path)
    table = sdata.tables["table"].to_memory()

    clusters = pd.read_csv(cluster_csv, dtype={"cell_id": str, OBS_COL: str})
    clusters = clusters.drop_duplicates(subset="cell_id").set_index("cell_id")[OBS_COL]

    # align to table obs order via the cell_id column (obs index is just 0..N)
    cell_ids = table.obs["cell_id"].astype(str)
    mapped = cell_ids.map(clusters)

    n_matched = mapped.notna().sum()
    print(f"[{sample}] total cells in zarr table: {len(mapped)}")
    print(f"[{sample}] annotated cells matched:    {n_matched}")
    print(f"[{sample}] unmatched cells (-> 'other'): {mapped.isna().sum()}")
    assert n_matched == len(clusters), (
        f"[{sample}] matched {n_matched} != exported {len(clusters)}, cell_id mismatch?"
    )

    # cells without an annotation (e.g. epithelial-only subsets) get "other"
    # instead of NaN, otherwise plotting renders them with the default color
    # and they show up as a bogus giant category.
    mapped = mapped.fillna("other")

    table.obs[OBS_COL] = pd.Categorical(mapped)
    print(f"[{sample}] celltype counts in zarr table:")
    print(table.obs[OBS_COL].value_counts(dropna=False))

    # write the updated table element back into the zarr (in place)
    # the on-disk element must be deleted first because the store is in use
    sdata.tables["table"] = table
    sdata.delete_element_from_disk("table")
    sdata.write_element("table")
    print(f"[{sample}] done: obs column written to", xenium_path)
