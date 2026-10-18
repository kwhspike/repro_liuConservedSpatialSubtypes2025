# Add epithelial SCT_snn_res.0.4 cluster labels (exported from Seurat) as an obs
# column of the Xenium zarr table, so they can be visualized with
# spatialdata_plot / squidpy.
# Raw data directory (containing the Xenium *.zarr stores) is set via the
# XENIUM_BASE environment variable; adjust the default below to your local download location.
import os

import pandas as pd
import spatialdata as sd

BASE = os.environ.get("XENIUM_BASE", "/path/to/Xenium")
xenium_path = f"{BASE}/Xenium_P2_CRC_Add_on_FFPE.zarr"
cluster_csv = "P2_seurat/P2_epi_SCT_snn_res_0_4.csv"

OBS_COL = "SCT_snn_res_0_4"

sdata = sd.read_zarr(xenium_path)
table = sdata.tables["table"].to_memory()

clusters = pd.read_csv(cluster_csv, dtype={"cell_id": str, OBS_COL: str})
clusters = clusters.drop_duplicates(subset="cell_id").set_index("cell_id")[OBS_COL]

# align to table obs order via the cell_id column (obs index is just 0..N)
cell_ids = table.obs["cell_id"].astype(str)
mapped = cell_ids.map(clusters)

n_epi = mapped.notna().sum()
print(f"total cells in zarr table: {len(mapped)}")
print(f"epithelial cells matched:  {n_epi}")
print(f"non-epithelial cells (-> 'other'): {mapped.isna().sum()}")
assert n_epi == len(clusters), (
    f"matched {n_epi} != exported {len(clusters)}, cell_id mismatch?"
)

# non-epithelial cells have no cluster label: assign them "other"
# instead of leaving NaN, otherwise plotting renders them with the
# default color and they show up as a bogus giant "cluster".
mapped = mapped.fillna("other")

table.obs[OBS_COL] = pd.Categorical(mapped)
print("cluster counts in zarr table:")
print(table.obs[OBS_COL].value_counts(dropna=False))

# write the updated table element back into the zarr (in place)
# the on-disk element must be deleted first because the store is in use
sdata.tables["table"] = table
sdata.delete_element_from_disk("table")
sdata.write_element("table")
print("done: obs column written to", xenium_path)
