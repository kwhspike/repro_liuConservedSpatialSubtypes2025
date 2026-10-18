# Merge epithelial tumor/normal annotation into the existing "celltype" obs
# column of each Xenium zarr table.
#
# Normal epithelial clusters (from spatial visualization of the epithelial
# SCT clustering, see P*_seurat/S1.7):
#   P1 (SCT_snn_res_0.3): 4, 11
#   P2 (SCT_snn_res_0.4): 14, 6
#   P5 (SCT_snn_res_0.3): 2, 8, 17
# Normal epithelial -> "Epithelial cell"; all other epithelial clusters ->
# "Tumor cell" (same labels as the non-epi annotation, S3).
# Cells not covered by the epithelial csv keep their current celltype.
import os

import pandas as pd
import spatialdata as sd

# Raw data directory (containing the Xenium *.zarr stores) is set via the
# XENIUM_BASE environment variable; adjust the default below to your local
# download location. The epithelial cluster csvs live inside the repository
# (written by P*_seurat/S1.5_export_epi_clusters.R).
BASE = os.environ.get("XENIUM_BASE", "/path/to/Xenium")
SAMPLES = {
    "P1": {
        "zarr": f"{BASE}/Xenium_P1_CRC_Add_on_FFPE.zarr",
        "csv": "P1_seurat/P1_epi_SCT_snn_res_0.3.csv",
        "col": "SCT_snn_res_0.3",
        "normal": {"4", "11"},
    },
    "P2": {
        "zarr": f"{BASE}/Xenium_P2_CRC_Add_on_FFPE.zarr",
        "csv": "P2_seurat/P2_epi_SCT_snn_res_0_4.csv",
        "col": "SCT_snn_res_0_4",
        "normal": {"6", "14"},
    },
    "P5": {
        "zarr": f"{BASE}/Xenium_P5_CRC_Add_on_FFPE.zarr",
        "csv": "P5_seurat/P5_epi_SCT_snn_res_0.3.csv",
        "col": "SCT_snn_res_0.3",
        "normal": {"2", "8", "17"},
    },
}

OBS_COL = "celltype"

for sample, cfg in SAMPLES.items():
    sdata = sd.read_zarr(cfg["zarr"])
    table = sdata.tables["table"].to_memory()

    epi = pd.read_csv(cfg["csv"], dtype={"cell_id": str, cfg["col"]: str})
    epi = epi.drop_duplicates(subset="cell_id").set_index("cell_id")[cfg["col"]]
    # normal epithelial -> "Epithelial cell", other epi clusters -> "Tumor cell"
    epi = epi.map(lambda c: "Epithelial cell" if c in cfg["normal"] else "Tumor cell")

    cell_ids = table.obs["cell_id"].astype(str)
    current = table.obs[OBS_COL].astype(str)
    is_epi = cell_ids.isin(epi.index)

    # sanity: epi cells must currently be unannotated ("other"), so we never
    # overwrite a non-epithelial annotation
    prev = current[is_epi]
    n_other = (prev == "other").sum()
    print(f"[{sample}] epi cells: {is_epi.sum()}, currently 'other': {n_other}")
    assert n_other == is_epi.sum(), (
        f"[{sample}] {is_epi.sum() - n_other} epi cells already have a "
        f"non-'other' celltype, refusing to overwrite"
    )
    assert is_epi.sum() == len(epi), (
        f"[{sample}] matched {is_epi.sum()} != exported {len(epi)}, cell_id mismatch?"
    )

    mapped = current.copy()
    mapped[is_epi] = cell_ids[is_epi].map(epi)

    table.obs[OBS_COL] = pd.Categorical(mapped)
    print(f"[{sample}] celltype counts after merge:")
    print(table.obs[OBS_COL].value_counts(dropna=False))

    # write the updated table element back into the zarr (in place)
    # the on-disk element must be deleted first because the store is in use
    sdata.tables["table"] = table
    sdata.delete_element_from_disk("table")
    sdata.write_element("table")
    print(f"[{sample}] done: obs column written to", cfg["zarr"])
