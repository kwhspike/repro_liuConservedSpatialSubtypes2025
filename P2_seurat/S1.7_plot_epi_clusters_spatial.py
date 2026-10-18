# Visualize the epithelial SCT_snn_res.0.4 clusters on the Xenium morphology image.
# Run with the squidpy_env conda environment.
# Raw data directory (containing the Xenium *.zarr stores) is set via the
# XENIUM_BASE environment variable; adjust the default below to your local download location.
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import spatialdata as sd
import spatialdata_plot  # noqa: F401

BASE = os.environ.get("XENIUM_BASE", "/path/to/Xenium")
xenium_path = f"{BASE}/Xenium_P2_CRC_Add_on_FFPE.zarr"
OBS_COL = "SCT_snn_res_0_4"
# non-epithelial cells are labeled "other" in the zarr (epithelial selection, S1.4).
# Give "other" its own color (gray) so the legend shows the "other" label
# (converting to NaN would make the legend show a bogus "NA" entry instead).
OUT_FIG = "P2_seurat/P2_epi_SCT_clusters_spatial_only.png"

sdata = sd.read_zarr(xenium_path)

cats = list(sdata.tables["table"].obs[OBS_COL].cat.categories)
palette = [mpl.colors.to_hex(mpl.colormaps["tab10"](i)) for i in range(len(cats))]
palette[cats.index("other")] = "#808080"  # gray for non-epithelial cells

# --- clusters only (no image), clearer boundaries ---
fig, ax = plt.subplots(figsize=(12, 12))
(
    sdata.pl.render_shapes("cell_circles", color=OBS_COL, groups=cats, palette=palette)
    .pl.show(ax=ax, title="Epithelial SCT_snn_res.0.4 clusters")
)
fig.savefig(OUT_FIG, dpi=200, bbox_inches="tight")
print("saved:", OUT_FIG)
