# Visualize the non-epithelial integration celltype annotations on the
# Xenium morphology images for P1, P2 and P5.
# Run with the squidpy_env conda environment.
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import spatialdata as sd
import spatialdata_plot  # noqa: F401

# Raw data directory (containing the Xenium *.zarr stores) is set via the
# XENIUM_BASE environment variable; adjust the default below to your local
# download location. Figures are written into the repository.
BASE = os.environ.get("XENIUM_BASE", "/path/to/Xenium")
OUT_DIR = "seurat_integration"
SAMPLES = {
    "P1": f"{BASE}/Xenium_P1_CRC_Add_on_FFPE.zarr",
    "P2": f"{BASE}/Xenium_P2_CRC_Add_on_FFPE.zarr",
    "P5": f"{BASE}/Xenium_P5_CRC_Add_on_FFPE.zarr",
}

OBS_COL = "celltype"
# cells without an annotation are labeled "other" in the zarr (see S3.2).
# Give "other" its own color (gray) so the legend shows the "other" label
# (converting to NaN would make the legend show a bogus "NA" entry instead).

for sample, xenium_path in SAMPLES.items():
    sdata = sd.read_zarr(xenium_path)

    cats = list(sdata.tables["table"].obs[OBS_COL].cat.categories)
    # tab20 indices 14/15 are grays (#7f7f7f / #c7c7c7); with 16 categories the
    # alphabetically 15th celltype ("Tumor cell") lands on a gray and becomes
    # indistinguishable from "other". Skip the gray entries of tab20.
    tab20_no_gray = [i for i in range(20) if i not in (14, 15)]
    palette = {
        c: mpl.colors.to_hex(mpl.colormaps["tab20"](tab20_no_gray[i]))
        for i, c in enumerate(c for c in cats if c != "other")
    }
    palette["other"] = "#808080"  # gray for unannotated cells
    palette = [palette[c] for c in cats]  # align with groups order

    # --- celltypes only (no image), clearer boundaries ---
    fig, ax = plt.subplots(figsize=(12, 12))
    (
        sdata.pl.render_shapes(
            "cell_circles", color=OBS_COL, groups=cats, palette=palette
        ).pl.show(ax=ax, title=f"{sample} non-epithelial celltype")
    )
    out_fig = f"{OUT_DIR}/{sample}_non_epi_celltype_spatial.png"
    fig.savefig(out_fig, dpi=200, bbox_inches="tight")
    print("saved:", out_fig)
    plt.close(fig)
