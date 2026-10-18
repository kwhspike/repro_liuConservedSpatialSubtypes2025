"""Convert official 10x Xenium output folders to SpatialData zarr stores.

For each sample:
  1. Delete redundant archives that duplicate the extracted parquet/h5 data
     (these are not used by spatialdata_io, which reads the .parquet/.h5 files).
  2. Parse the Xenium output with spatialdata_io.xenium() (boundaries, labels,
     morphology MIP + focus images, cell circles, expression table).
  3. Write to a .zarr store, read it back to verify, then delete the original
     outs folder to free disk space (user-approved).

Skipped: xenium_Colorectal (missing experiment.xenium).
"""

import os
import shutil
import sys
import time
from pathlib import Path

import spatialdata as sd
from spatialdata_io import xenium

# Directory containing the downloaded 10x Xenium *_outs folders; set via the
# XENIUM_BASE environment variable or adjust the default below.
BASE = Path(os.environ.get("XENIUM_BASE", "/path/to/Xenium"))

SAMPLES = [
    ("Xenium_V1_Human_Colon_Cancer_P1_CRC_Add_on_FFPE_outs", "Xenium_P1_CRC_Add_on_FFPE.zarr"),
    ("Xenium_V1_Human_Colon_Cancer_P2_CRC_Add_on_FFPE_outs", "Xenium_P2_CRC_Add_on_FFPE.zarr"),
    ("Xenium_V1_Human_Colon_Cancer_P5_CRC_Add_on_FFPE_outs", "Xenium_P5_CRC_Add_on_FFPE.zarr"),
]

# Files NOT used by spatialdata_io (reader uses *.parquet / cell_feature_matrix.h5).
# cell_feature_matrix/ (mtx dir), *.tar.gz, *.zarr.zip (except cells.zarr.zip),
# and *.csv.gz are all duplicates of the extracted data.
REDUNDANT_FILES = [
    "analysis.tar.gz",
    "aux_outputs.tar.gz",
    "cell_feature_matrix.tar.gz",
    "cell_feature_matrix.zarr.zip",
    "analysis.zarr.zip",
    "transcripts.zarr.zip",
    "transcripts.csv.gz",
    "cells.csv.gz",
    "cell_boundaries.csv.gz",
    "nucleus_boundaries.csv.gz",
]
REDUNDANT_DIRS = ["cell_feature_matrix"]


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def free_space() -> str:
    st = shutil.disk_usage(BASE)
    return f"{st.free / 1e9:.1f} G free"


def cleanup_redundant(src: Path) -> None:
    for name in REDUNDANT_FILES:
        f = src / name
        if f.exists():
            f.unlink()
            log(f"  removed redundant file: {name}")
    for name in REDUNDANT_DIRS:
        d = src / name
        if d.is_dir():
            shutil.rmtree(d)
            log(f"  removed redundant dir: {name}")


def convert_sample(src: Path, out: Path) -> None:
    log(f"=== {src.name} === ({free_space()})")
    cleanup_redundant(src)

    log("parsing the data ...")
    sdata = xenium(
        path=str(src),
        cells_boundaries=True,
        nucleus_boundaries=True,
        cells_as_circles=True,
        morphology_mip=True,
        morphology_focus=True,
        n_jobs=8,
    )
    log(f"parsed: {sdata}")

    if out.exists():
        shutil.rmtree(out)
    log("writing zarr ...")
    sdata.write(out)
    log(f"wrote {out} ({free_space()})")

    # verify by reading back
    log("verifying (reading zarr back) ...")
    s2 = sd.SpatialData.read(out)
    assert "table" in s2.tables, "table missing"
    n_cells = s2.tables["table"].n_obs
    n_vars = s2.tables["table"].n_vars
    elements = {
        "images": list(s2.images.keys()),
        "labels": list(s2.labels.keys()),
        "points": list(s2.points.keys()),
        "shapes": list(s2.shapes.keys()),
        "tables": list(s2.tables.keys()),
    }
    has_image = bool(s2.images)
    log(f"verified: cells={n_cells}, genes={n_vars}, images={has_image}, elements={elements}")
    assert n_cells > 0 and n_vars > 0 and has_image, "verification failed"

    # The MIP image (morphology.ome.tif, a compressed multi-page OME-TIFF) is not parsed by
    # spatialdata_io for this data version (analyzer 2.0.0). Preserve it as a standalone
    # OME-TIFF next to the zarr instead of inflating the zarr with ~30G of decompressed data.
    mip = src / "morphology.ome.tif"
    if mip.exists():
        mip_dest = BASE / (out.stem + ".morphology_mip.ome.tif")
        if mip_dest.exists():
            mip_dest.unlink()
        shutil.move(str(mip), str(mip_dest))
        log(f"preserved MIP image as {mip_dest.name}")

    log(f"deleting original folder: {src}")
    shutil.rmtree(src)
    log(f"DONE {out.name} ({free_space()})\n")


def main() -> None:
    only = sys.argv[1] if len(sys.argv) > 1 else None
    for src_name, out_name in SAMPLES:
        if only and only not in src_name:
            continue
        convert_sample(BASE / src_name, BASE / out_name)
    log("ALL DONE")


if __name__ == "__main__":
    main()
