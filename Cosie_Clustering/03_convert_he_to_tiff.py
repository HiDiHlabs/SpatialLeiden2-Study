"""
Convert the full-resolution H&E microscope images (JPEG) into pyramidal TIFFs that
lazyslide/wsidata can open.

Input: data/01_mouse_brain_MALDI_visium/he_fullres/*.jpg, downloaded with
downloads/01_mouse_brain_HE_fullres_figshare.py. These are the images Space Ranger was run
on, so pxl_row/pxl_col in tissue_positions_list.csv are pixel coordinates in them.

The JPEGs carry no physical resolution, so it is derived from the Visium grid: neighbouring
spots in the same array row (array_col +- 2) are 100 um apart. The result is written into the
TIFF resolution tags, from which openslide/tiffslide (and thus lazyslide) read the mpp.

Run with the lazyslide_gpu2 env (has pyvips) from the SpatialLeiden2-Study folder:
    /dh-projects/ag-ishaque/analysis/malta/envs/envs/lazyslide_gpu2/bin/python Cosie_Clustering/03_convert_he_to_tiff.py
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pyvips

data_dir = Path("data/01_mouse_brain_MALDI_visium")
jpg_dir = data_dir / "he_fullres"
out_dir = data_dir / "he_fullres_tiff"
out_dir.mkdir(parents=True, exist_ok=True)

slide_dir = data_dir / "sma" / "V11L12-109"
SECTIONS = {
    "V11L12-109_A1": slide_dir / "V11L12-109_A1/output_data/V11L12-109_A1_RNA/outs/spatial",
    "V11L12-109_B1": slide_dir / "V11L12-109_B1/output_data/V11L12-109_B1_RNA/outs/spatial",
    "V11L12-109_C1": slide_dir / "V11L12-109_C1/output_data/V11L12-109_C1_RNA/outs/spatial",
    "V11L12-109_D1": slide_dir / "V11L12-109_D1/output_data/V11L12-109_D1/outs/spatial",
}
spot_spacing_um = 100  # Visium centre-to-centre distance


def microns_per_pixel(spatial_dir):
    pos = pd.read_csv(
        spatial_dir / "tissue_positions_list.csv", header=None,
        names=["barcode", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col"],
    )
    # pairs of spots that are direct neighbours within an array row
    right = pos.assign(array_col=pos["array_col"] - 2)
    pairs = pos.merge(right, on=["array_row", "array_col"], suffixes=("", "_nb"))
    dist = np.hypot(pairs["pxl_row"] - pairs["pxl_row_nb"], pairs["pxl_col"] - pairs["pxl_col_nb"])
    return spot_spacing_um / np.median(dist), pos


for section, spatial_dir in SECTIONS.items():
    array_id, capture_area = section.split("_")
    jpgs = list(jpg_dir.glob(f"*{array_id}.{capture_area}-Spot*.jpg"))
    if len(jpgs) != 1:
        print(f"{section}: expected one JPEG in {jpg_dir}, found {len(jpgs)}; skipped")
        continue

    mpp, pos = microns_per_pixel(spatial_dir)
    img = pyvips.Image.new_from_file(str(jpgs[0]), access="sequential")

    # spot coordinates must fall inside the image, otherwise this is not the Space Ranger input image
    assert pos["pxl_col"].max() < img.width and pos["pxl_row"].max() < img.height, (
        f"{section}: spot coordinates exceed the image size {img.width} x {img.height}"
    )
    scalef = json.loads((spatial_dir / "scalefactors_json.json").read_text())["tissue_hires_scalef"]
    print(
        f"{section}: {img.width} x {img.height} px, {mpp:.4f} um/px "
        f"(hires image {mpp / scalef:.2f} um/px), {jpgs[0].name}"
    )

    out = out_dir / f"{section}.tiff"
    px_per_mm = 1000 / mpp
    img.tiffsave(
        str(out), tile=True, tile_width=256, tile_height=256, pyramid=True, bigtiff=True,
        compression="jpeg", Q=90, xres=px_per_mm, yres=px_per_mm, resunit="cm",
    )
    (out_dir / f"{section}.json").write_text(json.dumps({"mpp": mpp, "source": jpgs[0].name}, indent=1))
    print(f"  -> {out}")
