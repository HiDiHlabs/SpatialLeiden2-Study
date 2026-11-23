"""data_inventory.csv: one row per downloaded Kuppe 2022 file.

Columns: path, set, zenodo_record, modality, sample, zone, patient, n_obs, n_vars,
size_gb, loaded, note. "loaded" says whether the file opened cleanly; large R objects
outside the Visium scope are listed but not opened.

Usage: python 04_inventory.py <data/kuppe_2022 dir> <output dir>
"""

import re
import sys
from pathlib import Path

import anndata as ad
import pandas as pd
from PIL import Image

Image.MAX_IMAGE_PIXELS = None
DATA, OUT = Path(sys.argv[1]), Path(sys.argv[2])
RECORD = {
    "processed1": "6578047",
    "processed2": "7098004",
    "visium_raw": "6580069",
    "he_images": "6778620",
    "metadata": "6578617",
}

vis = pd.read_csv(DATA / "visium_raw" / "metadata-Visium.csv")
by_hca = vis.set_index("hca_sample_id")
by_region = vis.drop_duplicates("patient_region_id").set_index("patient_region_id")
he = pd.read_excel(DATA / "he_images" / "metadata.xlsx")
he["file"] = he.sample_id.map(
    lambda s: (
        f"{s.split('_')[1]}_2019{s.split('_')[0]}.tif"
        if s.startswith("AKK")
        else f"{s}.tif"
    )
)
he = he.set_index("file")


def region_info(region):
    if region in by_region.index:
        r = by_region.loc[region]
        return region, r.major_labl, r.patient
    return region, None, None


rows = []
for set_dir in sorted(p for p in DATA.iterdir() if p.is_dir() and p.name in RECORD):
    for f in sorted(set_dir.iterdir()):
        if f.name.endswith(".part"):
            rows.append(
                {
                    "path": str(f.relative_to(DATA)),
                    "set": set_dir.name,
                    "loaded": "incomplete download",
                }
            )
            continue
        row = {
            "path": str(f.relative_to(DATA)),
            "set": set_dir.name,
            "zenodo_record": RECORD[set_dir.name],
            "size_gb": round(f.stat().st_size / 1e9, 3),
        }
        try:
            if f.suffix == ".h5ad":
                a = ad.read_h5ad(f, backed="r")
                row.update(n_obs=a.n_obs, n_vars=a.n_vars, loaded="yes")
                a.file.close()
                if f.name.startswith("Visium_"):
                    reg = f.stem.removeprefix("Visium_")
                    match = [r for r in by_region.index if r.replace("/", "_") == reg]
                    row["sample"], row["zone"], row["patient"] = region_info(
                        match[0] if match else reg
                    )
                    row["modality"] = "Visium (processed h5ad)"
                else:
                    row["modality"] = "snRNA-seq" if "snRNA" in f.name else "snATAC-seq"
                    row["sample"] = "all"
            elif f.name.endswith(".tar.gz"):
                hca = re.sub(r"10X_Visium_(.*)\.tar\.gz", r"\1", f.name)
                spots = DATA / "visium_tables" / f"{hca}_spots.csv"
                row["modality"] = "Visium (Space Ranger + Seurat)"
                row["sample"], row["zone"], row["patient"] = region_info(
                    by_hca.loc[hca].patient_region_id
                )
                row["note"] = f"hca_sample_id {hca}"
                if spots.exists():
                    row.update(
                        n_obs=len(pd.read_csv(spots, usecols=["barcode"])),
                        loaded="yes (Seurat object read)",
                    )
                else:
                    row["loaded"] = "not yet extracted"
            elif f.suffix == ".tif":
                with Image.open(f) as im:
                    row.update(loaded="yes", note=f"{im.size[0]} x {im.size[1]} px")
                row["modality"] = "H&E image"
                if f.name in he.index:
                    row["sample"], row["zone"], row["patient"] = region_info(
                        he.loc[f.name].patient_region_id
                    )
            elif f.suffix == ".Rds":
                row.update(
                    modality="snRNA-seq + snATAC-seq (per cell type, R)",
                    sample="all",
                    loaded="not opened (outside Visium scope)",
                )
            elif f.suffix in (".xlsx", ".csv"):
                x = pd.read_excel(f) if f.suffix == ".xlsx" else pd.read_csv(f)
                row.update(
                    modality="metadata", n_obs=len(x), n_vars=x.shape[1], loaded="yes"
                )
            else:
                row.update(modality="other", loaded="not opened")
        except Exception as e:  # record the failure instead of stopping
            row["loaded"] = f"failed: {type(e).__name__}: {e}"[:200]
        rows.append(row)

inv = pd.DataFrame(rows)
cols = [
    "path",
    "set",
    "zenodo_record",
    "modality",
    "sample",
    "zone",
    "patient",
    "n_obs",
    "n_vars",
    "size_gb",
    "loaded",
    "note",
]
inv = inv.reindex(columns=cols)
inv.to_csv(OUT / "data_inventory.csv", index=False)
print(inv.groupby(["set", "loaded"], dropna=False).size().to_string())
