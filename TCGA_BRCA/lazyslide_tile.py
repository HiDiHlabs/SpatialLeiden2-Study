#!/usr/bin/env python3
"""Spot-centred tiling with lazyslide / wsidata for TCGA_BRCA / SpatialLeiden2 crops.

Crops (µm): spot 55, spot_margin 75, non_overlap ≈86.6, hest 112.
Tiles are read at TARGET_MPP (0.5 µm/px) so FOV is crop_um across slides.

Writes JPEG tiles and an AnnData with obs tile_path / tile_name (same layout
as the former colo_snake_tile / her2_snake_tile outputs).

Example:
  python TCGA_BRCA/lazyslide_tile.py \\
    -p Valdeolivas2024 -s ZEN48 --crop_size hest \\
    --tiles Valdeolivas2024_out/data/hest/tiles/ZEN48 \\
    --adata_out Valdeolivas2024_out/data/hest/tiles/ZEN48.h5ad
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from PIL import Image
from sklearn.neighbors import NearestNeighbors

# Visium hires / WSI-derived PNGs routinely exceed PIL's default ~179MP bomb limit.
Image.MAX_IMAGE_PIXELS = None

sys.path.insert(0, str(Path(__file__).resolve().parent))
from crop_sizes import CROP_UM, CROP_ORDER, TARGET_MPP, tile_px


def _read_visium_counts(visium: Path) -> ad.AnnData:
    """Load Visium counts + fullres coords; never open tissue_*_image.png.

    scanpy/squidpy read_visium require hires/lowres PNGs; tiling crops from
    WSI or hires separately and only needs the h5 + tissue_positions.
    """
    import scanpy as sc

    count_path = visium / "filtered_feature_bc_matrix.h5"
    if not count_path.is_file():
        raise SystemExit(f"missing {count_path}")
    adata = sc.read_10x_h5(str(count_path))
    adata.var_names_make_unique()
    xy = _positions_fullres(visium / "spatial", list(adata.obs_names))
    if xy is None:
        raise SystemExit(f"missing or incomplete tissue_positions under {visium / 'spatial'}")
    adata.obsm["spatial"] = xy
    return adata


def _load_visium(project: Path, sample: str) -> tuple[ad.AnnData, Path, float]:
    """Return AnnData (counts + labels), image path, slide mpp in image coords.

    Counts via 10x h5 + tissue_positions; tiles via WSI or hires PIL crop.
    """
    visium = project / "visium" / sample
    if not visium.is_dir():
        raise SystemExit(f"missing Visium folder {visium}")

    adata = _read_visium_counts(visium)

    anno_path = project / "Pathology_SpotAnnotations" / f"Pathologist_Annotations_{sample}.csv"
    if anno_path.is_file():
        anno = pd.read_csv(anno_path)
        anno = anno.rename(columns={anno.columns[0]: "barcode", anno.columns[1]: "label"})
        anno["barcode"] = anno["barcode"].astype(str)
        labels = anno.set_index("barcode")["label"]
        adata.obs["label"] = [str(labels.get(str(b), "") or "") for b in adata.obs_names]
    else:
        adata.obs["label"] = ""

    sf_path = visium / "spatial" / "scalefactors_json.json"
    sf = json.loads(sf_path.read_text())
    hires = visium / "spatial" / "tissue_hires_image.png"
    wsi_tif = project / "wsis" / f"{sample}.tif"
    wsi_tiff = project / "wsis" / f"{sample}.tiff"

    if wsi_tif.is_file() or wsi_tiff.is_file():
        image_path = wsi_tif if wsi_tif.is_file() else wsi_tiff
        xy = np.asarray(adata.obsm["spatial"], dtype=float)
        if "tissue_positions_list.csv" in {p.name for p in (visium / "spatial").iterdir()} or (
            visium / "spatial" / "tissue_positions.csv"
        ).is_file():
            full = _positions_fullres(visium / "spatial", list(adata.obs_names))
            if full is not None:
                xy = full
        if hires.is_file():
            xy = _maybe_scale_hires_to_wsi(xy, hires, image_path)
        adata.obsm["spatial"] = xy
        adata.obs["imagecol"] = xy[:, 0]
        adata.obs["imagerow"] = xy[:, 1]
        mpp = 55.0 / float(sf["spot_diameter_fullres"])
        return adata, image_path, mpp

    if not hires.is_file():
        raise SystemExit(f"no WSI under {project}/wsis/ and no {hires}")
    scale = float(sf.get("tissue_hires_scalef", 1.0))
    xy = np.asarray(adata.obsm["spatial"], dtype=float)
    # scanpy/squidpy store fullres spatial; convert to hires pixels when needed.
    if scale > 0 and abs(scale - 1.0) > 1e-6:
        xy = xy * scale
    adata.obsm["spatial"] = xy
    adata.obs["imagecol"] = xy[:, 0]
    adata.obs["imagerow"] = xy[:, 1]
    um_per_px_hires = 55.0 / (float(sf["spot_diameter_fullres"]) * scale)
    return adata, hires, um_per_px_hires


def _positions_fullres(spatial_dir: Path, barcodes: list[str]) -> np.ndarray | None:
    csv_list = spatial_dir / "tissue_positions_list.csv"
    csv_new = spatial_dir / "tissue_positions.csv"
    if csv_new.is_file():
        df = pd.read_csv(csv_new)
        cols = {c.lower(): c for c in df.columns}
        bcol = cols.get("barcode")
        xcol = cols.get("pxl_col_in_fullres")
        ycol = cols.get("pxl_row_in_fullres")
        if not (bcol and xcol and ycol):
            return None
        m = df.set_index(df[bcol].astype(str))
        rows = []
        for b in barcodes:
            if b not in m.index:
                return None
            rows.append([float(m.loc[b, xcol]), float(m.loc[b, ycol])])
        return np.asarray(rows, dtype=float)
    if csv_list.is_file():
        df = pd.read_csv(csv_list, header=None)
        # barcode, in_tissue, array_row, array_col, pxl_row, pxl_col
        df[0] = df[0].astype(str)
        m = df.set_index(0)
        rows = []
        for b in barcodes:
            if b not in m.index:
                return None
            rows.append([float(m.loc[b, 5]), float(m.loc[b, 4])])
        return np.asarray(rows, dtype=float)
    return None


def _maybe_scale_hires_to_wsi(xy: np.ndarray, hires: Path, wsi_path: Path) -> np.ndarray:
    from wsidata import open_wsi

    with Image.open(hires) as im:
        hires_w, hires_h = im.size
    wsi = open_wsi(str(wsi_path))
    shape = wsi.properties.shape
    # wsidata: shape is [height, width]
    wsi_h, wsi_w = int(shape[0]), int(shape[1])
    if hires_w <= 0 or hires_h <= 0:
        return xy
    sx = wsi_w / hires_w
    sy = wsi_h / hires_h
    # Already fullres if span is near WSI size.
    span_x = float(np.ptp(xy[:, 0])) if len(xy) else 0.0
    if span_x > 0.5 * wsi_w:
        return xy
    if abs(sx - 1.0) < 0.02 and abs(sy - 1.0) < 0.02:
        return xy
    return np.column_stack([xy[:, 0] * sx, xy[:, 1] * sy])


def _load_her2(project: Path, sample: str) -> tuple[ad.AnnData, Path, float]:
    """Old ST HER2 layout: spots×genes counts TSV + pixel coords (no stlearn).

    Andersson Zenodo matrices are [n_spots]×[n_genes] with obs ids like ``10x15``.
    Coords are joined the same way as stlearn.add.parsing (array x/y → pixel).
    """
    count_path = project / "count-matrices" / f"{sample}.tsv"
    spatial_path = project / "spot-selections" / f"{sample}_selection.tsv"
    image_path = project / "images" / "HE" / f"{sample}.jpg"
    meta_path = project / "meta" / f"{sample}_labeled_coordinates.tsv"
    for path in (count_path, spatial_path, image_path, meta_path):
        if not path.is_file():
            raise SystemExit(f"missing {path}")

    counts = pd.read_csv(count_path, sep="\t", index_col=0)
    # Already spots × genes (do not transpose).
    spatial = pd.read_csv(spatial_path, sep="\t")
    if "pixel_x" not in spatial.columns or "pixel_y" not in spatial.columns:
        raise SystemExit(f"HER2 {sample}: selection needs pixel_x / pixel_y columns")
    xcol = "x" if "x" in spatial.columns else "new_x"
    ycol = "y" if "y" in spatial.columns else "new_y"
    if xcol not in spatial.columns or ycol not in spatial.columns:
        raise SystemExit(f"HER2 {sample}: selection needs x/y (or new_x/new_y)")
    coord_map = {
        (int(r[xcol]), int(r[ycol])): (float(r["pixel_x"]), float(r["pixel_y"]))
        for _, r in spatial.iterrows()
    }

    imgcol: list[float] = []
    imgrow: list[float] = []
    keep: list[str] = []
    for spot_id in counts.index.astype(str):
        parts = spot_id.split("x")
        if len(parts) != 2:
            continue
        try:
            key = (int(parts[0]), int(parts[1]))
        except ValueError:
            continue
        if key not in coord_map:
            continue
        px, py = coord_map[key]
        keep.append(spot_id)
        imgcol.append(px)
        imgrow.append(py)
    if not keep:
        raise SystemExit(
            f"HER2 {sample}: no count rows matched selection "
            f"({counts.shape[0]} count rows, {len(spatial)} selection)"
        )
    adata = ad.AnnData(counts.loc[keep].copy())
    adata.obs["imagecol"] = np.asarray(imgcol, dtype=float)
    adata.obs["imagerow"] = np.asarray(imgrow, dtype=float)

    meta = pd.read_csv(meta_path, sep="\t").rename(columns={"x": "new_x", "y": "new_y"})
    labels_merged = pd.merge(spatial, meta, on=["new_x", "new_y"])
    # Both tables have pixel_*; pandas suffixes them (_x = selection).
    if "pixel_x_x" in labels_merged.columns:
        labels_merged = labels_merged.rename(
            columns={"pixel_x_x": "imagecol", "pixel_y_x": "imagerow"}
        )
    else:
        labels_merged = labels_merged.rename(
            columns={"pixel_x": "imagecol", "pixel_y": "imagerow"}
        )
    adata_obs_m = pd.merge(
        adata.obs, labels_merged, on=["imagecol", "imagerow"], how="left"
    )
    if len(adata_obs_m) != adata.n_obs:
        raise SystemExit(
            f"HER2 {sample}: label merge changed spot count "
            f"({adata.n_obs} -> {len(adata_obs_m)})"
        )
    adata.obs["label"] = adata_obs_m["label"].fillna("").astype(str).to_list()

    xy = np.column_stack(
        [adata.obs["imagecol"].astype(float), adata.obs["imagerow"].astype(float)]
    )
    adata.obsm["spatial"] = xy
    nbrs = NearestNeighbors(n_neighbors=2, algorithm="ball_tree").fit(xy)
    dist = nbrs.kneighbors(xy)[0]
    median_dist = float(np.median(dist[dist != 0]))
    # Old ST pitch is 200 µm.
    um_per_px = 200.0 / median_dist
    return adata, image_path, um_per_px


_WSI_SUFFIXES = {".tif", ".tiff", ".svs", ".ndpi", ".vms", ".vmu", ".scn", ".mrxs", ".bif"}


def _is_nullable_string(dtype) -> bool:
    text = str(dtype)
    return text == "string" or text.startswith("string[") or isinstance(dtype, pd.StringDtype)


def _as_plain_str_list(values) -> list[str]:
    out: list[str] = []
    for value in values:
        if value is None or (isinstance(value, float) and np.isnan(value)):
            out.append("")
        else:
            out.append(str(value))
    return out


def _frame_strings_as_object(df: pd.DataFrame) -> pd.DataFrame:
    """Rebuild a DataFrame without pd.StringArray / nullable-string columns."""
    cols: dict[str, object] = {}
    for col in df.columns:
        series = df[col]
        if (
            pd.api.types.is_string_dtype(series)
            or _is_nullable_string(series.dtype)
            or series.dtype == object
        ):
            vals = series.tolist()
            if series.dtype == object and not all(
                v is None or isinstance(v, (str, bytes)) or (isinstance(v, float) and np.isnan(v))
                for v in vals
            ):
                cols[col] = vals
            else:
                cols[col] = _as_plain_str_list(vals)
        else:
            cols[col] = series.to_numpy()
    index = pd.Index(_as_plain_str_list(df.index), name=df.index.name, dtype=object)
    out = pd.DataFrame(
        {
            k: pd.Series(v, index=index, dtype=object) if isinstance(v, list) else v
            for k, v in cols.items()
        },
        index=index,
    )
    for col in out.columns:
        if out[col].dtype != object and _is_nullable_string(out[col].dtype):
            out[col] = pd.Series(_as_plain_str_list(out[col].tolist()), index=index, dtype=object)
    return out


def _obs_strings_as_object(adata: ad.AnnData) -> None:
    """Rebuild obs/var without nullable strings (breaks anndata read < 0.11, e.g. gene_ids)."""
    adata.obs = _frame_strings_as_object(adata.obs)
    adata.var = _frame_strings_as_object(adata.var)


def _crop_window_pil(im: Image.Image, cx: int, cy: int, size: int) -> Image.Image:
    """Axis-aligned crop of ``size`` px centered on (cx, cy); pad with white if needed."""
    w, h = im.size
    half = size // 2
    x0, y0 = cx - half, cy - half
    x1, y1 = x0 + size, y0 + size
    out = Image.new("RGB", (size, size), (255, 255, 255))
    xa0, ya0 = max(0, x0), max(0, y0)
    xa1, ya1 = min(w, x1), min(h, y1)
    if xa0 >= xa1 or ya0 >= ya1:
        return out
    out.paste(im.crop((xa0, ya0, xa1, ya1)), (xa0 - x0, ya0 - y0))
    return out


def _tile_path_str(
    dest: Path,
    *,
    absolute: bool,
    project: str,
    sample: str,
    crop_name: str,
    name: str,
) -> str:
    if absolute:
        return str(dest.resolve())
    return f"{project}_out/data/{crop_name}/tiles/{sample}/{name}.jpeg"


def _write_tiles_raster(
    adata: ad.AnnData,
    image_path: Path,
    slide_mpp: float,
    crop_name: str,
    tiles_dir: Path,
    project: str,
    sample: str,
    *,
    absolute_paths: bool = False,
) -> tuple[list[str], list[str], np.ndarray, int, float]:
    """Crop from PNG/JPG (wsidata cannot open these). Match FOV, resize to TARGET_MPP."""
    crop_um = CROP_UM[crop_name]
    px_out = tile_px(crop_name, TARGET_MPP)
    px_src = max(1, int(round(crop_um / slide_mpp)))
    # PIL only — avoid np.asarray(full image), which doubles peak RAM on colo hires PNGs.
    im = Image.open(image_path)
    im = im.convert("RGB")
    xy = np.asarray(adata.obsm["spatial"], dtype=float)
    centres = np.rint(xy).astype(int)
    tiles_dir.mkdir(parents=True, exist_ok=True)
    names: list[str] = []
    paths: list[str] = []
    for col, row in centres:
        tile = _crop_window_pil(im, int(col), int(row), px_src).resize(
            (px_out, px_out), Image.Resampling.BILINEAR
        )
        name = f"{int(col)}-{int(row)}-{px_out}"
        dest = tiles_dir / f"{name}.jpeg"
        tile.save(dest, quality=95)
        names.append(name)
        paths.append(
            _tile_path_str(
                dest,
                absolute=absolute_paths,
                project=project,
                sample=sample,
                crop_name=crop_name,
                name=name,
            )
        )
    im.close()
    print(
        f"raster crop {image_path.name}: {px_src} px @ {slide_mpp:.4f} µm/px -> "
        f"{px_out} px @ {TARGET_MPP} µm/px",
        flush=True,
    )
    return names, paths, xy, px_out, crop_um


def _write_tiles_wsi(
    adata: ad.AnnData,
    image_path: Path,
    slide_mpp: float,
    crop_name: str,
    tiles_dir: Path,
    project: str,
    sample: str,
    *,
    absolute_paths: bool = False,
) -> tuple[list[str], list[str], np.ndarray, int, float]:
    from wsidata import TileSpec, open_wsi
    from wsidata.io import add_tiles

    px = tile_px(crop_name, TARGET_MPP)
    crop_um = CROP_UM[crop_name]
    wsi = open_wsi(str(image_path))
    tile_spec = TileSpec.from_wsidata(
        wsi,
        tile_px=px,
        mpp=TARGET_MPP,
        slide_mpp=slide_mpp,
    )
    xy = np.asarray(adata.obsm["spatial"], dtype=float)
    centres = np.rint(xy).astype(int)
    top_left = centres - px // 2
    add_tiles(
        wsi,
        key="spots",
        xys=top_left,
        tile_spec=tile_spec,
        tissue_ids=np.zeros(len(centres), dtype=int),
    )
    wsi.shapes["spots"]["barcode"] = np.asarray(adata.obs_names).astype(str)

    tiles_dir.mkdir(parents=True, exist_ok=True)
    dataset = wsi.ds.tile_images(tile_key="spots")
    names: list[str] = []
    paths: list[str] = []
    for i in range(len(dataset)):
        item = dataset[i]
        image = np.asarray(item["image"])
        if image.ndim == 2:
            image = np.stack([image, image, image], axis=-1)
        if image.shape[-1] == 4:
            image = image[..., :3]
        col, row = int(centres[i, 0]), int(centres[i, 1])
        name = f"{col}-{row}-{px}"
        dest = tiles_dir / f"{name}.jpeg"
        Image.fromarray(image.astype(np.uint8)).save(dest, quality=95)
        names.append(name)
        paths.append(
            _tile_path_str(
                dest,
                absolute=absolute_paths,
                project=project,
                sample=sample,
                crop_name=crop_name,
                name=name,
            )
        )
    return names, paths, xy, px, crop_um


def _write_tiles(
    adata: ad.AnnData,
    image_path: Path,
    slide_mpp: float,
    crop_name: str,
    tiles_dir: Path,
    adata_out: Path,
    project: str,
    sample: str,
    *,
    absolute_paths: bool = False,
) -> None:
    suffix = image_path.suffix.lower()
    if suffix in _WSI_SUFFIXES:
        names, paths, xy, px, crop_um = _write_tiles_wsi(
            adata,
            image_path,
            slide_mpp,
            crop_name,
            tiles_dir,
            project,
            sample,
            absolute_paths=absolute_paths,
        )
    else:
        # tissue_hires_image.png / HER2 .jpg: wsidata has no reader for these.
        names, paths, xy, px, crop_um = _write_tiles_raster(
            adata,
            image_path,
            slide_mpp,
            crop_name,
            tiles_dir,
            project,
            sample,
            absolute_paths=absolute_paths,
        )

    adata = adata.copy()
    adata.obs["imagecol"] = xy[:, 0].astype(np.float64)
    adata.obs["imagerow"] = xy[:, 1].astype(np.float64)
    adata.obs["tile_name"] = names
    adata.obs["tile_path"] = paths
    _obs_strings_as_object(adata)
    print("obs dtypes before write:", {c: str(adata.obs[c].dtype) for c in adata.obs.columns}, flush=True)
    print("obs index dtype:", adata.obs.index.dtype, type(adata.obs.index.dtype), flush=True)
    adata.uns["crop_size"] = crop_name
    adata.uns["crop_um"] = float(crop_um)
    adata.uns["tile_px"] = int(px)
    adata.uns["tile_mpp"] = float(TARGET_MPP)
    adata_out.parent.mkdir(parents=True, exist_ok=True)
    adata.write(adata_out)
    print(
        f"wrote {len(names)} tiles ({crop_name}: {crop_um:.1f} µm, "
        f"{px} px @ {TARGET_MPP} µm/px) -> {tiles_dir}",
        flush=True,
    )


def _load_coords_h5ad(coords_h5ad: Path, slide_mpp: float | None) -> tuple[ad.AnnData, float]:
    """Spot centres from an existing AnnData (e.g. DeepSpot-M atlas prediction)."""
    adata = ad.read_h5ad(coords_h5ad)
    if "spatial" not in adata.obsm:
        raise SystemExit(f"{coords_h5ad}: missing obsm['spatial']")
    if slide_mpp is not None:
        mpp = float(slide_mpp)
    else:
        he2st = adata.uns.get("he2st", {}) if isinstance(adata.uns, dict) else {}
        mpp = float(he2st.get("mpp", TARGET_MPP)) if isinstance(he2st, dict) else float(TARGET_MPP)
    xy = np.asarray(adata.obsm["spatial"], dtype=float)
    adata.obs["imagecol"] = xy[:, 0]
    adata.obs["imagerow"] = xy[:, 1]
    return adata, mpp


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("-p", default=None, help="Project directory (e.g. Valdeolivas2024 or HER2)")
    p.add_argument("-s", default=None, help="Sample id")
    p.add_argument("--crop_size", required=True, choices=CROP_ORDER)
    p.add_argument("--tiles", required=True, type=Path, help="Output tile directory")
    p.add_argument("--adata_out", required=True, type=Path)
    p.add_argument(
        "--layout",
        choices=("auto", "visium", "her2"),
        default="auto",
        help="auto: HER2 project name -> her2, else visium",
    )
    p.add_argument(
        "--wsi",
        type=Path,
        default=None,
        help="WSI path (SVS/TIFF); with --coords-h5ad skips project Visium/HER2 layout",
    )
    p.add_argument(
        "--coords-h5ad",
        type=Path,
        default=None,
        help="AnnData with obsm['spatial'] spot centres (DeepSpot-M atlas h5ad)",
    )
    p.add_argument(
        "--slide-mpp",
        type=float,
        default=None,
        help="Slide µm/px for --coords-h5ad (default: uns['he2st']['mpp'] or 0.5)",
    )
    p.add_argument(
        "--absolute-tile-paths",
        action="store_true",
        help="Write absolute paths in obs['tile_path'] (for non-project layouts)",
    )
    args = p.parse_args()

    if args.wsi is not None or args.coords_h5ad is not None:
        if args.wsi is None or args.coords_h5ad is None:
            raise SystemExit("Provide both --wsi and --coords-h5ad")
        if args.s is None:
            raise SystemExit("--s sample id is required with --wsi/--coords-h5ad")
        sample = args.s
        project_name = Path(args.p).name if args.p else "tcga"
        adata, slide_mpp = _load_coords_h5ad(Path(args.coords_h5ad), args.slide_mpp)
        image_path = Path(args.wsi)
        if not image_path.is_file():
            raise SystemExit(f"missing WSI: {image_path}")
        _write_tiles(
            adata,
            image_path,
            slide_mpp,
            args.crop_size,
            Path(args.tiles),
            Path(args.adata_out),
            project_name,
            sample,
            absolute_paths=True,
        )
        return

    if args.p is None or args.s is None:
        raise SystemExit("Provide -p/-s (project layout) or --wsi + --coords-h5ad")

    project = Path(args.p)
    sample = args.s
    layout = args.layout
    if layout == "auto":
        layout = "her2" if project.name == "HER2" or str(project) == "HER2" else "visium"

    if layout == "her2":
        adata, image_path, slide_mpp = _load_her2(project, sample)
    else:
        adata, image_path, slide_mpp = _load_visium(project, sample)

    _write_tiles(
        adata,
        image_path,
        slide_mpp,
        args.crop_size,
        Path(args.tiles),
        Path(args.adata_out),
        project.name if project.name else str(project),
        sample,
        absolute_paths=bool(args.absolute_tile_paths),
    )


if __name__ == "__main__":
    main()
