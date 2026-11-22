#!/usr/bin/env python3
"""Tile TRICARE 3D H&E planes → UNI2 features → SpatialLeiden2 Multileiden.

Input: ``TRICARE/mudata/{cohort}/{specimen}.h5mu`` from ``download_tricare.py``.

Multileiden follows HiDiHlabs/SpatialLeiden2-Study
(https://github.com/HiDiHlabs/SpatialLeiden2-Study): OpenST 3D pattern uses
``squidpy.gr.spatial_neighbors`` on 3D coords + ``distance2connectivity`` +
``spatialleiden.spatialleiden`` (latent UNI2 graph + spatial graph). VisiumHD
H&E tiles use the same ``spatialleiden`` API; multimodal ``spatialleiden_multimodal``
needs rna+HE and is not used here (TRICARE has no genes).

Envs (one per stage)
--------------------
  conda env create -f TRICARE_3D/download_tricare.yml              # download
  conda env create -f TRICARE_3D/tricare_tile.yml         # tile
  conda env create -f TRICARE_3D/tricare_univ2.yml        # univ2
  conda env create -f TRICARE_3D/tricare_multileiden.yml  # multileiden

Usage
-----
  export TRICARE=./data/TRICARE

  # local / interactive
  conda activate tricare_tile
  python TRICARE_3D/run_tricare_univ2_multileiden.py \\
    --tricare-dir "$TRICARE" --cohorts prostate_test --stage tile \\
    --tile-backend lazyslide --find-tissues

  # SLURM array (one specimen per task; tile=CPU, univ2=GPU)
  bash TRICARE_3D/submit_tricare_array.sh tile
  bash TRICARE_3D/submit_tricare_array.sh univ2   # after tile finishes

  conda activate tricare_multileiden
  python TRICARE_3D/run_tricare_univ2_multileiden.py \\
    --tricare-dir "$TRICARE" --cohorts prostate_test --stage multileiden
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image

Image.MAX_IMAGE_PIXELS = None

# Lateral µm/px from TRICARE / 3D CPath summaries (approx.).
DEFAULT_MPP = {
    "prostate": 1.0,
    "esophagus": 0.6,
}


def _log(msg: str) -> None:
    print(msg, flush=True)


def _discover_h5mu(tricare_dir: Path, cohorts: list[str]) -> list[Path]:
    root = tricare_dir / "mudata"
    paths: list[Path] = []
    for cohort in cohorts:
        d = root / cohort
        if not d.is_dir():
            raise SystemExit(f"Missing {d}; run download_tricare.py first")
        paths.extend(sorted(d.glob("*.h5mu")))
    return paths


def _organ_from_mdata(mdata) -> str:
    if "organ" in mdata.obs.columns and len(mdata.obs):
        return str(mdata.obs["organ"].iloc[0])
    uns = mdata.uns.get("tricare") or {}
    return str(uns.get("organ", "prostate"))


def _read_plane(path: Path) -> np.ndarray:
    """Return HxWx3 uint8 RGB."""
    try:
        import tifffile

        arr = tifffile.imread(path)
    except Exception:
        arr = np.asarray(Image.open(path))
    if arr.ndim == 2:
        arr = np.stack([arr, arr, arr], axis=-1)
    if arr.ndim == 3 and arr.shape[-1] > 3:
        arr = arr[..., :3]
    if arr.dtype != np.uint8:
        arr = np.clip(arr, 0, 255).astype(np.uint8) if arr.dtype != np.float32 else (np.clip(arr, 0, 1) * 255).astype(np.uint8)
    return arr


def _plane_meta(row: pd.Series, organ: str) -> tuple[int, float, str, bool]:
    depth = int(row.get("depth", 0))
    if "depth_um" in row.index and pd.notna(row["depth_um"]):
        z_um = float(row["depth_um"])
    else:
        z_step = 5.0 if organ == "prostate" else 3.0
        z_um = float(depth) * z_step
    label = str(row.get("pathologist_annotation", "") or "")
    is_ann = bool(row.get("is_annotated", bool(label)))
    return depth, z_um, label, is_ann


def _wsi_readable_tiff(plane_path: Path, cache_dir: Path) -> Path:
    """Rewrite flat/odd TRICARE TIFFs to RGB uint8 so wsidata/tiffslide can open them.

    Zenodo planes are often non-pyramidal; openslide-style readers then fail with
    \"Cannot open image ... using any of the readers\".
    """
    import tifffile

    cache_dir.mkdir(parents=True, exist_ok=True)
    dest = cache_dir / f"{plane_path.stem}_rgb.tiff"
    if dest.is_file() and dest.stat().st_size > 0:
        return dest
    arr = _read_plane(plane_path)
    if arr.ndim != 3 or arr.shape[-1] != 3:
        raise ValueError(f"expected HxWx3 RGB, got {arr.shape} for {plane_path}")
    # Contiguous RGB, no compression — readable by tiffslide / openslide readers
    tmp = dest.with_suffix(".tiff.part")
    tifffile.imwrite(
        tmp,
        np.ascontiguousarray(arr),
        photometric="rgb",
        metadata=None,
    )
    tmp.replace(dest)
    return dest


def _open_wsi_plane(plane_path: Path, cache_dir: Path):
    """open_wsi on the plane; on reader failure, convert to RGB TIFF and retry."""
    from wsidata import open_wsi

    try:
        return open_wsi(str(plane_path))
    except Exception as first:
        converted = _wsi_readable_tiff(plane_path, cache_dir)
        try:
            return open_wsi(str(converted))
        except Exception as second:
            raise RuntimeError(
                f"open_wsi failed on {plane_path.name} ({first}); "
                f"RGB rewrite {converted.name} also failed ({second})"
            ) from second


def _tile_plane_lazyslide(
    plane_path: Path,
    *,
    tile_dir: Path,
    cache_dir: Path,
    depth: int,
    z_um: float,
    mpp: float,
    tile_px: int,
    stride_px: int,
    find_tissues: bool,
    max_tiles_per_plane: int | None,
    spot: str,
    specimen_id: str,
    case_id: str,
    label: str,
    is_ann: bool,
    organ: str,
    cohort: str,
) -> list[dict]:
    """Non-overlapping tissue tiles via lazyslide (SpatialLeiden2 VisiumHD style).

    ``zs.pp.tile_tissues(..., stride_px=tile_px, edge=False)`` → no overlap, no
    partial edge tiles. See lazyslide docs: stride defaults to tile_px.
    """
    import lazyslide as zs

    wsi = _open_wsi_plane(plane_path, cache_dir)
    if find_tissues:
        try:
            zs.pp.find_tissues(wsi)
        except Exception as e:
            _log(f"  find_tissues failed ({e}); tiling full plane")
    # Non-overlap: stride_px == tile_px; edge=False drops incomplete border tiles
    zs.pp.tile_tissues(
        wsi,
        tile_px,
        stride_px=stride_px,
        edge=False,
        slide_mpp=mpp,
        mpp=mpp,
        background_filter=True,
        tissue_key="tissues" if find_tissues and "tissues" in getattr(wsi, "shapes", {}) else None,
        key_added="tiles",
    )
    if "tiles" not in wsi.shapes or len(wsi.shapes["tiles"]) == 0:
        return []

    dataset = wsi.ds.tile_images(tile_key="tiles")
    shapes = wsi.shapes["tiles"]
    rows: list[dict] = []
    n = len(dataset)
    if max_tiles_per_plane is not None:
        n = min(n, max_tiles_per_plane)
    for i in range(n):
        item = dataset[i]
        image = np.asarray(item["image"])
        if image.ndim == 2:
            image = np.stack([image, image, image], axis=-1)
        if image.shape[-1] == 4:
            image = image[..., :3]
        image = image.astype(np.uint8)
        # Centroid in slide/pixel coords from geometry
        geom = shapes.geometry.iloc[i]
        cx = float(geom.centroid.x)
        cy = float(geom.centroid.y)
        name = f"d{depth:04d}_y{int(round(cy)):05d}_x{int(round(cx)):05d}.jpeg"
        out_path = tile_dir / name
        if not out_path.exists():
            Image.fromarray(image).save(out_path, quality=90)
        rows.append(
            {
                "tile_name": name,
                "tile_path": str(out_path.resolve()),
                "plane_slide_id": str(spot),
                "specimen_id": specimen_id,
                "case_id": case_id,
                "depth": depth,
                "depth_um": z_um,
                "x_px": cx,
                "y_px": cy,
                "x_um": cx * mpp,
                "y_um": cy * mpp,
                "z_um": z_um,
                "pathologist_annotation": label,
                "is_annotated": is_ann,
                "organ": organ,
                "cohort": cohort,
            }
        )
    return rows


def _tile_plane_grid(
    plane_path: Path,
    *,
    tile_dir: Path,
    depth: int,
    z_um: float,
    mpp: float,
    tile_px: int,
    stride_px: int,
    min_mean: float,
    max_mean: float,
    max_tiles_per_plane: int | None,
    spot: str,
    specimen_id: str,
    case_id: str,
    label: str,
    is_ann: bool,
    organ: str,
    cohort: str,
) -> list[dict]:
    """Fallback non-overlapping grid (no padding; skip incomplete edge tiles)."""
    img = _read_plane(plane_path)
    h, w = img.shape[:2]
    # Only full tile_px×tile_px windows → no overlap when stride_px == tile_px
    ys = list(range(0, h - tile_px + 1, stride_px))
    xs = list(range(0, w - tile_px + 1, stride_px))
    rows: list[dict] = []
    kept = 0
    for y0 in ys:
        for x0 in xs:
            patch = img[y0 : y0 + tile_px, x0 : x0 + tile_px]
            mean = float(patch.mean())
            if mean < min_mean or mean > max_mean:
                continue
            name = f"d{depth:04d}_y{y0:05d}_x{x0:05d}.jpeg"
            out_path = tile_dir / name
            if not out_path.exists():
                Image.fromarray(patch).save(out_path, quality=90)
            cx = x0 + tile_px / 2.0
            cy = y0 + tile_px / 2.0
            rows.append(
                {
                    "tile_name": name,
                    "tile_path": str(out_path.resolve()),
                    "plane_slide_id": str(spot),
                    "specimen_id": specimen_id,
                    "case_id": case_id,
                    "depth": depth,
                    "depth_um": z_um,
                    "x_px": cx,
                    "y_px": cy,
                    "x_um": cx * mpp,
                    "y_um": cy * mpp,
                    "z_um": z_um,
                    "pathologist_annotation": label,
                    "is_annotated": is_ann,
                    "organ": organ,
                    "cohort": cohort,
                }
            )
            kept += 1
            if max_tiles_per_plane is not None and kept >= max_tiles_per_plane:
                return rows
    return rows


def tile_specimen(
    h5mu_path: Path,
    out_root: Path,
    *,
    tile_px: int,
    stride_px: int,
    backend: str,
    find_tissues: bool,
    min_mean: float,
    max_mean: float,
    max_tiles_per_plane: int | None,
) -> Path:
    import anndata as ad
    import mudata as md

    if stride_px != tile_px:
        _log(f"Warning: stride_px={stride_px} != tile_px={tile_px} → overlapping tiles")

    mdata = md.read(h5mu_path)
    organ = _organ_from_mdata(mdata)
    mpp = float(DEFAULT_MPP.get(organ, 1.0))
    specimen_id = h5mu_path.stem
    cohort = h5mu_path.parent.name
    tile_dir = out_root / "tiles" / cohort / specimen_id
    tile_dir.mkdir(parents=True, exist_ok=True)
    cache_dir = out_root / "wsi_cache" / cohort / specimen_id

    rows: list[dict] = []
    for spot, row in mdata.obs.iterrows():
        plane_path = Path(str(row["tile_path"]))
        if not plane_path.is_file():
            _log(f"  skip missing plane: {plane_path}")
            continue
        depth, z_um, label, is_ann = _plane_meta(row, organ)
        common = dict(
            tile_dir=tile_dir,
            depth=depth,
            z_um=z_um,
            mpp=mpp,
            tile_px=tile_px,
            stride_px=stride_px,
            max_tiles_per_plane=max_tiles_per_plane,
            spot=str(spot),
            specimen_id=specimen_id,
            case_id=str(row.get("case_id", "")),
            label=label,
            is_ann=is_ann,
            organ=organ,
            cohort=cohort,
        )
        used = backend
        plane_rows: list[dict] = []
        if backend == "lazyslide":
            try:
                plane_rows = _tile_plane_lazyslide(
                    plane_path,
                    find_tissues=find_tissues,
                    cache_dir=cache_dir,
                    **common,
                )
            except Exception as e:
                _log(f"  lazyslide failed on {plane_path.name} ({e}); falling back to grid")
                plane_rows = _tile_plane_grid(
                    plane_path, min_mean=min_mean, max_mean=max_mean, **common
                )
                used = "grid"
        else:
            plane_rows = _tile_plane_grid(
                plane_path, min_mean=min_mean, max_mean=max_mean, **common
            )
        rows.extend(plane_rows)
        _log(f"  {plane_path.name}: kept {len(plane_rows)} tiles ({used})")

    if not rows:
        raise SystemExit(f"No tiles for {h5mu_path}")

    obs = pd.DataFrame(rows)
    obs.index = pd.Index([f"{r['tile_name']}" for r in rows], name="spot")
    spatial = obs[["x_um", "y_um", "z_um"]].to_numpy(dtype=np.float64)
    adata = ad.AnnData(obs=obs)
    adata.obsm["spatial"] = spatial
    adata.uns["tricare"] = {
        "specimen_id": specimen_id,
        "cohort": cohort,
        "organ": organ,
        "mpp": mpp,
        "tile_px": tile_px,
        "stride_px": stride_px,
        "backend": backend,
        "non_overlapping": bool(stride_px == tile_px),
        "source_h5mu": str(h5mu_path),
    }
    out_h5ad = out_root / "tiles" / cohort / f"{specimen_id}.h5ad"
    out_h5ad.parent.mkdir(parents=True, exist_ok=True)
    adata.write(out_h5ad)
    _log(f"Wrote {out_h5ad} ({adata.n_obs} tiles)")
    return out_h5ad


def load_univ2(model_dir: str, device):
    import torch
    import timm
    from torchvision import transforms

    model_path = os.path.join(model_dir, "pytorch_model.bin")
    if not os.path.isfile(model_path):
        raise SystemExit(
            f"UNI2 weights not found at {model_path}. "
            "Pass --model-dir (directory with pytorch_model.bin) or set UNIV2_MODEL_DIR."
        )
    timm_kwargs = {
        "model_name": "vit_giant_patch14_224",
        "img_size": 224,
        "patch_size": 14,
        "depth": 24,
        "num_heads": 24,
        "init_values": 1e-5,
        "embed_dim": 1536,
        "mlp_ratio": 2.66667 * 2,
        "num_classes": 0,
        "no_embed_class": True,
        "mlp_layer": timm.layers.SwiGLUPacked,
        "act_layer": torch.nn.SiLU,
        "reg_tokens": 8,
        "dynamic_img_size": True,
    }
    model = timm.create_model(pretrained=False, **timm_kwargs)
    model.load_state_dict(torch.load(model_path, map_location="cpu"), strict=True)
    transform = transforms.Compose(
        [
            transforms.Resize((224, 224)),
            transforms.ToTensor(),
            transforms.Normalize(mean=(0.485, 0.456, 0.406), std=(0.229, 0.224, 0.225)),
        ]
    )
    model.to(device)
    model.eval()
    return model, transform


def run_univ2(
    tile_h5ad: Path,
    out_h5ad: Path,
    *,
    model,
    transform,
    device,
    batch_size: int,
) -> Path:
    import anndata as ad
    import torch
    from tqdm import tqdm

    adata = ad.read_h5ad(tile_h5ad)
    paths = adata.obs["tile_path"].tolist()
    feats = []
    for i in tqdm(range(0, len(paths), batch_size), desc=tile_h5ad.stem):
        batch_paths = paths[i : i + batch_size]
        imgs = []
        for p in batch_paths:
            imgs.append(transform(Image.open(p).convert("RGB")))
        batch = torch.stack(imgs).to(device)
        with torch.inference_mode():
            feats.append(model(batch).cpu().numpy())
    X = np.concatenate(feats, axis=0).astype(np.float32)
    out = ad.AnnData(X=X, obs=adata.obs.copy())
    out.obsm["spatial"] = np.asarray(adata.obsm["spatial"], dtype=np.float64)
    out.obsm["univ2"] = X
    out.uns = dict(adata.uns)
    out.uns["univ2"] = {"model": "UNI2", "n_features": int(X.shape[1])}
    out_h5ad.parent.mkdir(parents=True, exist_ok=True)
    out.write(out_h5ad)
    _log(f"Wrote {out_h5ad} HE={X.shape}")
    return out_h5ad


def run_multileiden(
    feature_h5ad: Path,
    out_h5mu: Path,
    *,
    n_neighs: int,
    n_neighs_spatial: int,
    resolution: float,
    layer_ratio: float,
    n_iterations: int,
    max_spatial_dist: float | None,
    random_state: int,
) -> Path:
    """SpatialLeiden2-Study style: UNI2 neighbors + 3D spatial → spatialleiden().

    Matches OpenST_3D_lymphnode.ipynb (3D generic neighbors + distance2connectivity
    + spatialleiden) and VisiumHD-CRC H&E tile SpatialLeiden.
    """
    import anndata as ad
    import mudata as md
    import scanpy as sc
    from spatialleiden import distance2connectivity, spatialleiden

    adata = ad.read_h5ad(feature_h5ad)
    X = np.asarray(adata.obsm["univ2"] if "univ2" in adata.obsm else adata.X, dtype=np.float32)
    if X.ndim != 2 or X.shape[0] != adata.n_obs:
        raise ValueError(f"{feature_h5ad}: bad UNI2 matrix shape {getattr(X, 'shape', None)} vs n_obs={adata.n_obs}")
    if "spatial" not in adata.obsm:
        raise ValueError(f"{feature_h5ad}: missing obsm['spatial']")
    spatial = np.asarray(adata.obsm["spatial"], dtype=np.float64)
    if spatial.shape != (adata.n_obs, 3):
        raise ValueError(f"{feature_h5ad}: expected spatial (n,3), got {spatial.shape}")
    if adata.n_obs < 3:
        raise ValueError(f"{feature_h5ad}: need >=3 tiles for Multileiden, got {adata.n_obs}")

    he = ad.AnnData(
        X=X,
        obs=adata.obs.copy(),
        var=pd.DataFrame(index=pd.Index([f"he_{i}" for i in range(X.shape[1])], name="he_feature")),
    )
    he.obsm["spatial"] = spatial
    he.obsm["X_he"] = X
    k_he = min(n_neighs, max(2, adata.n_obs - 1))
    k_sp = min(n_neighs_spatial, max(2, adata.n_obs - 1))
    # Latent (UNI2) knn — VisiumHD / OpenST pattern
    sc.pp.neighbors(he, use_rep="X_he", n_neighbors=k_he, random_state=random_state)

    # 3D spatial knn via scanpy (avoid squidpy→spatialdata→dask legacy error)
    sc.pp.neighbors(
        he,
        use_rep="spatial",
        n_neighbors=k_sp,
        key_added="spatial",
        random_state=random_state,
    )
    if "spatial_distances" in he.obsp:
        if max_spatial_dist is not None and max_spatial_dist > 0:
            # OpenST_3D: zero out edges beyond a µm cutoff before connectivity
            he.obsp["spatial_distances"].data[
                he.obsp["spatial_distances"].data > max_spatial_dist
            ] = 0
            he.obsp["spatial_distances"].eliminate_zeros()
        he.obsp["spatial_connectivities"] = distance2connectivity(he.obsp["spatial_distances"])

    # OpenST_3D: directed=(False, True), layer_ratio>1 emphasizes spatial
    spatialleiden(
        he,
        resolution=resolution,
        directed=(False, True),
        layer_ratio=layer_ratio,
        n_iterations=n_iterations,
        random_state=random_state,
        key_added="spatialleiden",
    )
    he.obs["multileiden"] = he.obs["spatialleiden"]

    mdata = md.MuData({"HE": he})
    mdata.obs = he.obs.copy()
    mdata.obsm["spatial"] = spatial
    mdata.obsp["spatial_connectivities"] = he.obsp["spatial_connectivities"]
    if "spatial_distances" in he.obsp:
        mdata.obsp["spatial_distances"] = he.obsp["spatial_distances"]
    mdata.uns["tricare"] = dict(adata.uns.get("tricare") or {})
    mdata.uns["multileiden"] = {
        "api": "spatialleiden.spatialleiden",
        "study": "HiDiHlabs/SpatialLeiden2-Study",
        "mode": "HE_latent+spatial_3d",
        "resolution": resolution,
        "layer_ratio": layer_ratio,
        "n_iterations": n_iterations,
        "n_neighs": n_neighs,
        "n_neighs_spatial": n_neighs_spatial,
        "max_spatial_dist": max_spatial_dist,
        "n_clusters": int(he.obs["spatialleiden"].nunique()),
        "note": "No rna modality; SpatialLeiden2 H&E/3D path (not spatialleiden_multimodal).",
    }
    out_h5mu.parent.mkdir(parents=True, exist_ok=True)
    mdata.write(out_h5mu)
    _log(f"Wrote {out_h5mu}; clusters={mdata.uns['multileiden']['n_clusters']}")
    return out_h5mu


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--tricare-dir", type=Path, required=True, help="Root from download_tricare.py")
    p.add_argument("--out-dir", type=Path, default=None, help="Default: {tricare-dir}/univ2_multileiden")
    p.add_argument(
        "--cohorts",
        nargs="+",
        default=["prostate_test"],
        help="Cohort folder names under mudata/ (default: prostate_test)",
    )
    p.add_argument(
        "--stage",
        choices=["all", "tile", "univ2", "multileiden"],
        default="all",
        help="Pipeline stage (default: all)",
    )
    p.add_argument(
        "--specimen-index",
        type=int,
        default=None,
        help="Process only h5mus[i] (for SLURM_ARRAY_TASK_ID)",
    )
    p.add_argument(
        "--specimen",
        type=str,
        default=None,
        help="Process only this specimen stem (e.g. PRT_P002)",
    )
    p.add_argument(
        "--list-specimens",
        action="store_true",
        help="Print mudata paths (one per line) and exit",
    )
    p.add_argument(
        "--skip-existing",
        action="store_true",
        help="Skip specimens whose stage output already exists",
    )
    p.add_argument("--max-specimens", type=int, default=None)
    p.add_argument("--tile-px", type=int, default=224, help="Tile size in pixels (UNI2=224)")
    p.add_argument(
        "--stride-px",
        type=int,
        default=None,
        help="Default: tile-px (non-overlapping). Smaller → overlap.",
    )
    p.add_argument(
        "--tile-backend",
        choices=["lazyslide", "grid"],
        default="lazyslide",
        help="lazyslide: zs.pp.tile_tissues (SpatialLeiden2); grid: numpy fallback",
    )
    p.add_argument(
        "--find-tissues",
        action="store_true",
        help="Run zs.pp.find_tissues before tiling (lazyslide only)",
    )
    p.add_argument("--min-mean", type=float, default=15.0, help="grid backend: skip dark tiles")
    p.add_argument("--max-mean", type=float, default=240.0, help="grid backend: skip blank tiles")
    p.add_argument("--max-tiles-per-plane", type=int, default=None, help="Cap tiles per depth (debug)")
    p.add_argument("--model-dir", default=None, help="UNI2 dir with pytorch_model.bin")
    p.add_argument("--device", default="auto", choices=("auto", "cuda", "cpu"))
    p.add_argument("--batch-size", type=int, default=32)
    p.add_argument("--n-neighs", type=int, default=15, help="UNI2 latent knn (sc.pp.neighbors)")
    p.add_argument(
        "--n-neighs-spatial",
        type=int,
        default=10,
        help="3D spatial knn (OpenST_3D default 10)",
    )
    p.add_argument(
        "--resolution",
        type=float,
        default=0.7,
        help="spatialleiden resolution (OpenST_3D default 0.7)",
    )
    p.add_argument(
        "--layer-ratio",
        type=float,
        default=2.0,
        help="spatialleiden layer_ratio; higher → more spatial (OpenST_3D default 2)",
    )
    p.add_argument("--n-iterations", type=int, default=10, help="Leiden iterations (OpenST_3D: 10)")
    p.add_argument(
        "--max-spatial-dist",
        type=float,
        default=None,
        help="Drop 3D spatial edges beyond this µm (OpenST_3D uses a cutoff; default: keep all knn)",
    )
    p.add_argument("--random-state", type=int, default=42)
    return p.parse_args()


def _select_h5mus(h5mus: list[Path], args: argparse.Namespace) -> list[Path]:
    if args.max_specimens is not None:
        h5mus = h5mus[: args.max_specimens]
    if args.specimen is not None:
        want = str(args.specimen)
        h5mus = [p for p in h5mus if p.stem == want]
        if not h5mus:
            raise SystemExit(f"No specimen matching --specimen {want!r}")
    if args.specimen_index is not None:
        i = int(args.specimen_index)
        if i < 0 or i >= len(h5mus):
            raise SystemExit(f"--specimen-index {i} out of range 0..{len(h5mus) - 1}")
        h5mus = [h5mus[i]]
    return h5mus


def main() -> None:
    args = parse_args()
    tricare = args.tricare_dir
    out = args.out_dir or (tricare / "univ2_multileiden")
    out.mkdir(parents=True, exist_ok=True)
    stride = args.stride_px if args.stride_px is not None else args.tile_px
    stage = args.stage

    h5mus = _discover_h5mu(tricare, args.cohorts)
    if args.list_specimens:
        if args.max_specimens is not None:
            h5mus = h5mus[: args.max_specimens]
        for p in h5mus:
            print(p)
        return

    h5mus = _select_h5mus(h5mus, args)
    _log(f"Specimens: {len(h5mus)}; out={out}; stage={stage}")

    # --- tile ---
    tile_h5ads: list[Path] = []
    if stage in ("all", "tile"):
        for h5mu in h5mus:
            tile_out = out / "tiles" / h5mu.parent.name / f"{h5mu.stem}.h5ad"
            if args.skip_existing and tile_out.is_file():
                _log(f"[tile] skip existing {tile_out}")
                tile_h5ads.append(tile_out)
                continue
            _log(f"[tile] {h5mu}")
            tile_h5ads.append(
                tile_specimen(
                    h5mu,
                    out,
                    tile_px=args.tile_px,
                    stride_px=stride,
                    backend=args.tile_backend,
                    find_tissues=args.find_tissues,
                    min_mean=args.min_mean,
                    max_mean=args.max_mean,
                    max_tiles_per_plane=args.max_tiles_per_plane,
                )
            )
    elif stage != "multileiden":
        for h5mu in h5mus:
            p = out / "tiles" / h5mu.parent.name / f"{h5mu.stem}.h5ad"
            if p.is_file():
                tile_h5ads.append(p)
            else:
                _log(f"  missing tiles: {p}")

    # --- univ2 ---
    feat_h5ads: list[Path] = []
    if stage in ("all", "univ2"):
        import torch

        if args.device == "auto":
            device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        else:
            device = torch.device(args.device)
        model_dir = args.model_dir or os.environ.get("UNIV2_MODEL_DIR", "models")
        _log(f"[univ2] device={device} model_dir={model_dir}")
        model, transform = load_univ2(model_dir, device)
        for tile_h5ad in tile_h5ads:
            cohort = tile_h5ad.parent.name
            feat_out = out / "features" / cohort / f"{tile_h5ad.stem}_univ2.h5ad"
            if args.skip_existing and feat_out.is_file():
                _log(f"[univ2] skip existing {feat_out}")
                feat_h5ads.append(feat_out)
                continue
            _log(f"[univ2] {tile_h5ad.name}")
            feat_h5ads.append(
                run_univ2(
                    tile_h5ad,
                    feat_out,
                    model=model,
                    transform=transform,
                    device=device,
                    batch_size=max(1, args.batch_size),
                )
            )
    elif stage != "multileiden":
        for tile_h5ad in tile_h5ads:
            p = out / "features" / tile_h5ad.parent.name / f"{tile_h5ad.stem}_univ2.h5ad"
            if p.is_file():
                feat_h5ads.append(p)

    # --- multileiden ---
    if stage in ("all", "multileiden"):
        # Discover features directly (do not require tile h5ads to still be present)
        if stage == "multileiden" or not feat_h5ads:
            feat_h5ads = []
            for cohort in args.cohorts:
                feat_h5ads.extend(sorted((out / "features" / cohort).glob("*_univ2.h5ad")))
            if args.specimen is not None or args.specimen_index is not None:
                stems = {p.stem for p in h5mus}
                feat_h5ads = [f for f in feat_h5ads if f.stem.replace("_univ2", "") in stems]
        if not feat_h5ads:
            raise SystemExit(f"No *_univ2.h5ad under {out / 'features'}")
        _log(f"Multileiden inputs: {len(feat_h5ads)} feature files")
        n_ok = 0
        n_fail = 0
        for feat in feat_h5ads:
            cohort = feat.parent.name
            stem = feat.stem.replace("_univ2", "")
            out_h5mu = out / "multileiden" / cohort / f"{stem}.h5mu"
            if args.skip_existing and out_h5mu.is_file():
                _log(f"[multileiden] skip existing {out_h5mu}")
                n_ok += 1
                continue
            _log(f"[multileiden] {feat.name}")
            try:
                run_multileiden(
                    feat,
                    out_h5mu,
                    n_neighs=args.n_neighs,
                    n_neighs_spatial=args.n_neighs_spatial,
                    resolution=args.resolution,
                    layer_ratio=args.layer_ratio,
                    n_iterations=args.n_iterations,
                    max_spatial_dist=args.max_spatial_dist,
                    random_state=args.random_state,
                )
                n_ok += 1
            except Exception as e:
                n_fail += 1
                _log(f"[multileiden] FAILED {feat.name}: {type(e).__name__}: {e}")
        _log(f"Multileiden done: ok={n_ok} fail={n_fail}")
        if n_fail and not n_ok:
            raise SystemExit(1)

    if args.specimen_index is None and args.specimen is None:
        meta = {
            "tricare_dir": str(tricare),
            "out_dir": str(out),
            "cohorts": args.cohorts,
            "stage": stage,
            "tile_px": args.tile_px,
            "stride_px": stride,
            "n_specimens": len(h5mus),
        }
        (out / "run_meta.json").write_text(json.dumps(meta, indent=2) + "\n")
    _log("Done.")

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.exit(130)
