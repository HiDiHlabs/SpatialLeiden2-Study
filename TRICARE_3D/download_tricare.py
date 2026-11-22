#!/usr/bin/env python3
"""Download TRICARE 3D pathology H&E-like TIFFs + annotations → MuData.

Paper: Gao et al., Nat. Biomed. Eng. (2026)
  Deep-learning triage of three-dimensional pathology datasets for
  comprehensive and efficient pathologist assessments
Zenodo: https://doi.org/10.5281/zenodo.20052262 (~191 GB full release)

Each specimen depth-stack becomes one SpatialLeiden2-style .h5mu with modality
``HE`` (pathologist labels in obs, depth as spatial, ``tile_path`` for FM
extractors). There is no gene expression in this release, so Multileiden needs
an ``rna`` modality added separately (e.g. predicted ST) before
``spatialleiden_mudata_run.py``.

Usage
-----
  export TRICARE=./data/TRICARE

  conda env create -f TRICARE_3D/download_tricare.yml
  conda activate download_tricare

  # annotations + prostate test only (~7 GB 7z), then MuData
  python TRICARE_3D/download_tricare.py \\
    --out-dir "$TRICARE" --cohorts prostate_test

  # full prostate development (~60+ GB of zips)
  python TRICARE_3D/download_tricare.py \\
    --out-dir "$TRICARE" --cohorts prostate_development

  # CSVs + MuData from already-extracted images
  python TRICARE_3D/download_tricare.py \\
    --out-dir "$TRICARE" --skip-download --skip-extract

Extra esophagus development volumes (Zenodo storage limit) are on Google Drive;
see Zenodo record description. Pass those TIFFs under ``--images-dir`` if needed.

Next: ``TRICARE_3D/run_tricare_univ2_multileiden.py`` tiles planes → UNI2 → Multileiden
(HE + 3D spatial).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
import urllib.request
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

ZENODO_RECORD = "20052262"
ZENODO_API = f"https://zenodo.org/api/records/{ZENODO_RECORD}"
ZENODO_DOI = "10.5281/zenodo.20052262"

# Depth spacing used when packing neighboring levels (Zenodo description).
DEPTH_UM = {
    "prostate": 5.0,
    "esophagus": 3.0,
}

COHORTS = {
    "prostate_development": {
        "organ": "prostate",
        "split": "development",
        "csv": "prostate_development_cohort.csv",
        "archive_prefix": "prostate_development_cohort",
    },
    "prostate_test": {
        "organ": "prostate",
        "split": "test",
        "csv": "prostate_test_cohort.csv",
        "archive_prefix": "prostate_test_cohort",
    },
    "esophagus_development": {
        "organ": "esophagus",
        "split": "development",
        "csv": "esophagus_development_cohort.csv",
        "archive_prefix": "esophagus_development_cohort",
    },
    "esophagus_test": {
        "organ": "esophagus",
        "split": "test",
        "csv": "esophagus_test_cohort.csv",
        "archive_prefix": "esophagus_test_cohort",
    },
}

SLIDE_RE = re.compile(r"^(?P<specimen>.+)_(?P<depth>\d+)$")


def _log(msg: str) -> None:
    print(msg, flush=True)


def _urlopen(url: str, *, timeout: int = 120):
    req = urllib.request.Request(url, headers={"User-Agent": "tricare-download"})
    return urllib.request.urlopen(req, timeout=timeout)


def fetch_record() -> dict:
    with _urlopen(ZENODO_API, timeout=60) as r:
        return json.loads(r.read().decode("utf-8"))


def file_map(record: dict) -> dict[str, dict]:
    return {f["key"]: f for f in record.get("files", [])}


def download_file(
    meta: dict,
    dest: Path,
    *,
    skip_existing: bool = True,
) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    url = meta["links"]["self"]
    size = int(meta.get("size") or 0)
    md5 = (meta.get("checksum") or "").removeprefix("md5:")

    if skip_existing and dest.exists() and dest.stat().st_size > 0:
        if size and dest.stat().st_size == size:
            _log(f"  skip (exists): {dest.name}")
            return dest
        if size and dest.stat().st_size != size:
            _log(f"  incomplete {dest.name} ({dest.stat().st_size}/{size}); re-downloading")

    _log(f"  downloading {dest.name} ({size / 1e9:.2f} GB)" if size else f"  downloading {dest.name}")
    tmp = dest.with_suffix(dest.suffix + ".part")
    h = hashlib.md5()
    with _urlopen(url) as r, open(tmp, "wb") as out:
        while True:
            chunk = r.read(8 * 1024 * 1024)
            if not chunk:
                break
            out.write(chunk)
            if md5:
                h.update(chunk)
    if md5 and h.hexdigest() != md5:
        tmp.unlink(missing_ok=True)
        raise SystemExit(f"MD5 mismatch for {dest.name}: got {h.hexdigest()}, expected {md5}")
    tmp.replace(dest)
    return dest


def extract_archive(archive: Path, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    _log(f"  extracting {archive.name} -> {out_dir}")
    if archive.suffix.lower() == ".zip":
        with zipfile.ZipFile(archive) as zf:
            zf.extractall(out_dir)
        return
    if archive.suffix.lower() == ".7z":
        try:
            import py7zr
        except ImportError as e:
            raise SystemExit("Install py7zr (conda env TRICARE_3D/download_tricare.yml) for .7z archives") from e
        with py7zr.SevenZipFile(archive, mode="r") as zf:
            zf.extractall(path=out_dir)
        return
    raise SystemExit(f"Unsupported archive type: {archive}")


def parse_slide_id(slide_id: str) -> tuple[str, int]:
    m = SLIDE_RE.match(str(slide_id).strip())
    if not m:
        raise ValueError(f"Cannot parse slide_id: {slide_id}")
    return m.group("specimen"), int(m.group("depth"))


def index_tiffs(images_dir: Path) -> dict[str, Path]:
    """Map stem (e.g. PRD_P000_A_335) -> path."""
    out: dict[str, Path] = {}
    for p in images_dir.rglob("*"):
        if p.suffix.lower() not in {".tif", ".tiff"}:
            continue
        out[p.stem] = p.resolve()
    return out


def load_annotations(csv_path: Path, cohort_key: str) -> pd.DataFrame:
    meta = COHORTS[cohort_key]
    df = pd.read_csv(csv_path)
    need = {"case_id", "slide_id", "label"}
    if not need.issubset(df.columns):
        raise SystemExit(f"{csv_path.name} missing columns {need - set(df.columns)}")
    df = df.copy()
    df["slide_id"] = df["slide_id"].astype(str).str.strip()
    df["case_id"] = df["case_id"].astype(str).str.strip()
    df["pathologist_annotation"] = df["label"].astype(str).str.strip()
    parsed = df["slide_id"].map(parse_slide_id)
    df["specimen_id"] = [p[0] for p in parsed]
    df["depth"] = [p[1] for p in parsed]
    df["organ"] = meta["organ"]
    df["split"] = meta["split"]
    df["cohort"] = cohort_key
    return df


def build_specimen_mudata(
    specimen_id: str,
    levels: pd.DataFrame,
    image_index: dict[str, Path],
    out_h5mu: Path,
    *,
    images_root: Path,
) -> Path | None:
    import anndata as ad
    import mudata as md

    rows = levels.sort_values("depth").reset_index(drop=True)
    organ = str(rows["organ"].iloc[0])
    depth_um = float(DEPTH_UM[organ])
    paths = []
    keep = []
    for i, row in rows.iterrows():
        slide_id = row["slide_id"]
        path = image_index.get(slide_id)
        if path is None:
            continue
        paths.append(path)
        keep.append(i)
    if not keep:
        return None
    rows = rows.loc[keep].reset_index(drop=True)

    labels = rows["pathologist_annotation"].astype("string").fillna("").astype(str)
    obs = pd.DataFrame(
        {
            "slide_id": rows["slide_id"].astype(str).values,
            "specimen_id": specimen_id,
            "case_id": rows["case_id"].astype(str).values,
            "depth": rows["depth"].astype(int).values,
            "depth_um": rows["depth"].astype(float).values * depth_um,
            "pathologist_annotation": labels.values,
            "label": labels.values,
            "is_annotated": rows["is_annotated"].astype(bool).values,
            "organ": organ,
            "split": rows["split"].astype(str).values,
            "cohort": rows["cohort"].astype(str).values,
            "tile_path": [str(p) for p in paths],
        },
        index=pd.Index(rows["slide_id"].astype(str), name="spot"),
    )
    spatial = np.column_stack(
        [
            np.zeros(len(obs), dtype=np.float64),
            obs["depth_um"].to_numpy(dtype=np.float64),
        ]
    )
    # Placeholder until FM embeddings are written into HE.X / HE.obsm.
    he_X = np.zeros((len(obs), 1), dtype=np.float32)
    he = ad.AnnData(
        X=he_X,
        obs=obs.copy(),
        var=pd.DataFrame(index=pd.Index(["placeholder"], name="he_feature")),
    )
    he.obsm["spatial"] = spatial.copy()
    mdata = md.MuData({"HE": he})
    mdata.obs = obs.copy()
    mdata.obsm["spatial"] = spatial
    mdata.uns["tricare"] = {
        "doi": ZENODO_DOI,
        "zenodo_record": ZENODO_RECORD,
        "specimen_id": specimen_id,
        "organ": organ,
        "depth_interval_um": depth_um,
        "n_levels": int(len(obs)),
        "n_annotated": int(obs["is_annotated"].sum()),
        "images_root": str(images_root),
        "note": (
            "HE.X is a placeholder. Run an FM extractor on obs['tile_path'] "
            "and replace HE.X (or HE.obsm) before Multileiden; add mod['rna'] for gene+HE Multileiden."
        ),
    }
    out_h5mu.parent.mkdir(parents=True, exist_ok=True)
    mdata.write(out_h5mu)
    return out_h5mu


def expand_levels_with_neighbors(
    anno: pd.DataFrame,
    image_index: dict[str, Path],
) -> pd.DataFrame:
    """For each annotated specimen, include every on-disk depth (neighbors + annotated)."""
    frames = []
    for specimen_id, sub in anno.groupby("specimen_id", sort=False):
        label_by_slide = dict(zip(sub["slide_id"], sub["pathologist_annotation"]))
        case_id = sub["case_id"].iloc[0]
        organ = sub["organ"].iloc[0]
        split = sub["split"].iloc[0]
        cohort = sub["cohort"].iloc[0]
        prefix = specimen_id + "_"
        slides = []
        for stem, path in image_index.items():
            if not stem.startswith(prefix):
                continue
            try:
                sp, depth = parse_slide_id(stem)
            except ValueError:
                continue
            if sp != specimen_id:
                continue
            slides.append(
                {
                    "case_id": case_id,
                    "slide_id": stem,
                    "specimen_id": specimen_id,
                    "depth": depth,
                    "pathologist_annotation": label_by_slide.get(stem, pd.NA),
                    "is_annotated": stem in label_by_slide,
                    "organ": organ,
                    "split": split,
                    "cohort": cohort,
                    "image_path": str(path),
                }
            )
        if not slides:
            # Annotation row with no images yet — keep annotated ids for manifest.
            for _, row in sub.iterrows():
                slides.append(
                    {
                        "case_id": case_id,
                        "slide_id": row["slide_id"],
                        "specimen_id": specimen_id,
                        "depth": int(row["depth"]),
                        "pathologist_annotation": row["pathologist_annotation"],
                        "is_annotated": True,
                        "organ": organ,
                        "split": split,
                        "cohort": cohort,
                        "image_path": "",
                    }
                )
        frames.append(pd.DataFrame(slides))
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out-dir", type=Path, required=True, help="Output root (archives, images, mudata)")
    p.add_argument(
        "--cohorts",
        nargs="+",
        choices=list(COHORTS) + ["all"],
        default=["prostate_test"],
        help="Which Zenodo cohorts to fetch (default: prostate_test)",
    )
    p.add_argument("--skip-download", action="store_true", help="Use existing archives/CSVs under out-dir")
    p.add_argument("--skip-extract", action="store_true", help="Do not extract archives")
    p.add_argument(
        "--images-dir",
        type=Path,
        default=None,
        help="Extra TIFF root to index (merged with extracted images)",
    )
    p.add_argument("--max-specimens", type=int, default=None, help="Cap specimens written to MuData")
    p.add_argument("--list-files", action="store_true", help="List Zenodo files and exit")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    out = args.out_dir
    out.mkdir(parents=True, exist_ok=True)
    ann_dir = out / "annotations"
    arc_dir = out / "archives"
    img_dir = out / "images"
    mu_dir = out / "mudata"
    ann_dir.mkdir(parents=True, exist_ok=True)
    arc_dir.mkdir(parents=True, exist_ok=True)
    img_dir.mkdir(parents=True, exist_ok=True)

    cohorts = list(COHORTS) if "all" in args.cohorts else list(dict.fromkeys(args.cohorts))
    _log(f"TRICARE Zenodo {ZENODO_DOI}; cohorts={cohorts}")

    files: dict[str, dict] = {}
    if args.list_files or not args.skip_download:
        files = file_map(fetch_record())

    if args.list_files:
        for key in sorted(files):
            sz = files[key].get("size") or 0
            _log(f"{key}\t{sz / 1e9:.3f} GB")
        return

    # --- download CSVs + archives ---
    if not args.skip_download:
        for key in cohorts:
            meta = COHORTS[key]
            csv_key = meta["csv"]
            if csv_key not in files:
                raise SystemExit(f"Zenodo missing {csv_key}")
            download_file(files[csv_key], ann_dir / csv_key)
            prefix = meta["archive_prefix"]
            matches = [k for k in files if k.startswith(prefix) and k.endswith((".zip", ".7z"))]
            if not matches:
                raise SystemExit(f"No archives for prefix {prefix}")
            for k in sorted(matches):
                download_file(files[k], arc_dir / k)
    else:
        _log("skip-download: using local annotations/archives")

    # --- extract ---
    if not args.skip_extract:
        for key in cohorts:
            prefix = COHORTS[key]["archive_prefix"]
            dest = img_dir / key
            for archive in sorted(arc_dir.glob(f"{prefix}*")):
                if archive.suffix.lower() not in {".zip", ".7z"}:
                    continue
                # Avoid re-extracting if marker exists with files.
                marker = dest / f".extracted_{archive.stem}"
                if marker.exists():
                    _log(f"  skip extract (marker): {archive.name}")
                    continue
                extract_archive(archive, dest)
                marker.write_text(archive.name + "\n")
    else:
        _log("skip-extract: using existing images/")

    # --- index images ---
    image_index = index_tiffs(img_dir)
    if args.images_dir is not None:
        extra = index_tiffs(args.images_dir)
        image_index.update(extra)
        _log(f"Indexed {len(extra)} TIFFs from --images-dir")
    _log(f"Indexed {len(image_index)} TIFF levels under {img_dir}")

    # --- annotations + MuData ---
    written = []
    for key in cohorts:
        csv_path = ann_dir / COHORTS[key]["csv"]
        if not csv_path.exists():
            raise SystemExit(f"Missing annotation CSV: {csv_path}")
        anno = load_annotations(csv_path, key)
        levels = expand_levels_with_neighbors(anno, image_index)
        if levels.empty:
            _log(f"[{key}] no levels to write")
            continue
        levels.to_csv(out / f"{key}_levels.csv", index=False)
        specimens = list(dict.fromkeys(levels["specimen_id"].tolist()))
        if args.max_specimens is not None:
            specimens = specimens[: args.max_specimens]
        n_ok = 0
        for specimen_id in specimens:
            sub = levels[levels["specimen_id"] == specimen_id].copy()
            # Ensure annotation columns exist when image-only neighbors present
            if "is_annotated" not in sub.columns:
                sub["is_annotated"] = sub["pathologist_annotation"].notna()
            out_h5mu = mu_dir / key / f"{specimen_id}.h5mu"
            path = build_specimen_mudata(
                specimen_id,
                sub,
                image_index,
                out_h5mu,
                images_root=img_dir,
            )
            if path is None:
                _log(f"  skip {specimen_id}: no images on disk")
                continue
            n_ok += 1
            written.append(
                {
                    "cohort": key,
                    "specimen_id": specimen_id,
                    "h5mu": str(path),
                    "n_levels": int(len(sub[sub["slide_id"].isin(image_index)])),
                    "n_annotated": int(sub["is_annotated"].sum()) if "is_annotated" in sub else 0,
                }
            )
        _log(f"[{key}] wrote {n_ok}/{len(specimens)} MuData files -> {mu_dir / key}")

    manifest = pd.DataFrame(written)
    if not manifest.empty:
        manifest_path = out / "mudata_manifest.csv"
        manifest.to_csv(manifest_path, index=False)
        _log(f"Manifest: {manifest_path} ({len(manifest)} specimens)")
    meta_path = out / "tricare_download.json"
    meta_path.write_text(
        json.dumps(
            {
                "doi": ZENODO_DOI,
                "record": ZENODO_RECORD,
                "cohorts": cohorts,
                "depth_um": DEPTH_UM,
                "n_mudata": len(written),
            },
            indent=2,
        )
        + "\n"
    )
    _log(f"Done. MuData under {mu_dir}")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.exit(130)
