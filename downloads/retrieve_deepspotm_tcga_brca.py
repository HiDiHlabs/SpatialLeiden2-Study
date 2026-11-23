#!/usr/bin/env python3
"""Retrieve DeepSpot-M TCGA-BRCA predictions for MultiLeiden (MuData).

Selects slides that (1) have TCGA-CDR OS and/or PFI (PFS-like) endpoints,
(2) are primary-tumor diagnostic (DX / FFPE) specimens at the diagnosis
timepoint used by those endpoints, and (3) by default are AJCC pathologic
Stage II or III (sublevels IIA/IIB/… and IIIA/IIIB/IIIC collapsed to II / III).
By default keeps all atlas genes (no PCC filter). Optional ``--pcc-filter`` /
``--pcc-cameo`` subsets to PCC > threshold. Assembles SpatialLeiden2-style
MuData with rna + HE embeddings.

PCC cohort (default: CAMEO-Breast)
----------------------------------
DeepSpot-M was trained on a pan-cancer *Visium* corpus (15 datasets / ~500
slides; Nonchev et al., medRxiv 10.64898/2026.06.19.26356060). Xenium,
including the public 10x breast slides packaged as theislab/CAMEO-Breast
(DOI 10.57967/hf/7909), was used only for transfer / test-time adaptation,
not Visium pretraining — so CAMEO is the preferred held-out PCC source.
Optional ``--pcc-visium-h5ad`` remains for a Visium-tile PCC path (may overlap
Visium training cohorts such as Andersson HER2).

Data
----
- Predictions: huggingface.co/datasets/ratschlab/TCGA_virtual_spatial_transcriptomics_atlas
  (gated; huggingface-cli login after accepting terms)
- Clinical: TCGA-CDR Supplemental Table S1 (Liu et al., Cell 2018)
  https://gdc.cancer.gov/about-data/publications/PanCan-Clinical-2018
- HE: GDC open-access diagnostic SVS (gdc_file_id in atlas metadata). Preferred
  SLURM path: ``--stage fetch`` then ``lazyslide_tile.py --crop_size non_overlap``
  (override with ``CROP_SIZE``) then ``univ2_snake.py`` (envs
  ``download_deepspotm_tcga`` / ``lazyslide`` / ``univ2``).
  Legacy ``--stage he`` still OpenSlide-crops in-process. Embeddings →
  ``he_embeddings/{sample_id}.npy`` → MuData modality ``HE``.
  While the SVS is open, also write ``he_overview/{sample_id}.jpg`` (WSI
  thumbnail for Multileiden plots).
- PCC: theislab/CAMEO-Breast niches (224×224 + spot_gexp) by default; or
  measured Visium BRCA tiles via --pcc-visium-h5ad

STP-Bench (NEXGEM/STP-Bench) uses PCC >= 0.4 as “well-predicted”; this script
uses the same threshold by default. CAMEO scores the 280-gene Xenium Human
Breast v1 panel only (genes outside the panel cannot be PCC-filtered this way).

Usage
-----
  conda env create -f downloads/download_deepspotm_tcga.yml
  conda activate download_deepspotm_tcga
  huggingface-cli login

  # Stage II/III + 1:1 PFI event:censored + UNI2 HE → MuData
  python downloads/retrieve_deepspotm_tcga_brca.py \\
    --out-dir ./data/TCGA_BRCA \\
    --download-svs \\
    --he-encoder univ2 \\
    --univ2-model-dir models \\
    --stages II,III \\
    --enrich-pfi 1:1 \\
    --max-patients 100

  # assemble only (predictions + HE npy already present)
  python downloads/retrieve_deepspotm_tcga_brca.py \\
    --out-dir ./data/TCGA_BRCA --stage mudata
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

# Sibling download helpers + TCGA_BRCA Multileiden builders.
_SCRIPTS_DIR = Path(__file__).resolve().parent
_REPO_ROOT = _SCRIPTS_DIR.parent
_TCGA_DIR = _REPO_ROOT / "TCGA_BRCA"
for _p in (_SCRIPTS_DIR, _TCGA_DIR):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

ATLAS_REPO = "ratschlab/TCGA_virtual_spatial_transcriptomics_atlas"
CAMEO_BREAST_REPO = "theislab/CAMEO-Breast"
# Official 10x Xenium Human Breast Gene Expression Panel (280 genes); column
# order of CAMEO spot_gexp / gexp matches this panel metadata row order.
XENIUM_HBREAST_V1_CSV = (
    "https://cdn.10xgenomics.com/raw/upload/v1684952547/software-support/"
    "Xenium-panels/hBreast_panel_files/Xenium_hBreast_v1_metadata.csv"
)
CDR_URL = "https://api.gdc.cancer.gov/data/1b5f413e-a8d1-4d10-92eb-7c4ae739ed81"
CDR_FILENAME = "TCGA-CDR-SupplementalTableS1.xlsx"
GDC_DATA = "https://api.gdc.cancer.gov/data"

# Primary solid tumor at diagnosis; DX = diagnostic FFPE slide (atlas FFPE).
PRIMARY_SAMPLE_TYPES = {"01"}
DIAGNOSTIC_SLIDE_RE = re.compile(r"-DX\d*\.", re.IGNORECASE)
# AJCC pathologic stage: Roman numeral after optional "Stage" (IIA → II, IIIC → III).
_AJCC_STAGE_RE = re.compile(
    r"^\s*(?:stage\s*)?(iv|iii|ii|i)(?:[abc]|is)?\b",
    re.IGNORECASE,
)
_AJCC_MISSING = {
    "",
    "nan",
    "none",
    "na",
    "[not available]",
    "[not applicable]",
    "[discrepancy]",
    "[unknown]",
    "stage x",
}


def collapse_ajcc_stage(raw) -> str | None:
    """Collapse AJCC pathologic stage sublevels to I / II / III / IV.

    Examples: ``Stage IIA`` → ``II``, ``IIIB`` → ``III``, ``Stage IV`` → ``IV``.
    Returns ``None`` for missing / unparsable values.
    """
    if raw is None or (isinstance(raw, float) and np.isnan(raw)):
        return None
    s = str(raw).strip()
    if s.lower() in _AJCC_MISSING:
        return None
    m = _AJCC_STAGE_RE.match(s)
    if not m:
        return None
    return m.group(1).upper()


def _log(msg: str) -> None:
    print(msg, flush=True)


def _download(url: str, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    _log(f"Downloading {url} -> {dest}")
    urllib.request.urlretrieve(url, dest)
    return dest


def load_tcga_cdr(path: Path) -> pd.DataFrame:
    """Load TCGA-CDR; keep BRCA rows with OS and/or PFI (PFS-like) endpoints."""
    # Force openpyxl: system xlrd>=2 cannot read .xlsx and pandas may pick it first.
    try:
        df = pd.read_excel(path, sheet_name=0, engine="openpyxl")
    except ImportError as e:
        raise SystemExit(
            "Reading TCGA-CDR .xlsx requires openpyxl. "
            "pip/conda install openpyxl (see downloads/download_deepspotm_tcga.yml)"
        ) from e
    # Column names vary slightly across mirrors; normalize.
    lower = {c.lower(): c for c in df.columns}
    type_col = lower.get("type") or lower.get("cancer type")
    barcode_col = lower.get("bcr_patient_barcode")
    if type_col is None or barcode_col is None:
        raise SystemExit(f"CDR missing type/barcode columns: {list(df.columns)}")
    df = df.rename(columns={barcode_col: "case_id", type_col: "type"})
    df["case_id"] = df["case_id"].astype(str).str.strip()
    brca = df[df["type"].astype(str).str.upper().eq("BRCA")].copy()

    def _col(*names: str) -> str | None:
        for n in names:
            if n in brca.columns:
                return n
            if n.lower() in lower and lower[n.lower()] in brca.columns:
                return lower[n.lower()]
        return None

    os_e = _col("OS")
    os_t = _col("OS.time", "OS_time")
    pfi_e = _col("PFI", "PFS")
    pfi_t = _col("PFI.time", "PFI_time", "PFS.time", "PFS_time")
    stage_c = _col("ajcc_pathologic_tumor_stage", "ajcc_pathologic_stage", "pathologic_stage")
    if os_e is None or os_t is None or pfi_e is None or pfi_t is None:
        raise SystemExit(f"CDR BRCA missing OS/PFI columns: {list(brca.columns)}")
    if stage_c is None:
        raise SystemExit(f"CDR BRCA missing AJCC stage column: {list(brca.columns)}")

    brca["OS"] = pd.to_numeric(brca[os_e], errors="coerce")
    brca["OS_time"] = pd.to_numeric(brca[os_t], errors="coerce")
    brca["PFI"] = pd.to_numeric(brca[pfi_e], errors="coerce")
    brca["PFI_time"] = pd.to_numeric(brca[pfi_t], errors="coerce")
    brca["ajcc_pathologic_tumor_stage"] = brca[stage_c].astype(str)
    brca["ajcc_stage"] = brca["ajcc_pathologic_tumor_stage"].map(collapse_ajcc_stage)
    has_os = brca["OS"].notna() & brca["OS_time"].notna()
    has_pfi = brca["PFI"].notna() & brca["PFI_time"].notna()
    brca = brca[has_os | has_pfi].copy()
    brca["has_OS"] = has_os.loc[brca.index]
    brca["has_PFI"] = has_pfi.loc[brca.index]
    return brca


def parse_slide_barcode(sample_id: str) -> dict:
    """Parse TCGA slide file stem for sample type and DX vs frozen section."""
    # e.g. TCGA-A1-A0SB-01Z-00-DX1.<uuid>
    stem = str(sample_id)
    parts = stem.split("-")
    case_id = "-".join(parts[:3]) if len(parts) >= 3 else stem[:12]
    sample_type = ""
    if len(parts) >= 4 and len(parts[3]) >= 2:
        sample_type = parts[3][:2]
    is_dx = bool(DIAGNOSTIC_SLIDE_RE.search(stem)) or ("-DX" in stem.upper())
    is_frozen_section = any(tag in stem.upper() for tag in ("-TS", "-BS", "-MS"))
    return {
        "case_id": case_id,
        "sample_type_code": sample_type,
        "is_diagnostic_dx": is_dx,
        "is_frozen_section": is_frozen_section,
    }


def select_endpoint_slides(meta: pd.DataFrame, cdr: pd.DataFrame) -> pd.DataFrame:
    """Primary-tumor DX (FFPE) slides joined to CDR OS/PFI at diagnosis timepoint."""
    meta = meta.copy()
    if "dataset" in meta.columns:
        meta = meta[meta["dataset"].astype(str).str.upper().isin(["TCGA_BRCA", "BRCA"])]
    parsed = meta["sample_id"].map(parse_slide_barcode).apply(pd.Series)
    for c in parsed.columns:
        if c not in meta.columns:
            meta[c] = parsed[c]
        else:
            # Prefer metadata.csv fields when present; fill gaps from parse.
            meta[c] = meta[c].fillna(parsed[c])

    if "case_id" not in meta.columns:
        meta["case_id"] = parsed["case_id"]
    meta["case_id"] = meta["case_id"].astype(str).str.strip()

    # Atlas: FFPE == DX diagnostic; FF == TS/BS/MS frozen. Endpoints are from
    # diagnosis → keep primary tumor (01) diagnostic FFPE/DX only.
    st = meta.get("sample_type_code", parsed["sample_type_code"]).astype(str).str.zfill(2)
    meta["sample_type_code"] = st
    slide_type = meta.get("slide_type", pd.Series(["FFPE"] * len(meta), index=meta.index))
    is_dx = meta["is_diagnostic_dx"].fillna(False).astype(bool)
    is_ffpe = slide_type.astype(str).str.upper().eq("FFPE") | is_dx
    is_primary = meta["sample_type_code"].isin(PRIMARY_SAMPLE_TYPES)
    keep = meta[is_ffpe & is_primary].copy()

    cdr_cols = [
        "OS",
        "OS_time",
        "PFI",
        "PFI_time",
        "has_OS",
        "has_PFI",
        "ajcc_pathologic_tumor_stage",
        "ajcc_stage",
    ]
    cdr_idx = cdr.set_index("case_id")
    keep = keep.join(cdr_idx[cdr_cols], on="case_id", how="inner")
    keep["endpoint_timepoint"] = "diagnosis"
    keep["endpoint_note"] = (
        "DX/FFPE primary-tumor slide at pathologic diagnosis; "
        "OS/PFI from TCGA-CDR (PFI is the recommended BRCA PFS-like endpoint); "
        "ajcc_stage collapses AJCC sublevels (e.g. IIA→II, IIIC→III)"
    )
    return keep.reset_index(drop=True)


def parse_stages_arg(raw: str | None) -> list[str] | None:
    """Parse ``--stages``: comma-separated I/II/III/IV, or ``all`` / empty to disable."""
    if raw is None:
        return ["II", "III"]
    s = str(raw).strip()
    if not s or s.lower() in {"all", "any", "*"}:
        return None
    allowed = {"I", "II", "III", "IV"}
    stages: list[str] = []
    for tok in s.split(","):
        t = tok.strip().upper()
        if not t:
            continue
        # Accept arabic 2/3 as aliases for II/III.
        t = {"1": "I", "2": "II", "3": "III", "4": "IV"}.get(t, t)
        if t not in allowed:
            raise SystemExit(f"Invalid --stages token {tok!r}; use I,II,III,IV or all")
        if t not in stages:
            stages.append(t)
    return stages or None


def filter_by_ajcc_stage(selected: pd.DataFrame, stages: list[str] | None) -> pd.DataFrame:
    """Keep slides whose collapsed ``ajcc_stage`` is in ``stages`` (None = no filter)."""
    if stages is None:
        return selected
    if "ajcc_stage" not in selected.columns:
        raise SystemExit("selected slides missing ajcc_stage; re-run clinical select")
    before = len(selected)
    keep = selected[selected["ajcc_stage"].isin(stages)].copy()
    _log(
        f"AJCC stage filter {stages}: {before} → {len(keep)} slides "
        f"(raw stage counts among kept: "
        f"{keep['ajcc_pathologic_tumor_stage'].value_counts().to_dict() if len(keep) else {}})"
    )
    return keep.reset_index(drop=True)


def parse_enrich_ratio(raw: str) -> tuple[int, int]:
    """Parse ``events:censored`` ratio (e.g. ``1:1`` → (1, 1))."""
    parts = str(raw).strip().split(":")
    if len(parts) != 2:
        raise SystemExit(f"Invalid --enrich-pfi {raw!r}; expected events:censored like 1:1")
    try:
        n_event, n_cens = int(parts[0]), int(parts[1])
    except ValueError as e:
        raise SystemExit(f"Invalid --enrich-pfi {raw!r}; integers required") from e
    if n_event < 1 or n_cens < 1:
        raise SystemExit(f"Invalid --enrich-pfi {raw!r}; both sides must be >= 1")
    return n_event, n_cens


def one_slide_per_patient(selected: pd.DataFrame) -> pd.DataFrame:
    """Keep the first DX slide per ``case_id`` (assumes sorted by case_id, sample_id)."""
    if "case_id" not in selected.columns:
        raise SystemExit("selected slides missing case_id")
    return selected.drop_duplicates(subset=["case_id"], keep="first").reset_index(drop=True)


def enrich_pfi_balanced(
    selected: pd.DataFrame,
    *,
    ratio: tuple[int, int] = (1, 1),
    max_patients: int | None = None,
    seed: int = 0,
) -> pd.DataFrame:
    """Patient-level PFI event:censored enrichment (one slide per case).

    Takes as many complete ratio blocks as fit under ``max_patients`` and the
    available event / censored pools. Example ``1:1`` with ``max_patients=100``
    → up to 50 events + 50 censored.
    """
    if "PFI" not in selected.columns or "has_PFI" not in selected.columns:
        raise SystemExit("enrich-pfi requires PFI columns; re-run clinical select")
    need = selected["has_PFI"].fillna(False).astype(bool) & selected["PFI"].isin([0, 1, 0.0, 1.0])
    pool = one_slide_per_patient(selected.loc[need].copy())
    pool["PFI"] = pd.to_numeric(pool["PFI"], errors="coerce")
    pool = pool[pool["PFI"].isin([0, 1])].copy()
    pool["pfi_event"] = pool["PFI"].astype(int)
    events = pool[pool["pfi_event"] == 1]
    censored = pool[pool["pfi_event"] == 0]
    r_e, r_c = ratio
    # Largest k such that k*r_e events and k*r_c censored fit pools and max_patients.
    k_by_pool = min(len(events) // r_e, len(censored) // r_c)
    if max_patients is not None:
        if max_patients < r_e + r_c:
            raise SystemExit(
                f"--max-patients {max_patients} too small for ratio {r_e}:{r_c} "
                f"(need at least {r_e + r_c})"
            )
        k_by_cap = max_patients // (r_e + r_c)
        k = min(k_by_pool, k_by_cap)
    else:
        k = k_by_pool
    n_e, n_c = k * r_e, k * r_c
    if n_e == 0 or n_c == 0:
        raise SystemExit(
            f"PFI enrich {r_e}:{r_c} yielded empty set "
            f"(events={len(events)}, censored={len(censored)}, max_patients={max_patients})"
        )
    rng = np.random.default_rng(seed)
    ev_idx = rng.choice(events.index.to_numpy(), size=n_e, replace=False)
    ce_idx = rng.choice(censored.index.to_numpy(), size=n_c, replace=False)
    out = pd.concat([events.loc[ev_idx], censored.loc[ce_idx]], axis=0)
    out["enrich_arm"] = np.where(out["pfi_event"] == 1, "event", "censored")
    out = out.sort_values(["pfi_event", "case_id", "sample_id"], ascending=[False, True, True])
    _log(
        f"PFI enrich {r_e}:{r_c} (seed={seed}): "
        f"{len(events)} events / {len(censored)} censored available → "
        f"kept {n_e} events + {n_c} censored = {len(out)} patients "
        f"(max_patients={max_patients})"
    )
    return out.reset_index(drop=True)


def download_atlas_metadata(out_dir: Path) -> pd.DataFrame:
    from huggingface_hub import hf_hub_download

    path = hf_hub_download(ATLAS_REPO, "metadata.csv", repo_type="dataset", local_dir=out_dir / "hf_cache")
    return pd.read_csv(path)


def download_slides(file_paths: list[str], out_dir: Path) -> None:
    from huggingface_hub import snapshot_download

    if not file_paths:
        raise SystemExit("No atlas file_paths to download")
    _log(f"Downloading {len(file_paths)} atlas slides…")
    snapshot_download(
        ATLAS_REPO,
        repo_type="dataset",
        local_dir=str(out_dir / "atlas"),
        allow_patterns=list(file_paths),
    )


def _is_tiff_like(path: Path, *, min_bytes: int = 1_000_000) -> bool:
    """Aperio SVS is TIFF-based; reject tiny/HTML/JSON error bodies from GDC."""
    if not path.is_file() or path.stat().st_size < min_bytes:
        return False
    with path.open("rb") as f:
        magic = f.read(4)
    return magic in (b"II*\x00", b"MM\x00*")


def save_he_overview(
    svs_path: Path,
    out_path: Path,
    *,
    max_side: int = 2048,
    quality: int = 85,
) -> Path:
    """Write a JPEG WSI overview (openslide thumbnail) for Multileiden plots."""
    try:
        import openslide
    except ImportError as e:
        raise SystemExit("openslide-python required for H&E overview thumbnails") from e

    slide = openslide.OpenSlide(str(svs_path))
    try:
        w, h = slide.dimensions
        scale = min(1.0, float(max_side) / float(max(w, h)))
        size = (max(1, int(round(w * scale))), max(1, int(round(h * scale))))
        thumb = slide.get_thumbnail(size).convert("RGB")
    finally:
        slide.close()
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    thumb.save(out_path, format="JPEG", quality=int(quality), optimize=True)
    _log(f"Wrote H&E overview {out_path} ({thumb.size[0]}x{thumb.size[1]})")
    return out_path


def he_overview_path(out_dir: Path, sample_id: str) -> Path:
    return Path(out_dir) / "he_overview" / f"{sample_id}.jpg"


def download_gdc_svs(file_id: str, dest: Path, *, retries: int = 2) -> Path:
    """Download an open-access GDC slide; validate TIFF magic before returning.

    Set ``GDC_TOKEN`` if the file is controlled-access.
    """
    import os

    dest.parent.mkdir(parents=True, exist_ok=True)
    if _is_tiff_like(dest):
        return dest
    if dest.exists():
        _log(f"Removing invalid/partial SVS ({dest.stat().st_size} bytes): {dest}")
        dest.unlink()

    url = f"{GDC_DATA}/{file_id}"
    token = os.environ.get("GDC_TOKEN") or os.environ.get("GDC_ACCESS_TOKEN")
    headers = {}
    if token:
        headers["X-Auth-Token"] = token.strip()

    last_err: Exception | None = None
    for attempt in range(1, retries + 1):
        tmp = dest.with_suffix(dest.suffix + f".part{attempt}")
        if tmp.exists():
            tmp.unlink()
        _log(f"GDC SVS {file_id} -> {dest} (attempt {attempt}/{retries})")
        try:
            req = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(req, timeout=600) as resp, tmp.open("wb") as out:
                ctype = (resp.headers.get("Content-Type") or "").lower()
                if "json" in ctype or "html" in ctype or "text/" in ctype:
                    body = resp.read(2000)
                    raise RuntimeError(
                        f"GDC returned {ctype!r} instead of a slide "
                        f"(first bytes: {body[:200]!r}). "
                        "Check gdc_file_id / open-access status; set GDC_TOKEN if needed."
                    )
                while True:
                    chunk = resp.read(8 * 1024 * 1024)
                    if not chunk:
                        break
                    out.write(chunk)
            if not _is_tiff_like(tmp):
                size = tmp.stat().st_size if tmp.exists() else 0
                tmp.unlink(missing_ok=True)
                raise RuntimeError(
                    f"Downloaded file is not a TIFF/SVS (size={size}). "
                    "Likely a truncated or error payload from GDC."
                )
            tmp.replace(dest)
            return dest
        except Exception as e:
            last_err = e
            tmp.unlink(missing_ok=True)
            _log(f"GDC download failed: {e}")
    raise SystemExit(f"Failed to download valid SVS for {file_id}: {last_err}")


def load_xenium_breast_gene_symbols(cache_csv: Path | None = None) -> list[str]:
    """Gene symbols for the 280-gene Xenium Human Breast v1 panel (CAMEO column order)."""
    if cache_csv is not None and cache_csv.exists() and cache_csv.stat().st_size > 0:
        df = pd.read_csv(cache_csv)
    else:
        dest = cache_csv or Path("_xenium_hbreast_v1_metadata.csv")
        _download(XENIUM_HBREAST_V1_CSV, dest)
        df = pd.read_csv(dest)
    gene_col = "Genes" if "Genes" in df.columns else df.columns[0]
    genes = df[gene_col].astype(str).str.strip().tolist()
    if len(genes) != 280:
        raise SystemExit(f"Expected 280 Xenium breast panel genes, got {len(genes)}")
    return genes


def compute_gene_pcc_cameo_breast(
    *,
    pcc_threshold: float,
    out_csv: Path,
    device: str = "cpu",
    max_niches: int | None = 2000,
    gene_panel_csv: Path | None = None,
    hf_cache_dir: Path | None = None,
) -> pd.DataFrame:
    """DeepSpot-M vs CAMEO-Breast niche pseudobulk; gene-wise PCC on the 280-gene panel.

    Streams ``theislab/CAMEO-Breast`` (avoids full ~40 GB download). Ground truth is
    ``spot_gexp`` log1p-CP10k; predictions use the niche H&E ``image`` (already 224×224).
    """
    import torch
    from scipy.stats import pearsonr

    try:
        from datasets import load_dataset
    except ImportError as e:
        raise SystemExit(
            "datasets is required for --pcc-cameo. "
            "conda env update -f downloads/download_deepspotm_tcga.yml"
        ) from e
    try:
        from deepspotm import DeepSpotM
    except ImportError as e:
        raise SystemExit(
            "deepspotm is required for --pcc-cameo. pip install deepspotm"
        ) from e

    panel_cache = gene_panel_csv or (out_csv.parent / "Xenium_hBreast_v1_metadata.csv")
    genes = load_xenium_breast_gene_symbols(panel_cache)
    model, image_processor = DeepSpotM.from_pretrained("ratschlab/DeepSpotM", source="scgpt")
    model.eval()
    model.to(device)
    model_genes = set(model.gene_names)
    shared = [g for g in genes if g in model_genes]
    if not shared:
        raise SystemExit("No shared gene symbols between CAMEO breast panel and DeepSpot-M")
    shared_idx = [genes.index(g) for g in shared]
    _log(
        f"CAMEO-Breast PCC: {len(shared)}/{len(genes)} panel genes in DeepSpot-M; "
        f"max_niches={max_niches}"
    )

    ds = load_dataset(
        CAMEO_BREAST_REPO,
        split="train",
        streaming=True,
        cache_dir=str(hf_cache_dir) if hf_cache_dir else None,
    )
    # Reservoir sample niches so we do not require a full local copy.
    n_keep = max_niches if max_niches is not None else None
    reservoir: list[tuple[object, np.ndarray]] = []
    rng = np.random.default_rng(42)
    seen = 0
    for ex in ds:
        img = ex["image"]
        y = np.asarray(ex["spot_gexp"], dtype=np.float64).reshape(-1)
        if y.size != len(genes):
            raise SystemExit(
                f"CAMEO spot_gexp length {y.size} != panel {len(genes)}; "
                "check gene-order assumption (Xenium Human Breast v1)"
            )
        if n_keep is None:
            reservoir.append((img, y))
        else:
            if len(reservoir) < n_keep:
                reservoir.append((img, y))
            else:
                j = int(rng.integers(0, seen + 1))
                if j < n_keep:
                    reservoir[j] = (img, y)
        seen += 1
        if n_keep is not None and seen % 5000 == 0:
            _log(f"  streamed {seen} niches (reservoir={len(reservoir)})")

    if not reservoir:
        raise SystemExit(f"No niches streamed from {CAMEO_BREAST_REPO}")
    _log(f"Scoring DeepSpot-M on {len(reservoir)} / {seen} CAMEO niches")

    preds = np.full((len(reservoir), len(shared)), np.nan, dtype=np.float64)
    Y = np.zeros((len(reservoir), len(shared)), dtype=np.float64)
    with torch.no_grad():
        for row, (img, y_counts) in enumerate(reservoir):
            lib = float(y_counts.sum())
            if lib <= 0:
                continue
            y_norm = np.log1p(y_counts / lib * 1e4)
            Y[row, :] = y_norm[shared_idx]
            tile = image_processor(img.convert("RGB")).unsqueeze(0).to(device)
            vals = model.predict_genes(tile, shared)
            preds[row, :] = vals.squeeze(0).detach().cpu().numpy()

    rows = []
    for j, g in enumerate(shared):
        x = preds[:, j]
        y = Y[:, j]
        mask = np.isfinite(x) & np.isfinite(y)
        if mask.sum() < 30:
            pcc = np.nan
        else:
            pcc, _ = pearsonr(x[mask], y[mask])
        rows.append({"gene": g, "pcc": pcc, "n_spots": int(mask.sum()), "source": "cameo_breast"})
    pcc_df = pd.DataFrame(rows).sort_values("pcc", ascending=False)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    pcc_df.to_csv(out_csv, index=False)
    keep = pcc_df[pcc_df["pcc"] > pcc_threshold]
    _log(f"CAMEO PCC>{pcc_threshold}: {len(keep)} / {len(pcc_df)} genes -> {out_csv}")
    return pcc_df


def compute_gene_pcc_deepspotm(
    visium_h5ad: Path,
    tile_dir: Path | None,
    *,
    pcc_threshold: float,
    out_csv: Path,
    device: str = "cpu",
    max_spots: int | None = 2000,
) -> pd.DataFrame:
    """Predict DeepSpot-M on Visium BRCA tiles; gene-wise PCC vs measured expression."""
    import anndata as ad
    import scanpy as sc
    import torch
    from PIL import Image
    from scipy.stats import pearsonr

    try:
        from deepspotm import DeepSpotM
    except ImportError as e:
        raise SystemExit(
            "deepspotm is required for --pcc-visium-h5ad. "
            "pip install deepspotm  (or use --gene-pcc-csv)"
        ) from e

    adata = sc.read_h5ad(visium_h5ad)
    X_true = adata.layers["counts"] if "counts" in adata.layers else adata.X
    if hasattr(X_true, "toarray"):
        X_true = X_true.toarray()
    X_true = np.asarray(X_true, dtype=np.float64)
    # Match DeepSpot-M scale: log1p CP10k
    lib = X_true.sum(axis=1, keepdims=True)
    lib[lib == 0] = 1.0
    Y = np.log1p(X_true / lib * 1e4)

    model, image_processor = DeepSpotM.from_pretrained("ratschlab/DeepSpotM", source="scgpt")
    model.eval()
    model.to(device)

    genes_model = list(model.gene_names)
    gene_to_idx = {g: i for i, g in enumerate(genes_model)}
    shared = [g for g in adata.var_names.astype(str) if g in gene_to_idx]
    if not shared:
        raise SystemExit("No shared gene symbols between Visium AnnData and DeepSpot-M")

    n = adata.n_obs if max_spots is None else min(adata.n_obs, max_spots)
    rng = np.random.default_rng(42)
    spot_idx = np.sort(rng.choice(adata.n_obs, size=n, replace=False)) if n < adata.n_obs else np.arange(adata.n_obs)

    preds = np.full((len(spot_idx), len(shared)), np.nan, dtype=np.float64)
    obs_names = adata.obs_names.astype(str)

    def _tile_path(i: int) -> Path | None:
        if tile_dir is None:
            return None
        name = obs_names[i]
        for cand in (
            tile_dir / f"{name}.jpg",
            tile_dir / f"{name}.jpeg",
            tile_dir / f"{name}.png",
            tile_dir / f"{name}.tif",
        ):
            if cand.exists():
                return cand
        # barcode without sample prefix
        short = name.split("_")[-1]
        for cand in tile_dir.glob(f"*{short}*.jpg"):
            return cand
        return None

    with torch.no_grad():
        for row, i in enumerate(spot_idx):
            path = _tile_path(int(i))
            if path is None:
                continue
            img = Image.open(path).convert("RGB").resize((224, 224))
            tile = image_processor(img).unsqueeze(0).to(device)
            vals = model.predict_genes(tile, shared)
            preds[row, :] = vals.squeeze(0).detach().cpu().numpy()

    true_idx = [list(adata.var_names.astype(str)).index(g) for g in shared]
    Y_sub = Y[spot_idx][:, true_idx]
    rows = []
    for j, g in enumerate(shared):
        x = preds[:, j]
        y = Y_sub[:, j]
        mask = np.isfinite(x) & np.isfinite(y)
        if mask.sum() < 30:
            pcc = np.nan
        else:
            pcc, _ = pearsonr(x[mask], y[mask])
        rows.append({"gene": g, "pcc": pcc, "n_spots": int(mask.sum())})
    pcc_df = pd.DataFrame(rows).sort_values("pcc", ascending=False)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    pcc_df.to_csv(out_csv, index=False)
    keep = pcc_df[pcc_df["pcc"] > pcc_threshold]
    _log(f"PCC>={pcc_threshold}: {len(keep)} / {len(pcc_df)} genes -> {out_csv}")
    return pcc_df


def _load_univ2_model(model_dir: Path, device: str):
    """Load UNI2 (vit_giant_patch14_224) from ``pytorch_model.bin`` — same as univ2_snake.py."""
    import os

    import timm
    import torch
    from torchvision import transforms

    model_dir = Path(model_dir or os.environ.get("UNIV2_MODEL_DIR", "models"))
    model_path = model_dir / "pytorch_model.bin"
    if not model_path.is_file():
        raise SystemExit(
            f"UNI2 weights not found at {model_path}. "
            "Pass --univ2-model-dir or set UNIV2_MODEL_DIR to the directory with pytorch_model.bin."
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
    if device == "auto":
        device_t = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    else:
        device_t = torch.device(device)
        if device_t.type == "cuda" and not torch.cuda.is_available():
            raise SystemExit(
                "Requested --he-device cuda but this PyTorch build has no CUDA "
                f"(torch={torch.__version__}). Use the univ2 conda env: "
                "export CONDA_ENV_HE=hest  (or install CUDA torch in the HE env)."
            )
    model.to(device_t)
    model.eval()
    return model, transform, device_t


def extract_he_embeddings_from_svs(
    adata_path: Path,
    svs_path: Path,
    out_npy: Path,
    *,
    encoder: str = "univ2",
    univ2_model_dir: Path | None = None,
    device: str = "auto",
    batch_size: int = 32,
) -> np.ndarray:
    """Crop 224 px patches at spot centres from the GDC SVS and embed.

    encoder:
      - univ2: Mahmood UNI2 (default; same setup as univ2_snake.py)
      - deepspotm_midnight: pooled Midnight backbone via deepspotm
      - mean_rgb: 3-d mean RGB (debug / smoke only)
    """
    import anndata as ad

    try:
        import openslide
    except ImportError as e:
        raise SystemExit("openslide-python required for HE from SVS") from e

    adata = ad.read_h5ad(adata_path)
    coords = np.asarray(adata.obsm["spatial"], dtype=np.float64)
    # atlas: spatial is (x=col, y=row) in level-0 pixels; patch 224 @ 0.5 µm/px
    mpp = float(adata.uns.get("he2st", {}).get("mpp", 0.5))
    patch_um = float(adata.uns.get("he2st", {}).get("patch_um", 112.0))
    level0_patch = int(round(patch_um / mpp))

    slide = openslide.OpenSlide(str(svs_path))
    patches = []
    for x, y in coords:
        left = int(round(x - level0_patch / 2))
        top = int(round(y - level0_patch / 2))
        region = slide.read_region((left, top), 0, (level0_patch, level0_patch)).convert("RGB")
        region = region.resize((224, 224))
        patches.append(region)
    slide.close()

    if encoder == "mean_rgb":
        emb = np.stack([np.asarray(p, dtype=np.float32).reshape(-1, 3).mean(0) for p in patches])
    elif encoder in ("univ2", "uni2"):
        import torch

        model, transform, device_t = _load_univ2_model(univ2_model_dir or Path("models"), device)
        bs = max(1, int(batch_size))
        vecs = []
        with torch.inference_mode():
            for i in range(0, len(patches), bs):
                batch = torch.stack([transform(p) for p in patches[i : i + bs]]).to(device_t)
                feat = model(batch)
                vecs.append(feat.detach().cpu())
        emb = torch.cat(vecs, dim=0).numpy().astype(np.float32)
    elif encoder == "deepspotm_midnight":
        import torch
        from deepspotm import DeepSpotM

        model, image_processor = DeepSpotM.from_pretrained("ratschlab/DeepSpotM", source="scgpt")
        model.eval()
        encode = getattr(model, "encode_image", None) or getattr(model, "get_image_features", None)
        backbone = getattr(model, "backbone", None) or getattr(model, "image_encoder", None)

        def _embed_one(t: "torch.Tensor") -> np.ndarray:
            if encode is not None:
                v = encode(t)
                if isinstance(v, (tuple, list)):
                    v = v[0]
                return v.squeeze(0).detach().cpu().numpy().reshape(-1)
            if backbone is not None:
                out = backbone(t)
                if isinstance(out, (tuple, list)):
                    out = out[0]
                # token grid -> mean pool
                if out.ndim == 3:
                    out = out.mean(dim=1)
                return out.squeeze(0).detach().cpu().numpy().reshape(-1)
            raise SystemExit(
                "deepspotm build exposes neither encode_image/get_image_features nor backbone; "
                "precompute HE embeddings and pass --he-npy-dir"
            )

        vecs = []
        with torch.no_grad():
            for p in patches:
                t = image_processor(p).unsqueeze(0)
                vecs.append(_embed_one(t))
        emb = np.stack(vecs).astype(np.float32)
    else:
        raise SystemExit(f"Unknown --he-encoder {encoder}")

    out_npy.parent.mkdir(parents=True, exist_ok=True)
    np.save(out_npy, emb)
    return emb


def assemble_mudata(
    adata_path: Path,
    he_npy: Path,
    gene_keep: list[str] | None,
    clinical_row: pd.Series,
    out_h5mu: Path,
) -> None:
    import anndata as ad

    try:
        from multileiden_mudata import build_gene_he_mudata
    except ImportError as e:
        raise SystemExit(
            f"Cannot import multileiden_mudata from {_TCGA_DIR} "
            "(needed only for --stage mudata / all). "
            f"Ensure TCGA_BRCA/multileiden_mudata.py is present. ({e})"
        ) from e

    adata = ad.read_h5ad(adata_path)
    adata.X = np.asarray(adata.X, dtype=np.float32)
    if gene_keep:
        keep = [g for g in gene_keep if g in set(adata.var_names.astype(str))]
        if not keep:
            raise SystemExit(f"No requested genes found in {adata_path}")
        adata = adata[:, keep].copy()
    n_genes = int(adata.n_vars)
    he = np.load(he_npy)
    if he.shape[0] != adata.n_obs:
        raise SystemExit(f"HE rows {he.shape[0]} != spots {adata.n_obs} for {adata_path}")

    obs = adata.obs.copy()
    for col in (
        "OS",
        "OS_time",
        "PFI",
        "PFI_time",
        "case_id",
        "sample_id",
        "endpoint_timepoint",
        "ajcc_pathologic_tumor_stage",
        "ajcc_stage",
        "pfi_event",
        "enrich_arm",
    ):
        if col in clinical_row.index:
            obs[col] = clinical_row[col]

    mdata = build_gene_he_mudata(
        gene_X=adata.X,
        gene_var_names=adata.var_names,
        he_X=he,
        spatial=adata.obsm["spatial"],
        obs=obs,
        obs_names=adata.obs_names,
    )
    clin_keys = [
        k
        for k in (
            "OS",
            "OS_time",
            "PFI",
            "PFI_time",
            "case_id",
            "sample_id",
            "endpoint_note",
            "ajcc_pathologic_tumor_stage",
            "ajcc_stage",
            "pfi_event",
            "enrich_arm",
        )
        if k in clinical_row.index
    ]
    mdata.uns["clinical"] = {
        k: (None if pd.isna(clinical_row[k]) else clinical_row[k]) for k in clin_keys
    }
    mdata.uns["deepspotm"] = {
        "atlas": ATLAS_REPO,
        "n_genes_kept": n_genes,
        "pcc_filtered": bool(gene_keep),
    }
    if "he_encoder" in clinical_row.index and not pd.isna(clinical_row["he_encoder"]):
        mdata.uns["he_encoder"] = clinical_row["he_encoder"]
    # mudata/{sample}.h5mu → he_overview/{sample}.jpg
    overview = out_h5mu.parent.parent / "he_overview" / f"{out_h5mu.stem}.jpg"
    if overview.is_file():
        mdata.uns["he_overview"] = str(overview.resolve())
    out_h5mu.parent.mkdir(parents=True, exist_ok=True)
    mdata.write(out_h5mu)
    _log(f"Wrote {out_h5mu} ({mdata.n_obs} spots, rna={mdata.mod['rna'].n_vars}, HE={mdata.mod['HE'].n_vars})")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out-dir", type=Path, required=True)
    p.add_argument(
        "--stage",
        choices=["all", "select", "download", "fetch", "pcc", "he", "mudata"],
        default="all",
        help="Pipeline stage (default: all)",
    )
    p.add_argument("--max-slides", type=int, default=None, help="Cap number of slides after clinical filter")
    p.add_argument(
        "--max-patients",
        type=int,
        default=None,
        help="Cap unique patients (used with --enrich-pfi; one DX slide kept per case)",
    )
    p.add_argument("--min-spots", type=int, default=500)
    p.add_argument("--require-endpoint", choices=["any", "OS", "PFI", "both"], default="any")
    p.add_argument(
        "--stages",
        default="II,III",
        help=(
            "Comma-separated collapsed AJCC stages to keep (default: II,III). "
            "Sublevels are collapsed (IIA→II, IIIC→III). Use 'all' to disable."
        ),
    )
    p.add_argument(
        "--enrich-pfi",
        default=None,
        metavar="E:C",
        help=(
            "Patient-level PFI event:censored enrichment (e.g. 1:1). "
            "Keeps one DX slide per case; pair with --max-patients."
        ),
    )
    p.add_argument(
        "--enrich-seed",
        type=int,
        default=0,
        help="RNG seed for --enrich-pfi sampling (default: 0)",
    )
    p.add_argument("--skip-download", action="store_true")
    p.add_argument(
        "--pcc-filter",
        action="store_true",
        help=(
            "Subset MuData genes to PCC > --pcc-threshold. Off by default (keep all genes). "
            "Requires --pcc-cameo, --gene-pcc-csv, or --pcc-visium-h5ad."
        ),
    )
    p.add_argument("--pcc-threshold", type=float, default=0.4)
    p.add_argument("--gene-pcc-csv", type=Path, default=None, help="Precomputed gene,pcc CSV")
    p.add_argument(
        "--pcc-cameo",
        action="store_true",
        help=(
            "Compute PCC on theislab/CAMEO-Breast niches (held-out Xenium; not in "
            "DeepSpot-M Visium training). Implies gene subsetting when used with mudata."
        ),
    )
    p.add_argument(
        "--pcc-cameo-genes-csv",
        type=Path,
        default=None,
        help="Optional local Xenium Human Breast v1 metadata CSV (Genes column order)",
    )
    p.add_argument("--pcc-visium-h5ad", type=Path, default=None, help="Measured BRCA Visium h5ad for PCC")
    p.add_argument("--pcc-tile-dir", type=Path, default=None, help="JPEG/PNG tiles aligned to Visium obs_names")
    p.add_argument("--pcc-device", default="cpu")
    p.add_argument(
        "--pcc-max-spots",
        type=int,
        default=2000,
        help="Max Visium spots or CAMEO niches for PCC (reservoir sample)",
    )
    p.add_argument("--download-svs", action="store_true", help="Download GDC SVS for HE extraction")
    p.add_argument(
        "--download-atlas",
        action="store_true",
        help="During HE, fetch missing atlas h5ad for the sample only (no bulk prefetch)",
    )
    p.add_argument(
        "--delete-svs",
        action="store_true",
        help="Delete each SVS after HE .npy is written successfully (saves disk)",
    )
    p.add_argument(
        "--sample-id",
        action="append",
        default=None,
        help="Restrict to this sample_id (repeatable); for parallel HE workers",
    )
    p.add_argument(
        "--reselect",
        action="store_true",
        help="Force re-running clinical select even if selected_slides.csv exists",
    )
    p.add_argument(
        "--he-encoder",
        default="univ2",
        choices=["univ2", "uni2", "deepspotm_midnight", "mean_rgb"],
        help="HE patch encoder for MuData modality HE (default: univ2)",
    )
    p.add_argument(
        "--univ2-model-dir",
        type=Path,
        default=None,
        help="Directory with UNI2 pytorch_model.bin (default: UNIV2_MODEL_DIR or models/)",
    )
    p.add_argument("--he-device", default="auto", choices=["auto", "cuda", "cpu"])
    p.add_argument("--he-batch-size", type=int, default=32, help="UNI2 batch size (default: 32)")
    p.add_argument("--he-npy-dir", type=Path, default=None, help="Directory of {sample_id}.npy HE embeddings")
    p.add_argument("--svs-dir", type=Path, default=None)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    out = args.out_dir
    out.mkdir(parents=True, exist_ok=True)
    stage = args.stage
    sel_path = out / "selected_slides.csv"

    # Downstream stages reuse selected_slides.csv so download/PCC/HE do not
    # overwrite a prior enrich / max-patients selection (unless --reselect).
    reuse_selected = (
        stage in ("download", "fetch", "pcc", "he", "mudata")
        and sel_path.exists()
        and not args.reselect
    )

    if not reuse_selected:
        cdr_path = _download(CDR_URL, out / "clinical" / CDR_FILENAME)
        cdr = load_tcga_cdr(cdr_path)
        cdr.to_csv(out / "clinical" / "TCGA_BRCA_CDR_endpoints.csv", index=False)
        _log(f"CDR BRCA cases with OS and/or PFI: {len(cdr)}")

    if stage in ("all", "select") or (
        stage in ("download", "fetch", "pcc", "he", "mudata") and not reuse_selected
    ):
        meta = download_atlas_metadata(out)
        selected = select_endpoint_slides(meta, cdr)
        if args.require_endpoint == "OS":
            selected = selected[selected["has_OS"]]
        elif args.require_endpoint == "PFI":
            selected = selected[selected["has_PFI"]]
        elif args.require_endpoint == "both":
            selected = selected[selected["has_OS"] & selected["has_PFI"]]
        stages = parse_stages_arg(args.stages)
        selected = filter_by_ajcc_stage(selected, stages)
        if "n_spots" in selected.columns:
            selected = selected[selected["n_spots"] >= args.min_spots]
        selected = selected.sort_values(["case_id", "sample_id"]).reset_index(drop=True)
        enrich_note = ""
        if args.enrich_pfi is not None:
            # Enrichment needs defined PFI; override weaker --require-endpoint.
            selected = selected[selected["has_PFI"]].copy()
            ratio = parse_enrich_ratio(args.enrich_pfi)
            selected = enrich_pfi_balanced(
                selected,
                ratio=ratio,
                max_patients=args.max_patients,
                seed=args.enrich_seed,
            )
            enrich_note = f" enrich-pfi={args.enrich_pfi} max_patients={args.max_patients}"
        elif args.max_patients is not None:
            selected = one_slide_per_patient(selected).head(args.max_patients)
            enrich_note = f" max_patients={args.max_patients}"
        if args.max_slides is not None:
            selected = selected.head(args.max_slides)
        selected.to_csv(sel_path, index=False)
        stage_note = f" stages={stages}" if stages else " stages=all"
        _log(
            f"Selected {len(selected)} DX/primary slides with clinical endpoints"
            f"{stage_note}{enrich_note} -> {sel_path}"
        )
        if stage == "select":
            return
    else:
        selected = pd.read_csv(sel_path)
        _log(f"Reusing {sel_path} ({len(selected)} slides)")

    if args.sample_id:
        keep_ids = set(args.sample_id)
        selected = selected[selected["sample_id"].astype(str).isin(keep_ids)].copy()
        if selected.empty:
            raise SystemExit(f"No rows in selected_slides.csv for --sample-id {sorted(keep_ids)}")
        _log(f"Restricted to {len(selected)} slide(s) via --sample-id")

    if stage in ("all", "download") and not args.skip_download:
        if "file_path" not in selected.columns:
            raise SystemExit("selected_slides.csv needs file_path from atlas metadata")
        download_slides(selected["file_path"].astype(str).tolist(), out)

    # Per-sample fetch for lazyslide → univ2 pipeline (no embedding here).
    if stage == "fetch":
        svs_dir = args.svs_dir or (out / "svs")
        atlas_root = out / "atlas"
        svs_dir.mkdir(parents=True, exist_ok=True)
        for _, row in selected.iterrows():
            sample_id = str(row["sample_id"])
            rel = str(row.get("file_path", "")) if "file_path" in row.index else ""
            adata_path = (atlas_root / rel) if rel and rel != "nan" else None
            if adata_path is None or not adata_path.exists():
                if rel and rel != "nan":
                    _log(f"Fetching atlas h5ad for {sample_id}…")
                    download_slides([rel], out)
                    adata_path = atlas_root / rel
                if adata_path is None or not adata_path.exists():
                    raise SystemExit(f"Missing atlas h5ad for {sample_id}")
            fid = row.get("gdc_file_id")
            if pd.isna(fid) or not fid:
                raise SystemExit(f"No gdc_file_id for {sample_id}")
            svs_path = svs_dir / f"{sample_id}.svs"
            download_gdc_svs(str(fid), svs_path)
            _log(f"Fetched atlas={adata_path} svs={svs_path}")
        return

    # HE does not need gene filtering. MuData keeps all genes unless PCC is requested.
    want_pcc = bool(
        stage == "pcc"
        or args.pcc_filter
        or args.pcc_cameo
        or args.pcc_visium_h5ad
        or args.gene_pcc_csv
    )
    gene_keep: list[str] | None = None
    pcc_csv = args.gene_pcc_csv or (out / "gene_pcc_deepspotm.csv")
    pcc_source = "none"

    if stage == "he":
        gene_keep = None
    elif want_pcc:
        if args.gene_pcc_csv and Path(args.gene_pcc_csv).exists():
            pcc_df = pd.read_csv(args.gene_pcc_csv)
            pcc_source = "gene_pcc_csv"
        elif args.pcc_cameo or (stage == "pcc" and not args.pcc_visium_h5ad and not pcc_csv.exists()):
            pcc_df = compute_gene_pcc_cameo_breast(
                pcc_threshold=args.pcc_threshold,
                out_csv=pcc_csv,
                device=args.pcc_device,
                max_niches=args.pcc_max_spots,
                gene_panel_csv=args.pcc_cameo_genes_csv,
                hf_cache_dir=out / "hf_cache",
            )
            pcc_source = "cameo_breast"
        elif args.pcc_visium_h5ad:
            pcc_df = compute_gene_pcc_deepspotm(
                args.pcc_visium_h5ad,
                args.pcc_tile_dir,
                pcc_threshold=args.pcc_threshold,
                out_csv=pcc_csv,
                device=args.pcc_device,
                max_spots=args.pcc_max_spots,
            )
            pcc_source = "visium"
        elif pcc_csv.exists():
            pcc_df = pd.read_csv(pcc_csv)
            pcc_source = "gene_pcc_csv_or_cache"
        else:
            raise SystemExit(
                "PCC requested but no source found. Pass --pcc-cameo, --gene-pcc-csv, "
                "or --pcc-visium-h5ad — or omit them to keep all genes."
            )
        gene_keep = pcc_df.loc[pcc_df["pcc"] > args.pcc_threshold, "gene"].astype(str).tolist()
        if not gene_keep:
            raise SystemExit(f"No genes with pcc > {args.pcc_threshold}")
        (out / "genes_pcc_gt_threshold.txt").write_text("\n".join(gene_keep) + "\n")
        _log(f"PCC gene filter: keeping {len(gene_keep)} genes (pcc > {args.pcc_threshold})")
        if stage == "pcc":
            return
    else:
        _log("No PCC gene filter: MuData will keep all atlas genes")

    he_dir = args.he_npy_dir or (out / "he_embeddings")
    svs_dir = args.svs_dir or (out / "svs")
    atlas_root = out / "atlas"

    if stage in ("all", "he"):
        he_dir.mkdir(parents=True, exist_ok=True)
        for _, row in selected.iterrows():
            sample_id = str(row["sample_id"])
            npy = he_dir / f"{sample_id}.npy"
            svs_path = svs_dir / f"{sample_id}.svs"
            overview = he_overview_path(out, sample_id)
            if npy.exists():
                if not overview.exists() and svs_path.exists():
                    save_he_overview(svs_path, overview)
                if args.delete_svs and svs_path.exists():
                    svs_path.unlink()
                    _log(f"Deleted leftover SVS (npy exists): {svs_path}")
                continue
            # Locate prediction h5ad; optionally fetch just this slide (disk-friendly).
            rel = str(row.get("file_path", "")) if "file_path" in row.index else ""
            adata_path = (atlas_root / rel) if rel and rel != "nan" else None
            if adata_path is None or not adata_path.exists():
                if rel and rel != "nan" and (args.download_atlas or stage == "all"):
                    _log(f"Fetching atlas h5ad for {sample_id}…")
                    download_slides([rel], out)
                    adata_path = atlas_root / rel
                if adata_path is None or not adata_path.exists():
                    candidates = list(atlas_root.glob(f"**/TCGA_BRCA/**/{sample_id}.h5ad"))
                    if not candidates:
                        candidates = list(atlas_root.glob(f"**/{sample_id}.h5ad"))
                    if not candidates:
                        _log(f"Skip HE (missing h5ad): {sample_id}")
                        continue
                    adata_path = candidates[0]
            fid = row.get("gdc_file_id")
            if pd.isna(fid) or not fid:
                _log(f"Skip HE (no gdc_file_id): {sample_id}")
                continue
            if not args.download_svs and not svs_path.exists():
                _log(f"Skip HE (pass --download-svs or place SVS in {svs_dir}): {sample_id}")
                continue
            # One SVS at a time: download (validated) → embed → delete on success only.
            try:
                svs = download_gdc_svs(str(fid), svs_path)
                if not overview.exists():
                    save_he_overview(svs, overview)
                extract_he_embeddings_from_svs(
                    adata_path,
                    svs,
                    npy,
                    encoder=args.he_encoder,
                    univ2_model_dir=args.univ2_model_dir,
                    device=args.he_device,
                    batch_size=args.he_batch_size,
                )
            except Exception:
                if svs_path.exists() and not _is_tiff_like(svs_path):
                    svs_path.unlink()
                    _log(f"Removed bad/partial SVS after error: {svs_path}")
                raise
            if args.delete_svs and svs_path.exists():
                svs_path.unlink()
                _log(f"Deleted SVS after HE: {svs_path}")
        if stage == "he":
            return

    if stage in ("all", "mudata"):
        mudata_dir = out / "mudata"
        mudata_dir.mkdir(parents=True, exist_ok=True)
        manifest = []
        for _, row in selected.iterrows():
            sample_id = str(row["sample_id"])
            rel = str(row.get("file_path", ""))
            adata_path = atlas_root / rel if rel else None
            if adata_path is None or not Path(adata_path).exists():
                candidates = list(atlas_root.glob(f"**/{sample_id}.h5ad"))
                if not candidates:
                    _log(f"Skip mudata (missing h5ad): {sample_id}")
                    continue
                adata_path = candidates[0]
            he_npy = he_dir / f"{sample_id}.npy"
            if not he_npy.exists():
                _log(f"Skip mudata (missing HE npy): {sample_id}")
                continue
            out_h5mu = mudata_dir / f"{sample_id}.h5mu"
            row = row.copy()
            row["he_encoder"] = args.he_encoder
            assemble_mudata(Path(adata_path), he_npy, gene_keep, row, out_h5mu)
            manifest.append({"sample_id": sample_id, "h5mu": str(out_h5mu), "case_id": row.get("case_id")})
        pd.DataFrame(manifest).to_csv(out / "mudata_manifest.csv", index=False)
        summary = {
            "n_slides": len(manifest),
            "n_genes": None if gene_keep is None else len(gene_keep),
            "pcc_filter": gene_keep is not None,
            "pcc_threshold": args.pcc_threshold if gene_keep is not None else None,
            "pcc_source": pcc_source,
            "endpoint": "OS and/or PFI (TCGA-CDR); DX/FFPE primary tumor at diagnosis",
            "ajcc_stages": parse_stages_arg(args.stages),
            "enrich_pfi": args.enrich_pfi,
            "max_patients": args.max_patients,
            "enrich_seed": args.enrich_seed if args.enrich_pfi else None,
            "he_encoder": args.he_encoder,
            "mudata_modalities": ["rna", "HE"],
        }
        (out / "summary.json").write_text(json.dumps(summary, indent=2))
        _log(json.dumps(summary))


if __name__ == "__main__":
    main()
