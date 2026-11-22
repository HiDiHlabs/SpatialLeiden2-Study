#!/usr/bin/env python3
"""SpatialLeiden2-style 3D figures for TRICARE Multileiden outputs.

Mirrors HiDiHlabs/SpatialLeiden2-Study OpenST_3D_lymphnode + gif/gif.py, with
TRICARE extras on the slice panel:
  - square markers (tile footprints)
  - H&E plane thumbnail next to each depth
  - pathologist risk score strip below

Input: ``{TRICARE}/univ2_multileiden/multileiden/{cohort}/{specimen}.h5mu``
H&E planes: ``{TRICARE}/images/**/{slide_id}.tif`` (or source mudata tile_path).

Usage
-----
  export TRICARE=./data/TRICARE
  conda env create -f TRICARE_3D/tricare_plot3d.yml
  conda activate tricare_plot3d

  python TRICARE_3D/plot_tricare_3d_spatialleiden.py \\
    --tricare-dir \"$TRICARE\" --cohort prostate_test --specimen PRT_P000
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import mudata as md
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.gridspec import GridSpec
from PIL import Image

try:
    import colorcet as cc
except ImportError as e:
    raise SystemExit("Install colorcet (conda env: TRICARE_3D/tricare_plot3d.yml)") from e


def _log(msg: str) -> None:
    print(msg, flush=True)


def _cluster_col(obs: pd.DataFrame) -> str:
    for c in ("spatialleiden", "multileiden"):
        if c in obs.columns:
            return c
    raise SystemExit(f"No spatialleiden/multileiden in obs columns: {list(obs.columns)}")


def annotation_to_risk(label: object, is_annotated: object = True) -> float:
    """Map TRICARE pathologist labels → risk in [0, 1] (nan if unlabeled)."""
    del is_annotated  # label presence is authoritative; neighbors are usually blank
    if label is None or (isinstance(label, float) and np.isnan(label)):
        return float("nan")
    s = str(label).strip().lower()
    if s in ("", "nan", "none", "<na>"):
        return float("nan")
    if any(k in s for k in ("tumor", "higher", "dysplas", "neoplas", "cancer", "adenocarcinoma")):
        return 1.0
    if any(k in s for k in ("normal", "low", "benign")):
        return 0.0
    return float("nan")


def labels_from_h5mu(h5mu: Path) -> pd.DataFrame:
    """Build SpatialLeiden2-style labels table (parquet-ready)."""
    m = md.read(h5mu)
    if "spatial" not in m.obsm:
        raise SystemExit(f"{h5mu}: missing obsm['spatial']")
    xyz = np.asarray(m.obsm["spatial"], dtype=np.float64)
    if xyz.shape[1] != 3:
        raise SystemExit(f"{h5mu}: expected spatial (n,3), got {xyz.shape}")

    col = _cluster_col(m.obs)
    df = pd.DataFrame(
        {
            "spatialleiden": m.obs[col].astype(str).values,
            "x": xyz[:, 0],
            "y": xyz[:, 1],
            "z": xyz[:, 2],
        },
        index=m.obs_names.astype(str),
    )
    if "depth" in m.obs.columns:
        df["depth"] = m.obs["depth"].astype(int).values
    else:
        z_levels = {z: i for i, z in enumerate(sorted(df["z"].unique()))}
        df["depth"] = df["z"].map(z_levels).astype(int)
    df["n_section"] = df["depth"]
    if "pathologist_annotation" in m.obs.columns:
        df["pathologist_annotation"] = m.obs["pathologist_annotation"].astype(str).values
    if "is_annotated" in m.obs.columns:
        df["is_annotated"] = m.obs["is_annotated"].astype(bool).values
    else:
        df["is_annotated"] = False
    if "plane_slide_id" in m.obs.columns:
        df["plane_slide_id"] = m.obs["plane_slide_id"].astype(str).values
    elif "slide_id" in m.obs.columns:
        df["plane_slide_id"] = m.obs["slide_id"].astype(str).values
    df["risk_score"] = [
        annotation_to_risk(a, ia)
        for a, ia in zip(df.get("pathologist_annotation", [""] * len(df)), df["is_annotated"])
    ]
    df["spatialleiden"] = pd.Categorical(df["spatialleiden"])
    return df


def _palette(labels: pd.DataFrame) -> dict:
    cats = list(labels["spatialleiden"].cat.categories)
    colors = list(cc.glasbey) if hasattr(cc, "glasbey") else sns.color_palette("husl", n_colors=max(len(cats), 1))
    if isinstance(colors[0], str):
        return {c: colors[i % len(colors)] for i, c in enumerate(cats)}
    return dict(zip(cats, colors[: len(cats)]))


def _index_plane_tiffs(images_root: Path) -> dict[str, Path]:
    out: dict[str, Path] = {}
    if not images_root.is_dir():
        return out
    for p in images_root.rglob("*"):
        if p.suffix.lower() in {".tif", ".tiff"} and p.is_file():
            out[p.stem] = p.resolve()
    return out


def _plane_path_for_depth(
    depth: int,
    labels: pd.DataFrame,
    plane_index: dict[str, Path],
    specimen: str,
) -> Path | None:
    sub = labels.loc[labels["n_section"] == depth]
    if "plane_slide_id" in sub.columns:
        for sid in sub["plane_slide_id"].dropna().unique():
            p = plane_index.get(str(sid))
            if p is not None and p.is_file():
                return p
    # naming: {specimen}_{depth}
    for key in (f"{specimen}_{int(depth)}", f"{specimen}_{int(depth):04d}"):
        p = plane_index.get(key)
        if p is not None and p.is_file():
            return p
    # fuzzy: stem endswith _{depth}
    suffix = f"_{int(depth)}"
    for stem, path in plane_index.items():
        if stem.endswith(suffix) and specimen in stem:
            return path
    return None


def _load_he_thumbnail(path: Path, max_side: int = 512) -> np.ndarray | None:
    try:
        import tifffile

        arr = tifffile.imread(path)
    except Exception:
        try:
            arr = np.asarray(Image.open(path))
        except Exception as e:
            _log(f"  warn: could not read H&E {path.name}: {e}")
            return None
    arr = np.asarray(arr)
    if arr.ndim == 2:
        arr = np.stack([arr, arr, arr], axis=-1)
    if arr.ndim == 3 and arr.shape[-1] > 3:
        arr = arr[..., :3]
    if arr.dtype != np.uint8:
        a = arr.astype(np.float64)
        if a.max() > 0:
            a = a / a.max() * 255.0
        arr = np.clip(a, 0, 255).astype(np.uint8)
    img = Image.fromarray(arr)
    img.thumbnail((max_side, max_side), Image.Resampling.LANCZOS)
    return np.asarray(img)


def _depth_risk_table(labels: pd.DataFrame) -> pd.DataFrame:
    """One row per depth: risk + pathologist label."""
    rows = []
    for depth, g in labels.groupby("n_section", sort=True):
        ann = g.loc[g["is_annotated"]] if "is_annotated" in g.columns else g
        if len(ann) == 0:
            ann = g
        label = ""
        if "pathologist_annotation" in ann.columns:
            vals = [str(v) for v in ann["pathologist_annotation"] if str(v) not in ("", "nan", "None")]
            label = vals[0] if vals else ""
        risk = float(np.nanmean(g["risk_score"].to_numpy(dtype=float))) if "risk_score" in g.columns else float("nan")
        if np.isnan(risk) and label:
            risk = annotation_to_risk(label, True)
        rows.append({"depth": int(depth), "risk_score": risk, "pathologist_annotation": label})
    return pd.DataFrame(rows)


def _merge_optional_risk_csv(depth_df: pd.DataFrame, risk_csv: Path | None) -> pd.DataFrame:
    if risk_csv is None:
        return depth_df
    extra = pd.read_csv(risk_csv)
    cols = {c.lower(): c for c in extra.columns}
    dcol = cols.get("depth") or cols.get("n_section")
    rcol = cols.get("risk_score") or cols.get("risk") or cols.get("score")
    if dcol is None or rcol is None:
        raise SystemExit(f"--risk-csv needs depth + risk_score columns; got {list(extra.columns)}")
    m = extra[[dcol, rcol]].rename(columns={dcol: "depth", rcol: "risk_score"})
    m["depth"] = m["depth"].astype(int)
    out = depth_df.drop(columns=["risk_score"], errors="ignore").merge(m, on="depth", how="left")
    return out


def _select_sections(sections: list, max_sections: int | None) -> list:
    if max_sections is None or len(sections) <= max_sections:
        return list(sections)
    pick = np.linspace(0, len(sections) - 1, max_sections).round().astype(int)
    return [sections[i] for i in pick]


def plot_facet_slices(
    labels: pd.DataFrame,
    out_pdf: Path,
    *,
    specimen: str,
    plane_index: dict[str, Path],
    col_wrap: int = 4,
    point_size: float = 28.0,
    max_sections: int | None = None,
    he_max_side: int = 512,
    risk_csv: Path | None = None,
) -> None:
    """Per-depth SpatialLeiden squares | H&E, with risk score strip below."""
    all_sections = sorted(int(x) for x in labels["n_section"].unique())
    sections = _select_sections(all_sections, max_sections)
    df = labels[labels["n_section"].isin(sections)].copy()
    palette = _palette(labels)
    depth_risk = _merge_optional_risk_csv(_depth_risk_table(labels), risk_csv)

    n = len(sections)
    ncols = min(col_wrap, max(1, n))
    nrows = int(np.ceil(n / ncols))
    # each depth uses 2 subplot columns: clusters | H&E
    fig_w = 3.2 * ncols * 2
    fig_h = 3.0 * nrows + 1.4
    fig = plt.figure(figsize=(fig_w, fig_h))
    gs = GridSpec(
        nrows + 1,
        ncols * 2,
        figure=fig,
        height_ratios=[1.0] * nrows + [0.45],
        hspace=0.35,
        wspace=0.08,
    )

    xlim = (df["x"].min(), df["x"].max())
    ylim = (df["y"].min(), df["y"].max())
    # flip y so image-like orientation matches H&E
    ylim_plot = (ylim[1], ylim[0]) if ylim[0] < ylim[1] else ylim

    for i, depth in enumerate(sections):
        r, c = divmod(i, ncols)
        ax_c = fig.add_subplot(gs[r, 2 * c])
        ax_h = fig.add_subplot(gs[r, 2 * c + 1])
        sdf = df[df["n_section"] == depth]

        for cl, g in sdf.groupby("spatialleiden", observed=True):
            ax_c.scatter(
                g["x"],
                g["y"],
                c=[palette[str(cl)]],
                s=point_size,
                marker="s",
                linewidths=0,
                rasterized=True,
            )
        ax_c.set(aspect="equal", xlim=xlim, ylim=ylim_plot, xticks=[], yticks=[])
        ax_c.set_title(f"depth {depth}", fontsize=9)
        ax_c.set_xlabel("SpatialLeiden", fontsize=8)

        plane = _plane_path_for_depth(depth, labels, plane_index, specimen)
        he = _load_he_thumbnail(plane, max_side=he_max_side) if plane is not None else None
        if he is not None:
            ax_h.imshow(he)
        else:
            ax_h.text(0.5, 0.5, "H&E missing", ha="center", va="center", transform=ax_h.transAxes)
        ax_h.set(xticks=[], yticks=[])
        ax_h.set_xlabel("H&E", fontsize=8)
        lab = ""
        if "pathologist_annotation" in sdf.columns:
            vals = [str(v) for v in sdf["pathologist_annotation"] if str(v) not in ("", "nan", "None")]
            lab = vals[0] if vals else ""
        risk_row = depth_risk.loc[depth_risk["depth"] == depth]
        risk_v = float(risk_row["risk_score"].iloc[0]) if len(risk_row) else float("nan")
        title = "H&E"
        if lab:
            title = lab
        if not np.isnan(risk_v):
            title = f"{title}  (risk={risk_v:.2f})" if lab else f"risk={risk_v:.2f}"
        ax_h.set_title(title, fontsize=8)

    # blank unused cells
    for i in range(n, nrows * ncols):
        r, c = divmod(i, ncols)
        for j in (0, 1):
            ax = fig.add_subplot(gs[r, 2 * c + j])
            ax.axis("off")

    # risk strip across full width
    ax_r = fig.add_subplot(gs[nrows, :])
    d = depth_risk.sort_values("depth")
    ax_r.plot(d["depth"], d["risk_score"], color="0.35", lw=1.2, marker="o", ms=3, label="pathologist risk")
    for depth in sections:
        ax_r.axvline(depth, color="C0", alpha=0.35, lw=1)
    ax_r.set_ylim(-0.05, 1.05)
    ax_r.set_yticks([0, 0.5, 1.0])
    ax_r.set_yticklabels(["low", "0.5", "high"])
    ax_r.set_xlabel("depth")
    ax_r.set_ylabel("risk score")
    ax_r.set_title("Risk score by depth (pathologist; tumor=1, normal=0)", fontsize=9)
    ax_r.legend(loc="upper right", frameon=False, fontsize=8)

    fig.suptitle(f"{specimen}  SpatialLeiden + H&E", fontsize=12, y=0.995)
    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_pdf, bbox_inches="tight", dpi=200)
    plt.close(fig)
    _log(f"Wrote {out_pdf}")


def plot_gif(
    labels: pd.DataFrame,
    out_gif: Path,
    *,
    fps: int = 2,
    width: float = 6.0,
    height: float = 6.0,
    point_size: float = 28.0,
) -> None:
    """SpatialLeiden2 gif/gif.py: one frame per depth (square markers)."""
    try:
        from celluloid import Camera
    except ImportError as e:
        raise SystemExit("Install celluloid (conda env: TRICARE_3D/tricare_plot3d.yml)") from e

    df = labels.copy()
    df["spatialleiden"] = df["spatialleiden"].astype("category")
    palette = _palette(df)
    xlim = (df["x"].min(), df["x"].max())
    ylim = (df["y"].min(), df["y"].max())
    ylim_plot = (ylim[1], ylim[0]) if ylim[0] < ylim[1] else ylim

    fig, ax = plt.subplots(figsize=(width, height))
    camera = Camera(fig)
    for depth, sdf in df.groupby("n_section", observed=True, sort=True):
        for cl, g in sdf.groupby("spatialleiden", observed=True):
            ax.scatter(
                g["x"],
                g["y"],
                c=[palette[str(cl)]],
                s=point_size,
                marker="s",
                linewidths=0,
            )
        ax.set(aspect="equal", xlim=xlim, ylim=ylim_plot, title=f"depth {depth}")
        ax.set(xticks=[], yticks=[], xlabel=None, ylabel=None)
        ax.tick_params(left=False, bottom=False)
        fig.tight_layout()
        camera.snap()

    out_gif.parent.mkdir(parents=True, exist_ok=True)
    camera.animate().save(out_gif, fps=fps, writer="pillow")
    plt.close(fig)
    _log(f"Wrote {out_gif}")


def plot_3d_scatter(
    labels: pd.DataFrame,
    out_pdf: Path,
    *,
    point_size: float = 8.0,
    elev: float = 20.0,
    azim: float = -60.0,
    max_points: int = 80_000,
) -> None:
    """Single 3D scatter coloured by SpatialLeiden (square markers)."""
    df = labels
    if len(df) > max_points:
        df = df.sample(max_points, random_state=0)
    palette = _palette(labels)
    colors = [palette[c] for c in df["spatialleiden"].astype(str)]

    fig = plt.figure(figsize=(7, 6))
    ax = fig.add_subplot(111, projection="3d")
    ax.scatter(
        df["x"],
        df["y"],
        df["z"],
        c=colors,
        s=point_size,
        marker="s",
        linewidths=0,
        alpha=0.85,
    )
    ax.set_xlabel("x (µm)")
    ax.set_ylabel("y (µm)")
    ax.set_zlabel("z (µm)")
    ax.view_init(elev=elev, azim=azim)
    ax.set_title("TRICARE SpatialLeiden 3D")
    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_pdf, bbox_inches="tight", dpi=200)
    plt.close(fig)
    _log(f"Wrote {out_pdf}")


def process_specimen(
    h5mu: Path,
    out_dir: Path,
    args: argparse.Namespace,
    plane_index: dict[str, Path],
) -> None:
    specimen = h5mu.stem
    fig_dir = out_dir / specimen
    fig_dir.mkdir(parents=True, exist_ok=True)
    _log(f"[{specimen}] {h5mu}")

    labels = labels_from_h5mu(h5mu)
    pq = fig_dir / f"{specimen}_spatialleiden.parquet"
    labels.to_parquet(pq)
    _log(f"Wrote {pq} (n={len(labels)}, depths={labels['n_section'].nunique()})")

    if not args.skip_facet:
        plot_facet_slices(
            labels,
            fig_dir / f"{specimen}_slices.pdf",
            specimen=specimen,
            plane_index=plane_index,
            col_wrap=args.col_wrap,
            point_size=args.point_size,
            max_sections=args.max_sections,
            he_max_side=args.he_max_side,
            risk_csv=args.risk_csv,
        )
    if not args.skip_gif:
        plot_gif(
            labels,
            fig_dir / f"{specimen}_depth.gif",
            fps=args.fps,
            width=args.width,
            height=args.height,
            point_size=args.point_size,
        )
    if not args.skip_3d:
        plot_3d_scatter(
            labels,
            fig_dir / f"{specimen}_3d.pdf",
            point_size=max(args.point_size / 3.0, 4.0),
            elev=args.elev,
            azim=args.azim,
            max_points=args.max_points,
        )


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--tricare-dir", type=Path, required=True)
    p.add_argument("--out-dir", type=Path, default=None, help="Default: {tricare}/univ2_multileiden/figures")
    p.add_argument("--cohort", default="prostate_test")
    p.add_argument("--specimen", default=None, help="Single specimen stem (default: all *.h5mu)")
    p.add_argument("--max-specimens", type=int, default=None)
    p.add_argument("--col-wrap", type=int, default=4)
    p.add_argument("--max-sections", type=int, default=16, help="Max depth panels in slices PDF")
    p.add_argument("--point-size", type=float, default=28.0, help="Marker area for square tiles")
    p.add_argument("--he-max-side", type=int, default=512, help="H&E thumbnail max side (px)")
    p.add_argument(
        "--risk-csv",
        type=Path,
        default=None,
        help="Optional CSV with columns depth,risk_score (overrides pathologist 0/1)",
    )
    p.add_argument("--images-dir", type=Path, default=None, help="Extra TIFF root (default: {tricare}/images)")
    p.add_argument("--fps", type=int, default=2)
    p.add_argument("--width", type=float, default=6.0)
    p.add_argument("--height", type=float, default=6.0)
    p.add_argument("--elev", type=float, default=20.0)
    p.add_argument("--azim", type=float, default=-60.0)
    p.add_argument("--max-points", type=int, default=80_000)
    p.add_argument("--skip-facet", action="store_true")
    p.add_argument("--skip-gif", action="store_true")
    p.add_argument("--skip-3d", action="store_true")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    root = args.tricare_dir
    ml_dir = root / "univ2_multileiden" / "multileiden" / args.cohort
    out = args.out_dir or (root / "univ2_multileiden" / "figures" / args.cohort)
    if not ml_dir.is_dir():
        raise SystemExit(f"Missing Multileiden dir: {ml_dir}")

    plane_index = _index_plane_tiffs(root / "images")
    if args.images_dir is not None:
        plane_index.update(_index_plane_tiffs(args.images_dir))
    # Prefer RGB cache from tiling if present
    cache = root / "univ2_multileiden" / "wsi_cache"
    if cache.is_dir():
        for p in cache.rglob("*_rgb.tiff"):
            stem = p.name.replace("_rgb.tiff", "").replace("_rgb.tif", "")
            plane_index.setdefault(stem, p.resolve())
    _log(f"Indexed {len(plane_index)} H&E planes")

    if args.specimen:
        paths = [ml_dir / f"{args.specimen}.h5mu"]
        if not paths[0].is_file():
            raise SystemExit(f"Missing {paths[0]}")
    else:
        paths = sorted(ml_dir.glob("*.h5mu"))
    if args.max_specimens is not None:
        paths = paths[: args.max_specimens]
    if not paths:
        raise SystemExit(f"No .h5mu under {ml_dir}")

    _log(f"Plotting {len(paths)} specimens → {out}")
    for h5mu in paths:
        process_specimen(h5mu, out, args, plane_index)
    _log("Done.")


if __name__ == "__main__":
    main()
