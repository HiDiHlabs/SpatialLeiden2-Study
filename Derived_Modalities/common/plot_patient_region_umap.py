"""Plot any spot embedding coloured by patient and by region, in the style of Kuppe et al. 2022,
Extended Data Fig. 3a-b (and 5a).

Two panels side by side, "patients" and "regions", with the authors' palettes and look:
  * patient colours: ArchR::paletteDiscrete default ("stallion"), which the authors use for all
    discrete UMAPs (saezlab/visium_heart, plot_figures/*.ipynb); P1 ... P20 in the paper's order
  * region colours: `cols.region` from the authors' plot_figures/fig4.ipynb
    (control #4CAF50, FZ #009DB2, IZ #B50084, RZ #FFC107, BZ #F57F17)
  * no frame or ticks; two short arrows labelled UMAP 1 / UMAP 2 meeting at the bottom left;
    legends outside the data; points drawn in random order so no group hides another

Input, one of:
  A) a table (.csv / .csv.gz / .tsv) with one row per spot: two coordinate columns plus the
     patient and region columns
         python plot_patient_region_umap.py --table ccc_umap.csv.gz --out ccc_patient_region
  B) an .h5ad: coordinates from an obsm key (default X_umap), labels from obs
         python plot_patient_region_umap.py --h5ad ccc_joint.h5ad --basis X_umap --out ccc_patient_region
  C) an .h5ad without a UMAP yet: add --compute-umap. The UMAP is built from --use-rep (an obsm
     key, e.g. X_pca or X_harmony; default: X), with 15 nearest neighbours and seed 0. If
     --use-rep is X and has more than 50 columns, X is z-scored per feature and reduced to 30 PCs
     first. Needs scanpy.
         python plot_patient_region_umap.py --h5ad ccc_joint.h5ad --compute-umap --use-rep X_pca --out ccc_patient_region

Region labels are matched case-insensitively; CTRL / ctrl / Control are shown as "control".
Labels other than the five regions or P1-P20 still get a colour (from the end of the palette),
and are listed in the console.

Options: --x-col/--y-col (default UMAP1/UMAP2), --patient-col (default patient),
--region-col (default region; for the Kuppe spot tables use major_labl), --point-size (default 1),
--title-prefix (e.g. "CCC: "), --seed (default 0). Output: <out>.png (300 dpi) and <out>.pdf.

Dependencies: numpy, pandas, matplotlib; anndata for --h5ad; scanpy for --compute-umap.
Written for SpaceHack 2026 track C (SpatialLeiden2-Study, Derived_Modalities).
"""
import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ArchR "stallion" palette, in the order the paper's legend assigns it to P1 ... P20
STALLION = ["#D51F26", "#272E6A", "#208A42", "#89288F", "#F47D2B", "#FEE500", "#8A9FD1", "#C06CAB",
            "#D8A767", "#90D5E4", "#89C75F", "#F37B7D", "#9983BD", "#D24B27", "#3BBCA8", "#6E4B9E",
            "#0C727C", "#7E1416", "#E6C2DC", "#3D3D3D"]
PATIENT_COLOURS = {f"P{i + 1}": c for i, c in enumerate(STALLION)}
REGION_COLOURS = {"control": "#4CAF50", "FZ": "#009DB2", "IZ": "#B50084", "RZ": "#FFC107", "BZ": "#F57F17"}
EXTRA = ["#7f7f7f", "#17becf", "#bcbd22", "#8c564b", "#e377c2", "#1f77b4", "#ff7f0e", "#2ca02c"]


def normalise_region(r: str) -> str:
    key = str(r).strip()
    lut = {k.lower(): k for k in REGION_COLOURS}
    lut.update({"ctrl": "control"})
    return lut.get(key.lower(), key)


def palette(labels, known: dict) -> dict:
    """Known labels get the paper's colour; anything else gets a fallback colour."""
    cols, extra = {}, iter(EXTRA * 10)
    unknown = [lab for lab in labels if lab not in known]
    if unknown:
        print(f"[note] labels without a paper colour (fallback colours used): {unknown}", file=sys.stderr)
    for lab in labels:
        cols[lab] = known.get(lab) or next(extra)
    return cols


def patient_order(labels):
    def key(p):
        s = str(p)
        return (0, int(s[1:])) if s[:1] == "P" and s[1:].isdigit() else (1, s)
    return sorted(labels, key=key)


def umap_arrows(ax, title: str) -> None:
    """Paper style: no frame or ticks; two arrows meeting at the bottom-left corner."""
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set(xticks=[], yticks=[], xlabel="", ylabel="")
    ax.set_title(title, fontsize=12)
    arrow = dict(arrowstyle="-|>", lw=1.3, color="black", mutation_scale=12, shrinkA=0, shrinkB=0,
                 joinstyle="miter", capstyle="projecting")
    ax.annotate("", xy=(0.15, 0.0), xytext=(0.0, 0.0), xycoords="axes fraction", arrowprops=arrow)
    ax.annotate("", xy=(0.0, 0.15), xytext=(0.0, 0.0), xycoords="axes fraction", arrowprops=arrow)
    ax.text(0.075, -0.035, "UMAP 1", transform=ax.transAxes, ha="center", va="top", fontsize=10)
    ax.text(-0.035, 0.075, "UMAP 2", transform=ax.transAxes, ha="right", va="center", rotation=90, fontsize=10)


def legend(ax, cols: dict, order, ncol: int = 1) -> None:
    handles = [plt.Line2D([], [], ls="", marker="o", markersize=6, color=cols[k], label=k) for k in order]
    ax.legend(handles=handles, fontsize=9, ncol=ncol, handletextpad=0.3, labelspacing=0.35,
              loc="center left", bbox_to_anchor=(0.98, 0.5), frameon=False)


def plot(xy: np.ndarray, patient, region, out: Path, point_size=1.0, title_prefix="", seed=0) -> None:
    patient = pd.Series(patient, dtype=str).to_numpy()
    region = pd.Series(region, dtype=str).map(normalise_region).to_numpy()
    pats = patient_order(pd.unique(patient))
    regs = [r for r in REGION_COLOURS if r in set(region)] + sorted(set(region) - set(REGION_COLOURS))
    pcol, rcol = palette(pats, PATIENT_COLOURS), palette(regs, REGION_COLOURS)

    order = np.random.default_rng(seed).permutation(len(xy))
    fig, ax = plt.subplots(1, 2, figsize=(11.5, 5.2))
    for a, lab, cols, keys, title, panel in [(ax[0], patient, pcol, pats, "patients", "a"),
                                             (ax[1], region, rcol, regs, "regions", "b")]:
        a.scatter(xy[order, 0], xy[order, 1], c=[cols[k] for k in lab[order]], s=point_size,
                  linewidths=0, rasterized=True)
        umap_arrows(a, f"{title_prefix}{title}")
        legend(a, cols, keys, ncol=2 if len(keys) > 22 else 1)
        a.text(-0.08, 1.04, panel, transform=a.transAxes, fontsize=15, fontweight="bold")
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    print(f"wrote {out.with_suffix('.png')} and {out.with_suffix('.pdf')} ({len(xy)} spots, "
          f"{len(pats)} patients, regions {regs})")


def compute_umap(a, use_rep: str, seed: int) -> np.ndarray:
    import scanpy as sc

    if use_rep == "X":
        X = a.X.toarray() if hasattr(a.X, "toarray") else np.asarray(a.X)
        if X.shape[1] > 50:
            X = (X - X.mean(0)) / np.where(X.std(0) > 0, X.std(0), 1)
            a.obsm["X_pca_for_umap"] = sc.pp.pca(X, n_comps=min(30, X.shape[1] - 1), random_state=seed)
            use_rep = "X_pca_for_umap"
        else:
            a.obsm["X_for_umap"] = X
            use_rep = "X_for_umap"
    sc.pp.neighbors(a, n_neighbors=15, use_rep=use_rep, random_state=seed)
    sc.tl.umap(a, random_state=seed)
    return np.asarray(a.obsm["X_umap"])


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--table", type=Path, help="csv/tsv with coordinates and labels, one row per spot")
    src.add_argument("--h5ad", type=Path, help="AnnData with labels in obs")
    p.add_argument("--out", type=Path, required=True, help="output path without extension")
    p.add_argument("--x-col", default="UMAP1")
    p.add_argument("--y-col", default="UMAP2")
    p.add_argument("--basis", default="X_umap", help="obsm key with the 2D coordinates (--h5ad)")
    p.add_argument("--compute-umap", action="store_true", help="build the UMAP first (--h5ad; needs scanpy)")
    p.add_argument("--use-rep", default="X", help="obsm key (or X) the UMAP is built from")
    p.add_argument("--patient-col", default="patient")
    p.add_argument("--region-col", default="region")
    p.add_argument("--point-size", type=float, default=1.0)
    p.add_argument("--title-prefix", default="")
    p.add_argument("--seed", type=int, default=0)
    args = p.parse_args()

    if args.table:
        sep = "\t" if args.table.name.endswith((".tsv", ".tsv.gz")) else ","
        df = pd.read_csv(args.table, sep=sep)
        missing = [c for c in (args.x_col, args.y_col, args.patient_col, args.region_col) if c not in df]
        if missing:
            p.error(f"columns not in {args.table}: {missing}; available: {list(df.columns)}")
        xy, obs = df[[args.x_col, args.y_col]].to_numpy(float), df
    else:
        import anndata as ad

        a = ad.read_h5ad(args.h5ad)
        missing = [c for c in (args.patient_col, args.region_col) if c not in a.obs]
        if missing:
            p.error(f"obs columns not in {args.h5ad}: {missing}; available: {list(a.obs.columns)}")
        if args.compute_umap:
            xy = compute_umap(a, args.use_rep, args.seed)
        elif args.basis in a.obsm:
            xy = np.asarray(a.obsm[args.basis])[:, :2]
        else:
            p.error(f"obsm['{args.basis}'] not found (have {list(a.obsm)}); pass --basis or --compute-umap")
        obs = a.obs
    plot(xy, obs[args.patient_col].to_numpy(), obs[args.region_col].to_numpy(), args.out,
         args.point_size, args.title_prefix, args.seed)


if __name__ == "__main__":
    main()
