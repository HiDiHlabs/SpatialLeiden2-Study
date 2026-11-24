"""Kuppe et al. 2022, Fig. 1h-i reproduced, with panel j replaced by SpatialLeiden domains,
on a control, a border-zone, an ischaemic-zone and a fibrotic-zone slide:
    h  cell2location proportions of vCMs, fibroblasts and myeloid cells (as published)
    i  PROGENy pathway activity, TGFb and NFkB (as published, recomputed)
    j  NEW: SpatialLeiden domains (K = 9) from one joint clustering of all 27 slides on
       Harmony-integrated RNA, at layer_ratio 0 / 0.5 / 1 / 1.5 (seed 0; FIG_SL_RATIOS to change). Replaces the paper's
       TF binding activity panel.
Writes kuppe_fig1h.{pdf,png} (h alone) and kuppe_fig1hij.{pdf,png} (h, i and new j).

Slides (rows of the published panel, mapped with visium_raw/metadata-Visium.csv):
    Ctrl tissue (P7) = control_P7, ACH003
    BZ tissue (P3)   = RZ/BZ_P3,   10X0026
    IZ tissue (P3)   = IZ_P3,      10X0017
    FZ tissue (P14)  = FZ_P14,     ACH005

Inputs:
  h  `prop_<celltype>` from visium_tables/<run>_spots.csv: the authors' own c2l_props values,
     extracted from their Seurat objects by
     01_data_preprocess/kuppe_visium_exploration/01_extract_seurat_tables.R.
     Proportions, not absolute abundance (c2l): abundance does not reproduce the IZ (P3) and
     FZ (P14) vCM patterns; proportions do.
  i  modalities/progeny/progeny_<run>.csv, recomputed with the authors' formula
     (02_modalities/progeny/progeny_kuppe.py). Their own progeny assay is in Seurat objects
     not on this VM.
  j  results/kuppe_2022/exp2_joint_harmony/exp2_labels.csv.gz
     (03_single_modality_clustering/exp2_joint_harmony/exp2_joint_sweep.py). Domains are
     shared across slides (one multi-sample clustering). They are numbered D1-D9 by size at
     ratio 0; at the other ratios each domain takes the number of the ratio-0 domain it
     overlaps most (Hungarian matching over all 88,704 spots), so a colour means the same
     domain in every row and column.
Colour scale (h, i): per panel, clipped at the 1st and 99th percentile (Seurat min.cutoff
"q1", max.cutoff "q99", as in st_snRNAseq/07_structure/plot_spatial_ex_panels.R). Palettes,
black background, dashed outline and arrows follow the published panel (styled after
plotting; the authors' repo has no script for Fig. 1h-j). Arrows are transcribed by eye.

Usage (from Derived_Modalities/):
    # version 1: panel j on RNA
    python 01_data_preprocess/kuppe_fig1/fig1hij.py data/kuppe_2022 figures/kuppe_fig1h \
        results/kuppe_2022/exp2_joint_harmony/exp2_labels.csv.gz kuppe_fig1hij_rna RNA
    # version 2: panel j on CollecTRI TF activity (h and i identical)
    python 01_data_preprocess/kuppe_fig1/fig1hij.py data/kuppe_2022 figures/kuppe_fig1h \
        results/kuppe_2022/exp2_joint_harmony_tf/exp2_labels.csv.gz kuppe_fig1hij_tf \
        "TF activity"
"""

import json
import os
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import EllipseCollection
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from scipy import ndimage
from scipy.optimize import linear_sum_assignment
from scipy.spatial import cKDTree

DATA, OUT, SL_LABELS = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
# version of panel j: output name and the layer named in its title
NAME = sys.argv[4] if len(sys.argv) > 4 else "kuppe_fig1hij"
LAYER = sys.argv[5] if len(sys.argv) > 5 else "RNA"
OUT.mkdir(parents=True, exist_ok=True)
MM = 1 / 25.4

mpl.rcParams.update(
    {
        "font.family": "Arial",
        "font.size": 6,
        "axes.linewidth": 0.5,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)

ROWS = [
    ("ACH003", "Ctrl tissue (P7)"),
    ("10X0026", "BZ tissue (P3)"),
    ("10X0017", "IZ tissue (P3)"),
    ("ACH005", "FZ tissue (P14)"),
]
CMAP = {
    "h": mpl.colormaps["inferno"],
    # blue -> magenta -> yellow, matched by eye to the published panel i
    "i": LinearSegmentedColormap.from_list(
        "fig1i", ["#2a22d8", "#7b2fd8", "#c13fc0", "#f0698f", "#fca44b", "#fdec3a"]),
}
TITLE = {
    "h": "Cell-type deconvolution (cell2location)",
    "i": "Pathway activity (PROGENy)",
    "j": f"SpatialLeiden domains ({LAYER}, Harmony, all 27 slides jointly; "
         + ("K free" if os.environ.get("FIG_FREE_K") else "K = 9") + ")",
}
SL_RATIOS = [float(r) if float(r) else 0 for r in os.environ.get("FIG_SL_RATIOS", "0,0.5,1.0,1.5").split(",")]
SL_SEED = 0
K = 9
# 9 domain colours, distinct on black (Okabe-Ito + two)
DOMAIN_COL = ["#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00",
              "#CC79A7", "#FFFFFF", "#8C6BB1",
              # extra domains when K is not fixed (FIG_FREE_K): further distinct colours
              "#A6761D", "#66C2A5", "#E7298A", "#B3DE69", "#FB8072", "#80B1D3", "#BEBADA",
              "#FFED6F", "#1B9E77", "#7570B3", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FDB462"]

# arrows: tip (x, y) as fractions of the slide's spot extent (y downwards) and direction
# (dx, dy) from tail to tip as a fraction of the frame; transcribed from the published panel
R, D, U = (0.15, 0.06), (0.12, 0.12), (0.04, -0.15)
# one entry per column: (block, {run: (label, source, feature, arrows)})
COLUMNS = [
    ("h", {
        "ACH003": ("vCMs", "c2l", "CM", [(0.57, 0.21, *R), (0.41, 0.53, *R)]),
        "10X0026": ("vCMs", "c2l", "CM", [(0.48, 0.55, *R)]),
        "10X0017": ("vCMs", "c2l", "CM", [(0.42, 0.51, *R)]),
        "ACH005": ("vCMs", "c2l", "CM", [(0.44, 0.35, *D)]),
    }),
    ("h", {
        "ACH003": ("Fibroblast", "c2l", "Fib", [(0.66, 0.33, *R), (0.42, 0.71, *R)]),
        "10X0026": ("Fibroblast", "c2l", "Fib", [(0.53, 0.55, *R)]),
        "10X0017": ("Fibroblast", "c2l", "Fib", [(0.23, 0.30, *R)]),
        "ACH005": ("Fibroblast", "c2l", "Fib", [(0.35, 0.56, *D)]),
    }),
    ("h", {
        "ACH003": ("Myeloid", "c2l", "Myeloid", [(0.65, 0.33, *R), (0.43, 0.69, *R)]),
        "10X0026": ("Myeloid", "c2l", "Myeloid", [(0.57, 0.54, *R)]),
        "10X0017": ("Myeloid", "c2l", "Myeloid", [(0.97, 0.46, *R)]),
        "ACH005": ("Myeloid", "c2l", "Myeloid", [(0.42, 0.55, *D)]),
    }),
    ("i", {
        "ACH003": ("TGFβ", "progeny", "TGFb", [(0.64, 0.36, *R), (0.39, 0.69, *R)]),
        "10X0026": ("TGFβ", "progeny", "TGFb", [(0.31, 0.51, *R), (0.59, 0.51, *R)]),
        "10X0017": ("TGFβ", "progeny", "TGFb", [(0.24, 0.28, *R), (0.85, 0.25, *R)]),
        "ACH005": ("TGFβ", "progeny", "TGFb", [(0.77, 0.40, *R), (0.28, 0.56, *R)]),
    }),
    ("i", {
        "ACH003": ("NFκB", "progeny", "NFkB", [(0.43, 0.68, *R)]),
        "10X0026": ("NFκB", "progeny", "NFkB", [(0.59, 0.51, *R)]),
        "10X0017": ("NFκB", "progeny", "NFkB", [(0.85, 0.22, *R), (0.90, 0.35, *R)]),
        "ACH005": ("NFκB", "progeny", "NFkB", [(0.34, 0.61, *D), (0.75, 0.12, *U)]),
    }),
] + [
    ("j", {sid: (f"ratio {r:g}", "sl", f"r{r}_s{SL_SEED}", []) for sid, _ in ROWS})
    for r in SL_RATIOS
]


def matched_domains(path):
    """SpatialLeiden labels of the plotted runs as D1, D2, ..., consistent across ratios.
    D1.. by size at the first ratio; at the other ratios each domain takes the number of the
    first-ratio domain it overlaps most (Hungarian). With K free, domains left without a
    partner get new numbers (largest first) after the first ratio's."""
    lab = pd.read_csv(path, index_col=0, dtype=str)
    cols = [f"r{r}_s{SL_SEED}" for r in SL_RATIOS]
    ref = lab[cols[0]]
    order = ref.value_counts().index  # D1 = largest domain at ratio 0
    names = {c: f"D{i + 1}" for i, c in enumerate(order)}
    out = pd.DataFrame({cols[0]: ref.map(names)}, index=lab.index)
    nxt = len(order)
    for c in cols[1:]:
        tab = pd.crosstab(lab[c], out[cols[0]])
        rows, cols_ = linear_sum_assignment(-tab.to_numpy())
        m = {tab.index[r]: tab.columns[k] for r, k in zip(rows, cols_)}
        for extra in lab[c][~lab[c].isin(m)].value_counts().index:
            nxt += 1
            m[extra] = f"D{nxt}"
        out[c] = lab[c].map(m)
    return out


SL = matched_domains(SL_LABELS)
DOMAINS = sorted(set(SL.to_numpy().ravel()), key=lambda d: int(d[1:]))
if len(DOMAINS) > len(DOMAIN_COL):  # free K at high weights: many small domains
    DOMAIN_COL = DOMAIN_COL + [mpl.colors.to_hex(c) for c in
                               mpl.colormaps["tab20b"].colors + mpl.colormaps["tab20c"].colors]
DOMAIN_COL = DOMAIN_COL + ["#999999"] * max(0, len(DOMAINS) - len(DOMAIN_COL))
DCOL = dict(zip(DOMAINS, DOMAIN_COL))
# legend: only domains present on the plotted slides
_shown = SL.loc[SL.index.str.split("_", n=1).str[0].isin([sid for sid, _ in ROWS])]
DOMAINS = [d for d in DOMAINS if (_shown == d).to_numpy().any()]


def load_values(sid, barcodes):
    """All features of one slide, as frames aligned to `barcodes`."""
    spots = pd.read_csv(DATA / "visium_tables" / f"{sid}_spots.csv", index_col="barcode")
    c2l = spots.filter(like="prop_").rename(columns=lambda c: c.removeprefix("prop_"))
    prog = pd.read_csv(DATA / "modalities" / "progeny" / f"progeny_{sid}.csv",
                       index_col="barcode")
    sl_ = SL.loc[SL.index.str.startswith(f"{sid}_")]
    sl_.index = sl_.index.str.removeprefix(f"{sid}_")
    out = {"c2l": c2l, "progeny": prog, "sl": sl_}
    for k, v in out.items():
        assert set(v.index) == set(barcodes), (sid, k)
        out[k] = v.loc[barcodes]
    return out


def outline(ax, x, y, pitch):
    """Dashed white contour around the tissue: rasterise spots, close holes, trace edge."""
    step = pitch / 4
    gx = ((x - x.min()) / step).round().astype(int) + 8
    gy = ((y - y.min()) / step).round().astype(int) + 8
    mask = np.zeros((gy.max() + 9, gx.max() + 9), bool)
    mask[gy, gx] = True
    mask = ndimage.binary_dilation(mask, iterations=int(np.ceil(0.9 * pitch / step)))
    mask = ndimage.binary_fill_holes(mask)
    xs = x.min() + (np.arange(mask.shape[1]) - 8) * step
    ys = y.min() + (np.arange(mask.shape[0]) - 8) * step
    ax.contour(xs, ys, mask.astype(float), levels=[0.5],
               colors="white", linewidths=0.6, linestyles=[(0, (2.5, 1.5))])


SLIDES = {}
for sid, _ in ROWS:
    df = pd.read_csv(DATA / "visium_tables" / f"{sid}_spots.csv")
    x, y = df.pxl_col_fullres.to_numpy(float), df.pxl_row_fullres.to_numpy(float)
    pitch = np.median(cKDTree(np.c_[x, y]).query(np.c_[x, y], k=2)[0][:, 1])
    SLIDES[sid] = (x, y, pitch, load_values(sid, df.barcode.to_numpy()))


def draw(columns, name, width_mm):
    ncol, left = len(columns), 0.06 * 89 / width_mm
    cell_mm = width_mm * (0.995 - left) / (ncol + 0.03 * (ncol - 1))  # square panels
    height_mm = (len(ROWS) + 0.03 * (len(ROWS) - 1)) * cell_mm / (0.955 - 0.05)
    fig, axes = plt.subplots(
        len(ROWS), ncol, figsize=(width_mm * MM, height_mm * MM),
        gridspec_kw={"wspace": 0.03, "hspace": 0.03, "left": left,
                     "right": 0.995, "top": 0.955, "bottom": 0.05},
    )
    summary = {}
    for r, (sid, row_label) in enumerate(ROWS):
        x, y, pitch, vals = SLIDES[sid]
        pad = 1.6 * pitch
        x0, x1, y0, y1 = x.min() - pad, x.max() + pad, y.min() - pad, y.max() + pad
        cx, cy, half = (x0 + x1) / 2, (y0 + y1) / 2, max(x1 - x0, y1 - y0) / 2
        summary[sid] = {"row": row_label, "spots": int(len(x))}
        for c, (block, spec) in enumerate(columns):
            label, source, feat, arrows = spec[sid]
            ax = axes[r, c]
            v = vals[source][feat].to_numpy()
            if source == "sl":
                colors = [DCOL[d] for d in v]
                summary[sid][f"{block}:{label}"] = {
                    d: round(float((v == d).mean()), 3) for d in DOMAINS if (v == d).any()}
            else:
                lo, hi = np.quantile(v, [0.01, 0.99])  # Seurat min/max.cutoff q1 / q99
                colors = CMAP[block](np.clip((v - lo) / (hi - lo), 0, 1))
                summary[sid][f"{block}:{label}"] = {
                    "source": source, "q01": round(float(lo), 4), "q99": round(float(hi), 4)}
            ax.set_facecolor("black")
            ax.add_collection(EllipseCollection(
                np.full(len(x), pitch), np.full(len(x), pitch), np.zeros(len(x)),
                units="xy", offsets=np.c_[x, y], offset_transform=ax.transData,
                facecolors=colors, linewidths=0,
            ))
            outline(ax, x, y, pitch)
            for fx, fy, dx, dy in arrows:
                tip = (x.min() + fx * (x.max() - x.min()), y.min() + fy * (y.max() - y.min()))
                tail = (tip[0] - dx * 2 * half, tip[1] - dy * 2 * half)
                ax.annotate("", xy=tip, xytext=tail, arrowprops={
                    "arrowstyle": "-|>,head_length=0.5,head_width=0.25",
                    "lw": 1.0, "color": "white"})
            ax.set_xlim(cx - half, cx + half)
            ax.set_ylim(cy + half, cy - half)
            ax.set_aspect("equal")
            ax.set_xticks([])
            ax.set_yticks([])
            for s in ax.spines.values():
                s.set_visible(False)
            ax.text(0.03, 0.97, label, transform=ax.transAxes, color="white",
                    ha="left", va="top", fontsize=5,
                    bbox={"facecolor": "black", "edgecolor": "none", "pad": 0.8})
            if c == 0:
                ax.set_ylabel(row_label, fontsize=5.5, labelpad=2)

    # Min-Max colour bar under each continuous column; one domain legend under block j
    for c, (block, _) in enumerate(columns):
        if block == "j":
            continue
        pos = axes[-1, c].get_position()
        cax = fig.add_axes([pos.x0 + 0.55 * pos.width, 0.022, 0.4 * pos.width, 0.008])
        cax.imshow(np.linspace(0, 1, 256)[None], cmap=CMAP[block], aspect="auto")
        cax.set_xticks([0, 255], ["Min", "Max"], fontsize=4.5)
        cax.tick_params(length=0, pad=1)
        cax.set_yticks([])
        for s in cax.spines.values():
            s.set_visible(False)
    for block in dict.fromkeys(b for b, _ in columns):
        idx = [c for c, (b, _) in enumerate(columns) if b == block]
        left_ax, right_ax = axes[0, idx[0]].get_position(), axes[0, idx[-1]].get_position()
        first = idx[0] == 0
        fig.text((left_ax.x0 + right_ax.x1) / 2, 0.975, TITLE[block], ha="center", va="bottom")
        fig.lines.append(plt.Line2D([left_ax.x0 + (0 if first else 0.015), right_ax.x1],
                                    [0.972, 0.972], color="black", lw=0.5,
                                    transform=fig.transFigure))
        fig.text(0.0 if first else left_ax.x0, 0.995, block,
                 fontsize=8, fontweight="bold", va="top")
        if block == "j":
            bottom = axes[-1, idx[-1]].get_position()
            handles = [Line2D([], [], ls="", marker="o", ms=4, mfc=DCOL[d], mec="black",
                              mew=0.3, label=d) for d in DOMAINS]
            fig.legend(handles=handles, ncol=min(len(DOMAINS), 12), frameon=False, fontsize=5,
                       loc="lower center", handletextpad=0.1, columnspacing=0.9,
                       bbox_to_anchor=((left_ax.x0 + bottom.x1) / 2, 0.0))
    fig.savefig(OUT / f"{name}.pdf")
    fig.savefig(OUT / f"{name}.png", dpi=450)
    plt.close(fig)
    return summary


draw(COLUMNS[:3], "kuppe_fig1h", 89)
summary = draw(COLUMNS, NAME, 183 * 9 / 7)
(OUT / f"{NAME}_summary.json").write_text(json.dumps(summary, indent=1, ensure_ascii=False))
print(json.dumps(summary, indent=1, ensure_ascii=False))
