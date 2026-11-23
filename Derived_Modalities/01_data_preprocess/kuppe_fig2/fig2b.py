"""Reproduce Kuppe et al. 2022, Fig. 2b: cell-type niches on a control and an ischaemic slide.

Uses only shipped values: the authors' `celltype_niche` labels and spot coordinates,
extracted from their Seurat objects by
01_data_preprocess/kuppe_visium_exploration/01_extract_seurat_tables.R.

Slides (identified from the published panel):
    P1 (control) = control_P1, Visium run 10X001
    P9 (IZ)      = GT/IZ_P9_rep2, Visium run ACH0013 (not GT/IZ_P9 = ACH0012)
Niche colours are the authors' (saezlab/visium_heart,
st_snRNAseq/05.1_nicheanalysis/plot_niche_example.R). Arrows point at niche 8
(left) and niche 4 (right), as in the published legend; their positions are
transcribed from the published panel, relative to the slide's spot extent.

Usage (from Derived_Modalities/):
    python 01_data_preprocess/kuppe_fig2/fig2b.py data/kuppe_2022 figures/kuppe_fig2b
"""

import json
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import EllipseCollection
from matplotlib.lines import Line2D
from scipy.spatial import cKDTree

DATA, OUT = Path(sys.argv[1]), Path(sys.argv[2])
OUT.mkdir(parents=True, exist_ok=True)
TABLES = DATA / "visium_tables"
MM = 1 / 25.4

mpl.rcParams.update(
    {
        "font.family": "Arial",
        "font.size": 6,
        "axes.titlesize": 6,
        "legend.fontsize": 6,
        "axes.linewidth": 0.5,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)

NICHE_COL = {
    1: "#D51F26",
    2: "#272E6A",
    3: "#208A42",
    4: "#89288F",
    5: "#F47D2B",
    6: "#FEE500",
    7: "#8A9FD1",
    8: "#C06CAB",
    9: "#D8A767",
}
PANELS = [
    (
        "10X001",
        "P1 (control)",
        8,
        [(0.25, 0.46, 0.05, 0.10), (0.50, 0.48, -0.12, 0.02), (0.46, 0.68, 0.0, -0.12)],
    ),
    (
        "ACH0013",
        "P9 (IZ)",
        4,
        [(0.62, 0.42, 0.01, 0.12), (0.71, 0.60, -0.12, 0.02), (0.53, 0.80, 0.13, 0.03)],
    ),
]
# arrows: tip (x, y) as fractions of the slide's spot extent (y downwards), and the arrow's
# direction (dx, dy) from tail to tip as a fraction of the frame. Each tip is then moved to
# the nearest spot of the named niche, so the arrow points at that niche.

fig, axes = plt.subplots(
    1,
    2,
    figsize=(89 * MM, 52 * MM),
    gridspec_kw={
        "wspace": 0.06,
        "left": 0.02,
        "right": 0.98,
        "top": 0.90,
        "bottom": 0.27,
    },
)
summary = {}
for ax, (sid, title, arrow_niche, arrows) in zip(axes, PANELS):
    df = pd.read_csv(TABLES / f"{sid}_spots.csv")
    niche = df.celltype_niche.str.removeprefix("ctniche_").astype(int)
    x, y = df.pxl_col_fullres.to_numpy(float), df.pxl_row_fullres.to_numpy(float)
    # spots drawn touching (105 % of the centre-to-centre distance), as in the published panel
    pitch = np.median(cKDTree(np.c_[x, y]).query(np.c_[x, y], k=2)[0][:, 1])
    d = 1.05 * pitch
    ax.add_collection(
        EllipseCollection(
            np.full(len(x), d),
            np.full(len(x), d),
            np.zeros(len(x)),
            units="xy",
            offsets=np.c_[x, y],
            offset_transform=ax.transData,
            facecolors=niche.map(NICHE_COL).to_list(),
            linewidths=0,
        )
    )
    pad = 1.5 * pitch  # frame hugs the tissue, as in the published panel
    x0, x1, y0, y1 = x.min() - pad, x.max() + pad, y.min() - pad, y.max() + pad
    side = max(x1 - x0, y1 - y0)
    ax.set_xlim(x0, x1)
    ax.set_ylim(y1, y0)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(title, pad=3)
    for fx, fy, dx, dy in arrows:
        guess = (x.min() + fx * (x.max() - x.min()), y.min() + fy * (y.max() - y.min()))
        target = np.flatnonzero(niche.to_numpy() == arrow_niche)
        nearest = target[
            np.argmin((x[target] - guess[0]) ** 2 + (y[target] - guess[1]) ** 2)
        ]
        tip = (x[nearest], y[nearest])
        tail = (tip[0] - 1.4 * dx * side, tip[1] - 1.4 * dy * side)
        ax.annotate(
            "",
            xy=tip,
            xytext=tail,
            arrowprops={
                "arrowstyle": "-|>,head_length=0.8,head_width=0.4",
                "lw": 1.1,
                "color": "black",
            },
        )
    counts = niche.value_counts().sort_index()
    summary[sid] = {
        "title": title,
        "spots": int(len(df)),
        "niche_fraction": {int(k): round(v / len(df), 3) for k, v in counts.items()},
        "arrow_niche": arrow_niche,
    }

handles = [
    Line2D(
        [],
        [],
        ls="",
        marker="o",
        ms=4.5,
        mfc=NICHE_COL[k],
        mec="none",
        label=f"Niche {k}",
    )
    for k in range(1, 10)
]
# matplotlib fills legend columns top to bottom: 1-3 | 4-6 | 7-9, as in the paper
fig.legend(
    handles=handles,
    loc="lower center",
    ncol=3,
    frameon=False,
    bbox_to_anchor=(0.5, 0.0),
    handletextpad=0.2,
    columnspacing=1.6,
    labelspacing=0.35,
)
fig.text(0.0, 0.985, "b", fontsize=8, fontweight="bold", va="top")
fig.savefig(OUT / "kuppe_fig2b.pdf")
fig.savefig(OUT / "kuppe_fig2b.png", dpi=450)
(OUT / "fig2b_summary.json").write_text(json.dumps(summary, indent=1))
print(json.dumps(summary, indent=1))
