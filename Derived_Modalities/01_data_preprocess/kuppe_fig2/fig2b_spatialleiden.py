"""Kuppe et al. 2022, Fig. 2b with SpatialLeiden: the authors' cell-type niches on P1 (control)
and P9 (IZ), next to the joint SpatialLeiden domains (all 27 slides, K = 9) at increasing
spatial weight.

Columns: authors' cell-type niches (as in fig2b.py, their colours), then SpatialLeiden seed 0
at the weights in FIG_SL_RATIOS (default 0, 0.5, 1, 1.5). Rows: P1 control (10X001) and
P9 IZ rep2 (ACH0013), the published slides (see fig2b.py).

Domain numbering and colours are the same as in 01_data_preprocess/kuppe_fig1/fig1hij.py:
D1-D9 by size at weight 0 over all 88,704 spots; the other weights matched to weight 0 by
maximum overlap (Hungarian). So a domain has the same colour in Fig. 1 and Fig. 2.

Usage (from Derived_Modalities/):
    python 01_data_preprocess/kuppe_fig2/fig2b_spatialleiden.py data/kuppe_2022 figures/kuppe_fig2b \
        results/kuppe_2022/exp2_joint_harmony/exp2_labels.csv.gz kuppe_fig2b_spatialleiden_rna RNA
"""

import os
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import EllipseCollection
from matplotlib.lines import Line2D
from scipy.optimize import linear_sum_assignment
from scipy.spatial import cKDTree

DATA, OUT, SL_LABELS = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
NAME = sys.argv[4] if len(sys.argv) > 4 else "kuppe_fig2b_spatialleiden"
LAYER = sys.argv[5] if len(sys.argv) > 5 else "RNA"
OUT.mkdir(parents=True, exist_ok=True)
MM = 1 / 25.4
mpl.rcParams.update({"font.family": "Arial", "font.size": 6, "axes.titlesize": 6,
                     "legend.fontsize": 5.5, "pdf.fonttype": 42, "ps.fonttype": 42})

SL_RATIOS = [float(r) if float(r) else 0
             for r in os.environ.get("FIG_SL_RATIOS", "0,0.5,1.0,1.5").split(",")]
SEED = 0
ROWS = [("10X001", "P1 (control)"), ("ACH0013", "P9 (IZ)")]
NICHE_COL = {1: "#D51F26", 2: "#272E6A", 3: "#208A42", 4: "#89288F", 5: "#F47D2B",
             6: "#FEE500", 7: "#8A9FD1", 8: "#C06CAB", 9: "#D8A767"}  # authors' (fig2b.py)
DOMAIN_COL = ["#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00",
              "#CC79A7", "#FFFFFF", "#8C6BB1"]  # as fig1hij.py
K = 9
DOMAINS = [f"D{i + 1}" for i in range(K)]
DCOL = dict(zip(DOMAINS, DOMAIN_COL))


def matched_domains(path):
    """Same numbering as fig1hij.py: D1..D9 by size at weight 0, others matched to it."""
    lab = pd.read_csv(path, index_col=0, dtype=str)
    cols = [f"r{r}_s{SEED}" for r in SL_RATIOS]
    ref = lab[cols[0]]
    names = {c: f"D{i + 1}" for i, c in enumerate(ref.value_counts().index)}
    out = pd.DataFrame({cols[0]: ref.map(names)}, index=lab.index)
    for c in cols[1:]:
        tab = pd.crosstab(lab[c], out[cols[0]])
        r_, c_ = linear_sum_assignment(-tab.to_numpy())
        out[c] = lab[c].map({tab.index[a]: tab.columns[b] for a, b in zip(r_, c_)})
        assert out[c].notna().all(), c
    return out


SL = matched_domains(SL_LABELS)
ncol = 1 + len(SL_RATIOS)
fig, axes = plt.subplots(len(ROWS), ncol, figsize=(183 * MM, 183 * MM * len(ROWS) / ncol + 22 * MM),
                         gridspec_kw={"wspace": 0.04, "hspace": 0.06, "left": 0.04, "right": 0.995,
                                      "top": 0.88, "bottom": 0.16})
for r, (sid, label) in enumerate(ROWS):
    df = pd.read_csv(DATA / "visium_tables" / f"{sid}_spots.csv")
    x, y = df.pxl_col_fullres.to_numpy(float), df.pxl_row_fullres.to_numpy(float)
    pitch = np.median(cKDTree(np.c_[x, y]).query(np.c_[x, y], k=2)[0][:, 1])
    niche = df.celltype_niche.str.removeprefix("ctniche_").astype(int)
    sl_ = SL.loc[sid + "_" + df.barcode]
    colsets = [("Authors' cell-type niches", niche.map(NICHE_COL).to_list())] + [
        (f"SpatialLeiden, weight {w:g}", [DCOL[d] for d in sl_[f"r{w}_s{SEED}"]]) for w in SL_RATIOS]
    pad = 1.5 * pitch
    x0, x1, y0, y1 = x.min() - pad, x.max() + pad, y.min() - pad, y.max() + pad
    cx, cy, half = (x0 + x1) / 2, (y0 + y1) / 2, max(x1 - x0, y1 - y0) / 2
    for c, (title, colors) in enumerate(colsets):
        ax = axes[r, c]
        d = 1.05 * pitch
        ax.add_collection(EllipseCollection(np.full(len(x), d), np.full(len(x), d), np.zeros(len(x)),
                                            units="xy", offsets=np.c_[x, y],
                                            offset_transform=ax.transData, facecolors=colors,
                                            linewidths=0))
        ax.set_xlim(cx - half, cx + half)
        ax.set_ylim(cy + half, cy - half)
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])
        if c > 0:
            ax.set_facecolor("#d9d9d9")  # light grey so the white domain stays visible
        if r == 0:
            ax.set_title(title, pad=3)
        if c == 0:
            ax.set_ylabel(label)

h1 = [Line2D([], [], ls="", marker="o", ms=4, mfc=NICHE_COL[k], mec="none", label=f"Niche {k}")
      for k in range(1, 10)]
h2 = [Line2D([], [], ls="", marker="o", ms=4, mfc=DCOL[d], mec="black", mew=0.3, label=d)
      for d in DOMAINS]
p0, p1 = axes[-1, 0].get_position(), axes[-1, 1].get_position()
fig.legend(handles=h1, loc="lower center", ncol=3, frameon=False, columnspacing=0.8,
           handletextpad=0.1, bbox_to_anchor=((p0.x0 + p0.x1) / 2, 0.0))
fig.legend(handles=h2, loc="lower center", ncol=9, frameon=False, columnspacing=1.0,
           handletextpad=0.1, bbox_to_anchor=((p1.x0 + axes[-1, -1].get_position().x1) / 2, 0.04))
fig.suptitle(f"Cell-type niches (Kuppe et al.) and SpatialLeiden domains ({LAYER}, Harmony, "
             f"all 27 slides jointly; K = 9; seed {SEED})", y=0.97)
fig.text(0.0, 0.995, "b", fontsize=8, fontweight="bold", va="top")
fig.savefig(OUT / f"{NAME}.pdf")
fig.savefig(OUT / f"{NAME}.png", dpi=450)
print("written", OUT / f"{NAME}.png")
