"""Maps of cell-cell communication (CCC) scores on P1 and the three P9 sections.

Input: spot x ligand-receptor score tables provided by the team (LIANA+-style local
scores, columns "LIGAND^RECEPTOR", complexes joined by "_", values 0-1), stored as
data/kuppe_2022/modalities/ccc_<run>.csv. Spots match the authors' spot tables exactly.
How the scores were computed (metric, bandwidth, pair filtering) is not recorded here yet.

Rows: P1 control (10X001), P9 remote zone (ACH0023), P9 IZ (ACH0012), P9 IZ rep2 (ACH0013).
Columns: selected pairs relevant to the infarct biology, plus the mean over the pairs shared
by all four slides. One colour scale per column, shared across slides (0 to the 99th
percentile). A pair absent from a slide's table (filtered upstream) is marked, not imputed.

Usage (from Derived_Modalities/): python 02_modalities/ccc/ccc_maps.py data/kuppe_2022 results/kuppe_2022/ccc
"""

import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import EllipseCollection
from scipy.spatial import cKDTree

DATA, OUT = Path(sys.argv[1]), Path(sys.argv[2])
(OUT / "figures").mkdir(parents=True, exist_ok=True)
MM = 1 / 25.4
mpl.rcParams.update({"font.family": "Arial", "font.size": 6, "axes.titlesize": 6.5, "pdf.fonttype": 42})
SLIDES = {"10X001": "P1 control\n(10X001)", "ACH0023": "P9 RZ\n(ACH0023)",
          "ACH0012": "P9 IZ\n(ACH0012)", "ACH0013": "P9 IZ rep2\n(ACH0013)"}
SLIDE_COL = {"10X001": "#0072B2", "ACH0023": "#009E73", "ACH0012": "#E69F00", "ACH0013": "#D55E00"}
PAIRS = [("SPP1^CD44", "SPP1 → CD44\nSPP1+ macrophages"),
         ("FN1^ITGAV_ITGB1", "FN1 → ITGAV/ITGB1\nwound matrix"),
         ("THBS1^CD36", "THBS1 → CD36\ninjury response"),
         ("APP^CD74", "APP → CD74\nMHC-II / macrophages"),
         ("C3^LRP1", "C3 → LRP1\ncomplement"),
         ("VEGFB^NRP1", "VEGFB → NRP1\nvascular / myocardial")]

tabs = {r: pd.read_csv(DATA / "modalities" / f"ccc_{r}.csv", index_col=0) for r in SLIDES}
common = sorted(set.intersection(*(set(t.columns) for t in tabs.values())))
for r, t in tabs.items():
    t["__mean_shared__"] = t[common].mean(axis=1)
cols = PAIRS + [("__mean_shared__", f"Mean of {len(common)}\nshared pairs")]
vmax = {p: float(np.quantile(np.concatenate([t[p].to_numpy() for t in tabs.values() if p in t]), 0.99))
        for p, _ in cols}

fig, axes = plt.subplots(len(SLIDES), len(cols), figsize=(183 * MM, 34 * MM * len(SLIDES)),
                         gridspec_kw={"wspace": 0.05, "hspace": 0.1})
summary = []
for i, (run, label) in enumerate(SLIDES.items()):
    sp = pd.read_csv(DATA / "visium_tables" / f"{run}_spots.csv").set_index("barcode").loc[tabs[run].index]
    xy = sp[["pxl_col_fullres", "pxl_row_fullres"]].to_numpy(float)
    pitch = np.median(cKDTree(xy).query(xy, k=2)[0][:, 1])
    d = np.full(len(xy), 1.05 * pitch)
    for j, (pair, title) in enumerate(cols):
        ax = axes[i, j]
        ax.set_xlim(xy[:, 0].min() - pitch, xy[:, 0].max() + pitch)
        ax.set_ylim(xy[:, 1].max() + pitch, xy[:, 1].min() - pitch)
        ax.set_aspect("equal")
        ax.axis("off")
        if pair in tabs[run]:
            v = tabs[run][pair].to_numpy()
            norm = mpl.colors.Normalize(0, vmax[pair], clip=True)
            colours = mpl.colormaps["viridis"](norm(v))
            summary.append({"slide": run, "pair": pair, "mean": round(float(v.mean()), 3),
                            "share_spots_gt0": round(float((v > 0).mean()), 3)})
        else:
            colours = [(0.88, 0.88, 0.88, 1)] * len(xy)
            ax.text(0.5, 0.5, "not in this\nslide's results", transform=ax.transAxes, ha="center", va="center",
                    fontsize=5.5, color="#444444")
            summary.append({"slide": run, "pair": pair, "mean": np.nan, "share_spots_gt0": np.nan})
        coll = EllipseCollection(d, d, np.zeros(len(xy)), units="xy", offsets=xy, offset_transform=ax.transData,
                                 facecolors=colours, linewidths=0)
        coll.set_rasterized(True)
        ax.add_collection(coll)
        if i == 0:
            ax.set_title(title, fontsize=6)
        if j == 0:
            ax.text(-0.06, 0.5, label, transform=ax.transAxes, rotation=90, ha="right", va="center",
                    fontweight="bold", color=SLIDE_COL[run], fontsize=6.5)
for j, (pair, _) in enumerate(cols):
    cax = axes[-1, j].inset_axes([0.1, -0.12, 0.8, 0.05])
    cb = fig.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(0, vmax[pair]), cmap="viridis"), cax=cax,
                      orientation="horizontal")
    cb.ax.tick_params(labelsize=4.5, length=1.5, width=0.4)
    cb.outline.set_linewidth(0.4)
fig.savefig(OUT / "figures" / "ccc_selected_pairs.png", dpi=250, bbox_inches="tight")
fig.savefig(OUT / "figures" / "ccc_selected_pairs.pdf", dpi=300, bbox_inches="tight")
pd.DataFrame(summary).to_csv(OUT / "ccc_selected_pairs_summary.csv", index=False)
print(pd.DataFrame(summary).pivot(index="pair", columns="slide", values="mean").round(3).to_string())
