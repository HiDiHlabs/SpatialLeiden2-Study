"""Maps of selected TF activities, z-scored within each slide (spatial pattern), from the
saved collectri_ulm_<run>.h5ad files. Raw ULM scores are mostly positive on every slide,
so a shared raw scale hides within-slide structure; the z-score shows where a TF is more or
less active than elsewhere on the same slide.

Usage (from Derived_Modalities/): python 02_modalities/collectri/collectri_maps_zscore.py data/kuppe_2022 results/kuppe_2022/collectri
"""

import sys
from pathlib import Path

import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import EllipseCollection
from scipy.spatial import cKDTree

DATA, RES = Path(sys.argv[1]), Path(sys.argv[2])
SLIDES = {
    "10X001": "P1 control (10X001)",
    "ACH0013": "P9 IZ rep2 (ACH0013)",
    "ACH0012": "P9 IZ (ACH0012)",
}
SHOW = ["NFKB1", "STAT3", "SMAD3", "HIF1A", "SPI1", "MEF2C", "GATA4", "TEAD1"]
mpl.rcParams.update(
    {"font.family": "Arial", "font.size": 6, "axes.titlesize": 6.5, "pdf.fonttype": 42}
)
objs = {
    r: ad.read_h5ad(DATA / "modalities" / f"collectri_ulm_{r}.h5ad") for r in SLIDES
}
fig, axes = plt.subplots(
    len(SLIDES),
    len(SHOW),
    figsize=(183 / 25.4, 1.7 * len(SLIDES)),
    gridspec_kw={"wspace": 0.05, "hspace": 0.15},
)
norm = mpl.colors.TwoSlopeNorm(0, -2.5, 2.5)
for i, (run, label) in enumerate(SLIDES.items()):
    o = objs[run]
    x, y = o.obsm["spatial"].T
    pitch = np.median(cKDTree(np.c_[x, y]).query(np.c_[x, y], k=2)[0][:, 1])
    d = np.full(len(x), 1.05 * pitch)
    for j, tf in enumerate(SHOW):
        ax = axes[i, j]
        v = o[:, tf].X.ravel()
        z = (v - v.mean()) / v.std()
        coll = EllipseCollection(
            d,
            d,
            np.zeros(len(x)),
            units="xy",
            offsets=np.c_[x, y],
            offset_transform=ax.transData,
            facecolors=mpl.colormaps["RdBu_r"](norm(z)),
            linewidths=0,
        )
        coll.set_rasterized(True)
        ax.add_collection(coll)
        ax.set_xlim(x.min() - pitch, x.max() + pitch)
        ax.set_ylim(y.max() + pitch, y.min() - pitch)
        ax.set_aspect("equal")
        ax.axis("off")
        if i == 0:
            ax.set_title(tf, fontweight="bold")
        if j == 0:
            ax.text(
                -0.08,
                0.5,
                label,
                transform=ax.transAxes,
                rotation=90,
                ha="right",
                va="center",
            )
cax = fig.add_axes([0.35, 0.06, 0.3, 0.012])
cb = fig.colorbar(
    mpl.cm.ScalarMappable(norm=norm, cmap="RdBu_r"), cax=cax, orientation="horizontal"
)
cb.set_label("TF activity, z-score within slide (decoupler ULM, CollecTRI)")
fig.savefig(
    RES / "figures" / "collectri_selected_tfs_zscore.png", dpi=250, bbox_inches="tight"
)
fig.savefig(
    RES / "figures" / "collectri_selected_tfs_zscore.pdf", dpi=300, bbox_inches="tight"
)
print("written")
