"""P9 IZ rep2 (ACH0013): the diagonal domain found by SpatialLeiden (RNA + TF + space) is an
artery wall. One row: H&E | domains (band highlighted) | MYH11 | GJA5 | MYOCD activity |
JAG1 -> NOTCH3 CCC | authors' cell-type niches. The band is outlined on every panel.

Band = domain 6 of mm_proto_rna_tf, label sp0.75_tf0.5_s0 (RNA 1, TF 0.5, spatial 0.75,
K = 9 authors' cell-type niches, seed 0); it matches the pink RNA + space domain in
fig_rna_tf_combo. Genes: log1p(CP10k). MYOCD: decoupler ULM + CollecTRI score. JAG1^NOTCH3:
team CCC local score. Colour scales 1st-99th percentile of the slide.

Usage (from Derived_Modalities/):
    python 04_multimodal_weight_grid/proto_rna_tf/fig_vessel_band.py data/kuppe_2022 results/kuppe_2022
"""

import json
import sys
from pathlib import Path

import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import EllipseCollection, LineCollection
from PIL import Image
from scipy.spatial import cKDTree

DATA, RES = Path(sys.argv[1]), Path(sys.argv[2])
OUT = RES / "mm_proto_rna_tf" / "figures"
RUN, KEY, BAND = "ACH0013", "sp0.75_tf0.5_s0", "6"
MM = 1 / 25.4
mpl.rcParams.update({"font.family": "Arial", "font.size": 6, "axes.titlesize": 6.5, "pdf.fonttype": 42})
NICHE_COL = {f"ctniche_{i}": c for i, c in enumerate(
    ["#D51F26", "#272E6A", "#208A42", "#89288F", "#F47D2B", "#FEE500", "#8A9FD1", "#C06CAB", "#D8A767"], start=1)}

lab = pd.read_csv(RES / "mm_proto_rna_tf" / f"labels_{RUN}.csv.gz", index_col=0, dtype=str)[KEY]
idx = lab.index
sp = pd.read_csv(DATA / "visium_tables" / f"{RUN}_spots.csv").set_index("barcode").loc[idx]
xy = sp[["pxl_col_fullres", "pxl_row_fullres"]].to_numpy(float)
band = (lab == BAND).to_numpy()

# outline: hexagon edges between a band spot and any grid neighbour outside the band
grid = np.c_[sp.array_col, sp.array_row].astype(float)
A = np.c_[grid, np.ones(len(grid))]
coef, *_ = np.linalg.lstsq(A, xy, rcond=None)  # array (col, row) -> pixel, affine
pos = {(int(r), int(c)): n for n, (r, c) in enumerate(zip(sp.array_row, sp.array_col))}
segs = []
for n in np.flatnonzero(band):
    r, c = int(sp.array_row.iloc[n]), int(sp.array_col.iloc[n])
    for dr, dc in [(0, 2), (0, -2), (1, 1), (1, -1), (-1, 1), (-1, -1)]:
        m = pos.get((r + dr, c + dc))
        if m is not None and band[m]:
            continue
        v = np.array([dc, dr]) @ coef[:2]          # pixel vector to that neighbour
        mid = xy[n] + v / 2
        perp = np.array([-v[1], v[0]]) / np.linalg.norm(v) * np.linalg.norm(v) / (2 * np.sqrt(3))
        segs.append([mid - perp, mid + perp])
segs = np.array(segs)

# features
h = ad.read_h5ad(DATA / "processed1" / "Visium_GT_IZ_P9_rep2.h5ad")
counts = h.raw.X
tot = np.asarray(counts.sum(1)).ravel()
names = list(h.raw.var_names)
def gene(g):
    v = np.log1p(np.asarray(counts[:, names.index(g)].todense()).ravel() / tot * 1e4)
    return pd.Series(v, index=h.obs_names).loc[idx].to_numpy()
tf = ad.read_h5ad(DATA / "modalities" / f"collectri_ulm_{RUN}.h5ad")
myocd = pd.Series(np.asarray(tf[:, "MYOCD"].X).ravel(), index=tf.obs_names).loc[idx].to_numpy()
ccc = pd.read_csv(DATA / "modalities" / f"ccc_{RUN}.csv", index_col=0).loc[idx, "JAG1^NOTCH3"].to_numpy()
FEAT = [("MYH11\nsmooth muscle", gene("MYH11"), "log1p(CP10k)", "magma"),
        ("GJA5 (Cx40)\nartery", gene("GJA5"), "log1p(CP10k)", "magma"),
        ("MYOCD\nTF activity", myocd, "ULM score", "RdBu_r"),
        ("JAG1 → NOTCH3\nCCC", ccc, "local score", "viridis")]

PAL = [c for i, c in enumerate(mpl.colormaps["tab20"].colors) if i not in (14, 15)]


def draw(ax, colours):
    pitch = np.median(cKDTree(xy).query(xy, k=2)[0][:, 1])
    d = np.full(len(xy), 1.05 * pitch)
    coll = EllipseCollection(d, d, np.zeros(len(xy)), units="xy", offsets=xy, offset_transform=ax.transData,
                             facecolors=colours, linewidths=0)
    coll.set_rasterized(True)
    ax.add_collection(coll)
    ax.set_xlim(xy[:, 0].min() - pitch, xy[:, 0].max() + pitch)
    ax.set_ylim(xy[:, 1].max() + pitch, xy[:, 1].min() - pitch)
    ax.set_aspect("equal")
    ax.axis("off")


def outline(ax, scale=1.0, colour="black"):
    ax.add_collection(LineCollection(segs * scale, colors=colour, linewidths=0.6))


fig, axes = plt.subplots(1, 7, figsize=(183 * MM, 40 * MM), gridspec_kw={"wspace": 0.05})
ax = axes[0]
img = np.asarray(Image.open(DATA / "visium_tables" / f"{RUN}_lowres.png").convert("RGB"))
f = json.loads((DATA / "visium_tables" / f"{RUN}_scalefactors.json").read_text())["tissue_lowres_scalef"]
ax.imshow(img)
ax.set_xlim(xy[:, 0].min() * f - 8, xy[:, 0].max() * f + 8)
ax.set_ylim(xy[:, 1].max() * f + 8, xy[:, 1].min() * f - 8)
ax.axis("off")
outline(ax, f, "black")
ax.set_title("H&E", fontweight="bold")
# domains: band in its fig_rna_tf_combo colour, others grey
draw(axes[1], [PAL[7] if b else (0.82, 0.82, 0.82, 1) for b in band])
outline(axes[1])
axes[1].set_title("RNA + TF\n+ space", fontweight="bold")
for j, (title, v, unit, cm) in enumerate(FEAT, start=2):
    if cm == "RdBu_r":
        lim = float(np.quantile(np.abs(v), 0.99)); norm = mpl.colors.Normalize(-lim, lim, clip=True)
    else:
        norm = mpl.colors.Normalize(float(np.quantile(v, 0.01)), float(np.quantile(v, 0.99)), clip=True)
    draw(axes[j], mpl.colormaps[cm](norm(v)))
    outline(axes[j], colour="white" if cm in ("magma", "viridis") else "black")
    axes[j].set_title(title, fontweight="bold")
    cax = axes[j].inset_axes([0.1, -0.1, 0.8, 0.05])
    cb = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cm), cax=cax, orientation="horizontal")
    cb.ax.tick_params(labelsize=5, length=1.5, width=0.4)
    cb.set_label(unit, fontsize=5.5)
    cb.outline.set_linewidth(0.4)
draw(axes[6], [NICHE_COL[v] for v in sp.celltype_niche])
outline(axes[6])
axes[6].set_title("Authors'\nniches", fontweight="bold")
axes[0].text(-0.06, 0.5, "P9 IZ rep2\n(ACH0013)", transform=axes[0].transAxes, rotation=90, ha="right",
             va="center", fontweight="bold", color="#D55E00", fontsize=7)
fig.savefig(OUT / "fig_p9rep2_vessel_band.png", dpi=300, bbox_inches="tight")
fig.savefig(OUT / "fig_p9rep2_vessel_band.pdf", dpi=300, bbox_inches="tight")
stats = pd.DataFrame([{"feature": t.split("\n")[0], "mean_in_band": v[band].mean(), "mean_outside": v[~band].mean()}
                      for t, v, _, _ in FEAT]).round(3)
stats.to_csv(OUT / "fig_p9rep2_vessel_band_stats.csv", index=False)
print(f"band spots {band.sum()}, outline segments {len(segs)}")
print(stats.to_string(index=False))
