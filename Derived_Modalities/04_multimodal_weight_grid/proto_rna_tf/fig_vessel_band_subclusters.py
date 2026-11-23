"""P9 IZ rep2: nested SpatialLeiden inside the vessel band (vessel_band_subclusters.py).

One row, zoomed on the diagonal band (largest connected piece of domain 6, with margin):
H&E | sub-clusters K = 2 | K = 3 | MYH11 | ACTA2 | COL1A1 | POSTN.
The K = 2 split is artery wall vs a mixed vessel margin: the margin has the same fibroblast
markers as the surrounding granulation tissue (COL1A1, COL1A2, COL3A1, DCN, LUM, POSTN all
equal within 0.2 log units) and differs only by fading smooth-muscle signal (MYH11 wall 2.9,
margin 1.0, ring outside 0.3, rest 0.1). Not a distinct fibrotic layer. Sub-clusters from seed 0
(5/5 seeds give the same split, ARI 0.96). At K = 3 the third cluster is the small islands
off the main band; inside the band K = 3 keeps the K = 2 wall / sheath split. Genes: log1p(CP10k), scale 1st-99th percentile of
the zoomed spots. The band is outlined on every panel; spots outside it are grey in the
sub-cluster panels.

Usage (from Derived_Modalities/):
    python 04_multimodal_weight_grid/proto_rna_tf/fig_vessel_band_subclusters.py data/kuppe_2022 results/kuppe_2022
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
from matplotlib.patches import Patch
from PIL import Image
from scipy.spatial import cKDTree

DATA, RES = Path(sys.argv[1]), Path(sys.argv[2])
OUT = RES / "mm_proto_rna_tf" / "figures"
RUN, KEY, BAND = "ACH0013", "sp0.75_tf0.5_s0", "6"
MM = 1 / 25.4
mpl.rcParams.update({"font.family": "Arial", "font.size": 6, "axes.titlesize": 6.5, "pdf.fonttype": 42})
STEPS = [(0, 2), (0, -2), (1, 1), (1, -1), (-1, 1), (-1, -1)]

lab = pd.read_csv(RES / "mm_proto_rna_tf" / f"labels_{RUN}.csv.gz", index_col=0, dtype=str)[KEY]
idx = lab.index
sp = pd.read_csv(DATA / "visium_tables" / f"{RUN}_spots.csv").set_index("barcode").loc[idx]
xy = sp[["pxl_col_fullres", "pxl_row_fullres"]].to_numpy(float)
band = (lab == BAND).to_numpy()
sub = pd.read_csv(RES / "mm_proto_rna_tf" / "vessel_band" / "vessel_band_sublabels.csv.gz", index_col=0, dtype=str)
pos = {(int(r), int(c)): n for n, (r, c) in enumerate(zip(sp.array_row, sp.array_col))}

# name sub-clusters by their top genes (vessel_band_subclusters.csv) for the legend
summ = pd.read_csv(RES / "mm_proto_rna_tf" / "vessel_band" / "vessel_band_subclusters.csv", dtype={"sub": str})


def name(k, g):
    genes = summ.query("K == @k and sub == @g").top_genes.iloc[0].split(";")
    if "MYH11" in genes or "ACTG2" in genes:
        return "Artery wall (MYH11, ACTG2)"
    if "POSTN" in genes:
        return "Vessel margin (mixed SMC + tissue)"
    if "NPPB" in genes or "TNNT2" in genes:
        return "Margin, K = 3 (minus islands)"
    return "Islands off the band"


SUBCOL = {"Artery wall (MYH11, ACTG2)": "#C2185B", "Vessel margin (mixed SMC + tissue)": "#1E88E5",
          "Margin, K = 3 (minus islands)": "#43A047", "Islands off the band": "#FB8C00"}

# largest connected piece of the band, for the zoom window
seen, comps = set(), []
for n in np.flatnonzero(band):
    if n in seen:
        continue
    stack, comp = [n], []
    seen.add(n)
    while stack:
        q = stack.pop()
        comp.append(q)
        r, c = int(sp.array_row.iloc[q]), int(sp.array_col.iloc[q])
        for dr, dc in STEPS:
            m = pos.get((r + dr, c + dc))
            if m is not None and band[m] and m not in seen:
                seen.add(m)
                stack.append(m)
    comps.append(comp)
main = max(comps, key=len)
pitch = np.median(cKDTree(xy).query(xy, k=2)[0][:, 1])
lo, hi = xy[main].min(0) - 4 * pitch, xy[main].max(0) + 4 * pitch
win = np.all((xy >= lo) & (xy <= hi), axis=1)

# outline segments (as in fig_vessel_band.py)
grid = np.c_[sp.array_col, sp.array_row].astype(float)
coef, *_ = np.linalg.lstsq(np.c_[grid, np.ones(len(grid))], xy, rcond=None)
segs = []
for n in np.flatnonzero(band):
    r, c = int(sp.array_row.iloc[n]), int(sp.array_col.iloc[n])
    for dr, dc in STEPS:
        m = pos.get((r + dr, c + dc))
        if m is not None and band[m]:
            continue
        v = np.array([dc, dr]) @ coef[:2]
        mid = xy[n] + v / 2
        perp = np.array([-v[1], v[0]]) / (2 * np.sqrt(3))
        segs.append([mid - perp, mid + perp])
segs = np.array(segs)

h = ad.read_h5ad(DATA / "processed1" / "Visium_GT_IZ_P9_rep2.h5ad")
counts, tot, names = h.raw.X, np.asarray(h.raw.X.sum(1)).ravel(), list(h.raw.var_names)


def gene(g):
    v = np.log1p(np.asarray(counts[:, names.index(g)].todense()).ravel() / tot * 1e4)
    return pd.Series(v, index=h.obs_names).loc[idx].to_numpy()


GENES = [("MYH11", "smooth muscle"), ("ACTA2", "smooth muscle"), ("COL1A1", "fibroblasts"), ("POSTN", "activated fibroblasts")]


def draw(ax, colours):
    d = np.full(len(xy), 1.05 * pitch)
    coll = EllipseCollection(d, d, np.zeros(len(xy)), units="xy", offsets=xy, offset_transform=ax.transData,
                             facecolors=colours, linewidths=0)
    coll.set_rasterized(True)
    ax.add_collection(coll)
    ax.set_xlim(lo[0], hi[0])
    ax.set_ylim(hi[1], lo[1])
    ax.set_aspect("equal")
    ax.axis("off")


def outline(ax, colour="black", scale=1.0):
    ax.add_collection(LineCollection(segs * scale, colors=colour, linewidths=0.5))


fig, axes = plt.subplots(1, 7, figsize=(183 * MM, 62 * MM), gridspec_kw={"wspace": 0.06})
ax = axes[0]
img = np.asarray(Image.open(DATA / "visium_tables" / f"{RUN}_hires.png").convert("RGB")) \
    if (DATA / "visium_tables" / f"{RUN}_hires.png").exists() else None
sfk = "tissue_hires_scalef" if img is not None else "tissue_lowres_scalef"
if img is None:
    img = np.asarray(Image.open(DATA / "visium_tables" / f"{RUN}_lowres.png").convert("RGB"))
f = json.loads((DATA / "visium_tables" / f"{RUN}_scalefactors.json").read_text())[sfk]
ax.imshow(img)
ax.set_xlim(lo[0] * f, hi[0] * f)
ax.set_ylim(hi[1] * f, lo[1] * f)
ax.axis("off")
outline(ax, "black", f)
ax.set_title("H&E", fontweight="bold")
for j, k in [(1, 2), (2, 3)]:
    s0 = sub[f"k{k}_s0"]
    nm = {g: name(k, g) for g in s0.unique()}
    colours = [SUBCOL[nm[s0[b]]] if b in s0.index else (0.85, 0.85, 0.85, 1) for b in idx]
    draw(axes[j], colours)
    outline(axes[j])
    axes[j].set_title(f"Sub-clusters\nK = {k}", fontweight="bold")
for j, (g, role) in enumerate(GENES, start=3):
    v = gene(g)
    norm = mpl.colors.Normalize(float(np.quantile(v[win], 0.01)), float(np.quantile(v[win], 0.99)), clip=True)
    draw(axes[j], mpl.colormaps["magma"](norm(v)))
    outline(axes[j], "white")
    axes[j].set_title(f"{g}\n{role}", fontweight="bold")
    cax = axes[j].inset_axes([0.1, -0.06, 0.8, 0.03])
    cb = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap="magma"), cax=cax, orientation="horizontal")
    cb.ax.tick_params(labelsize=5, length=1.5, width=0.4)
    cb.set_label("log1p(CP10k)", fontsize=5.5)
    cb.outline.set_linewidth(0.4)
axes[1].legend(handles=[Patch(color=c, label=n) for n, c in SUBCOL.items()] + [Patch(color=(0.85, 0.85, 0.85), label="outside vessel domain")],
               loc="upper left", bbox_to_anchor=(0.0, -0.02), ncol=1, frameon=False, fontsize=5.5, handlelength=1)
axes[0].text(-0.06, 0.5, "P9 IZ rep2\n(ACH0013)", transform=axes[0].transAxes, rotation=90, ha="right",
             va="center", fontweight="bold", color="#D55E00", fontsize=7)
fig.savefig(OUT / "fig_p9rep2_vessel_subclusters.png", dpi=300, bbox_inches="tight")
fig.savefig(OUT / "fig_p9rep2_vessel_subclusters.pdf", dpi=300, bbox_inches="tight")
print(f"zoom: {win.sum()} spots, main band piece {len(main)} of {band.sum()} band spots")
