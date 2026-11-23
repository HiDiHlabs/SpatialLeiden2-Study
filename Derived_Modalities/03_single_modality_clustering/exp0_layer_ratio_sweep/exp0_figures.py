"""Figures for experiment 0 (layer_ratio sweep).

Fig. 1  coherence gap (observed - shuffled) vs layer_ratio, mean and seed range, one line per slide
Fig. 2  ARI to the ratio-0 partition (same seed) vs layer_ratio
Fig. 3  seed stability (pairwise ARI among the 5 seeds) vs layer_ratio
Fig. 5  how hard the search works: spatial resolution found, and the largest cluster's share
Fig. 4  per slide: H&E and shipped cell2location maps beside the domains at selected ratios (seed 0);
        cluster colours are matched to the ratio-0 partition so the same domain keeps its colour

Usage (from Derived_Modalities/):
    python 03_single_modality_clustering/exp0_layer_ratio_sweep/exp0_figures.py data/kuppe_2022 results/kuppe_2022/exp0
"""

import json
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import EllipseCollection
from PIL import Image
from scipy.optimize import linear_sum_assignment
from scipy.spatial import cKDTree

DATA, RES = Path(sys.argv[1]), Path(sys.argv[2])
FIG = RES / "figures"
FIG.mkdir(exist_ok=True)
MM = 1 / 25.4
mpl.rcParams.update(
    {
        "font.family": "Arial",
        "font.size": 6,
        "axes.titlesize": 6.5,
        "axes.labelsize": 6,
        "xtick.labelsize": 5.5,
        "ytick.labelsize": 5.5,
        "legend.fontsize": 5.5,
        "axes.linewidth": 0.5,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "pdf.fonttype": 42,
    }
)

df = pd.read_csv(RES / "exp0_sweep.csv")
KS = df.groupby("slide").target_k.first().to_dict()
K = int(df.target_k.iloc[0])
SLIDE_COL = {
    "ACH0013": "#D55E00",
    "ACH0012": "#E69F00",
    "10X001": "#0072B2",
}  # IZ deep, IZ typical, control
LABEL = dict(zip(df.slide, df.label))
SLIDES = [s for s in ["ACH0013", "ACH0012", "10X001"] if s in set(df.slide)]


def sweep_plot(ax, col, ylabel, lo=None, hi=None):
    for s in SLIDES:
        d = df[df.slide == s]
        g = d.groupby("layer_ratio")
        m = g[col].mean()
        a, b = (g[lo].mean(), g[hi].mean()) if lo else (g[col].min(), g[col].max())
        ax.fill_between(m.index, a, b, color=SLIDE_COL[s], alpha=0.2, lw=0)
        ax.plot(
            m.index, m.values, "-o", color=SLIDE_COL[s], ms=2.5, lw=1, label=LABEL[s]
        )
    ax.set_xlabel("layer_ratio")
    ax.set_ylabel(ylabel)


fig, axes = plt.subplots(1, 3, figsize=(183 * MM, 55 * MM))
sweep_plot(
    axes[0],
    "coherence_gap",
    "Coherence above shuffled\n(share of neighbour pairs, same label)",
)
sweep_plot(axes[1], "ARI_to_ratio0", "ARI to ratio-0 partition (same seed)")
axes[1].set_ylim(0, 1.02)
sweep_plot(
    axes[2],
    "seed_ARI_mean",
    "Seed stability\n(mean pairwise ARI, 5 seeds; band = min-max)",
    "seed_ARI_min",
    "seed_ARI_max",
)
axes[2].set_ylim(0, 1.02)
for ax, t in zip(
    axes,
    [
        "1  Spatial coherence",
        "2  Distance from the no-space baseline",
        "3  Reproducibility",
    ],
):
    ax.set_title(t, loc="left", fontweight="bold")
axes[0].legend(frameon=False, loc="upper left")
ktxt = (
    f"K = {K} clusters"
    if len(set(KS.values())) == 1
    else "K per slide: " + ", ".join(f"{LABEL[s]} {KS[s]}" for s in SLIDES)
)
fig.text(
    0.5,
    -0.02,
    f"{ktxt}, held at every layer_ratio. Line: mean over 5 seeds; band: range over seeds.",
    ha="center",
    fontsize=5.5,
)
fig.tight_layout()
fig.savefig(FIG / "fig1-3_sweep.pdf", bbox_inches="tight")
fig.savefig(FIG / "fig1-3_sweep.png", dpi=300, bbox_inches="tight")
plt.close(fig)

fig, axes = plt.subplots(1, 3, figsize=(183 * MM, 50 * MM))
sweep_plot(axes[0], "res_spatial", "Spatial resolution found")
sweep_plot(axes[1], "largest_cluster_frac", "Largest cluster (share of spots)")
sweep_plot(
    axes[2],
    "ARI_celltype_niche",
    "ARI to authors' cell-type niches\n(context, not accuracy)",
)
axes[0].legend(frameon=False)
for ax, t in zip(
    axes, ["5a  Search effort", "5b  Cluster sizes", "5c  Published niches"]
):
    ax.set_title(t, loc="left", fontweight="bold")
fig.tight_layout()
fig.savefig(FIG / "fig5_resolution_sizes_niches.pdf", bbox_inches="tight")
fig.savefig(FIG / "fig5_resolution_sizes_niches.png", dpi=300, bbox_inches="tight")
plt.close(fig)

# ---------------------------------------------------------------- Fig. 4: domain maps
SHOW = [0, 0.5, 1.0, 2.0]
C2L = [
    ("prop_CM", "Cardiomyocytes"),
    ("prop_Fib", "Fibroblasts"),
    ("prop_Myeloid", "Myeloid"),
    ("prop_vSMCs", "vSMCs"),
]
PAL = [
    c for i, c in enumerate(mpl.colormaps["tab20"].colors) if i not in (14, 15)
]  # drop greys
PAL += list(mpl.colormaps["Dark2"].colors)


def draw(ax, x, y, colours, title):
    pitch = np.median(cKDTree(np.c_[x, y]).query(np.c_[x, y], k=2)[0][:, 1])
    d = np.full(len(x), 1.05 * pitch)
    coll = EllipseCollection(
        d,
        d,
        np.zeros(len(x)),
        units="xy",
        offsets=np.c_[x, y],
        offset_transform=ax.transData,
        facecolors=colours,
        linewidths=0,
    )
    coll.set_rasterized(True)
    ax.add_collection(coll)
    ax.set_xlim(x.min() - pitch, x.max() + pitch)
    ax.set_ylim(y.max() + pitch, y.min() - pitch)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title(title)


def match_to(ref, lab):
    """Relabel `lab` so each cluster takes the id of the ratio-0 cluster it overlaps most."""
    ct = pd.crosstab(lab, ref)
    r, c = linear_sum_assignment(-ct.to_numpy())
    mapping = {ct.index[i]: ct.columns[j] for i, j in zip(r, c)}
    extra = [k for k in ct.index if k not in mapping]
    nxt = max(ref.astype(int)) + 1
    for k in extra:
        mapping[k] = str(nxt)
        nxt += 1
    return lab.map(mapping)


for s in SLIDES:
    lab = pd.read_csv(RES / f"labels_{s}.csv.gz", index_col=0, dtype=str)
    spots = (
        pd.read_csv(DATA / "visium_tables" / f"{s}_spots.csv")
        .set_index("barcode")
        .loc[lab.index]
    )
    x, y = spots.pxl_col_fullres.to_numpy(float), spots.pxl_row_fullres.to_numpy(float)
    fig, axes = plt.subplots(
        2, 5, figsize=(183 * MM, 80 * MM), gridspec_kw={"wspace": 0.05, "hspace": 0.25}
    )
    # row 1: H&E + cell2location
    img = np.asarray(
        Image.open(DATA / "visium_tables" / f"{s}_lowres.png").convert("RGB")
    )
    f = json.loads((DATA / "visium_tables" / f"{s}_scalefactors.json").read_text())[
        "tissue_lowres_scalef"
    ]
    ax = axes[0, 0]
    ax.imshow(img)
    ax.set_xlim(x.min() * f - 8, x.max() * f + 8)
    ax.set_ylim(y.max() * f + 8, y.min() * f - 8)
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_title("H&E")
    for ax, (col, name) in zip(axes[0, 1:], C2L):
        v = spots[col].to_numpy()
        norm = mpl.colors.Normalize(0, np.quantile(v, 0.99))
        draw(ax, x, y, mpl.colormaps["magma"](norm(v)), name)
    # row 2: domains
    ref = lab["r0_s0"]
    for ax, ratio in zip(axes[1, :4], SHOW):
        key = f"r{ratio}_s0" if ratio != 0 else "r0_s0"
        if key not in lab:
            key = f"r{float(ratio)}_s0"
        m = match_to(ref, lab[key]) if ratio != 0 else ref
        draw(ax, x, y, [PAL[int(v) % len(PAL)] for v in m], f"layer_ratio {ratio}")
    ct = spots.celltype_niche.str.removeprefix("ctniche_").astype(int)
    niche_col = {
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
    draw(axes[1, 4], x, y, ct.map(niche_col).to_list(), "Authors' niches")
    fig.suptitle(
        f"{LABEL[s]}. Top: H&E and shipped cell2location proportions. Bottom: SpatialLeiden domains at K = {KS[s]} (seed 0; colours matched to the ratio-0 partition) and the authors' cell-type niches",
        fontsize=7,
        y=0.99,
    )
    fig.savefig(FIG / f"fig4_domains_{s}.pdf", bbox_inches="tight", dpi=300)
    fig.savefig(FIG / f"fig4_domains_{s}.png", dpi=250, bbox_inches="tight")
    plt.close(fig)
print("figures written to", FIG)
