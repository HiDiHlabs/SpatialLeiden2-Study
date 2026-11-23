"""Figures for the RNA + TF + space multimodal prototype (mm_proto.py).

fig_mm_curves.png/.pdf  coherence gap, ARI to RNA + space (TF weight 0), seed stability and
                        ARI to the authors' cell-type niches against the TF weight; one colour
                        per slide, solid = spatial weight 0.75, dashed = 0.5. Means over runs
                        that reached K; the number of such seeds is printed under each point
                        when it is below 5.
fig_mm_maps_<run>.png   domains at TF weight 0 / 0.25 / 0.5 / 1 / 2 (spatial weight 0.75,
                        seed 0; runs that missed K are marked), beside the authors' niches.

Usage (from Derived_Modalities/):
    python 04_multimodal_weight_grid/proto_rna_tf/mm_proto_figures.py data/kuppe_2022 results/kuppe_2022/mm_proto_rna_tf
"""

import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import EllipseCollection
from scipy.optimize import linear_sum_assignment
from scipy.spatial import cKDTree

DATA, RES = Path(sys.argv[1]), Path(sys.argv[2])
FIG = RES / "figures"
FIG.mkdir(exist_ok=True)
MM = 1 / 25.4
mpl.rcParams.update({"font.family": "Arial", "font.size": 6, "axes.titlesize": 6.5, "axes.labelsize": 6,
                     "xtick.labelsize": 5.5, "ytick.labelsize": 5.5, "legend.fontsize": 5.5,
                     "axes.linewidth": 0.5, "axes.spines.top": False, "axes.spines.right": False,
                     "pdf.fonttype": 42})
SLIDES = {"ACH0013": "P9 IZ rep2 (ACH0013)", "ACH0012": "P9 IZ (ACH0012)", "10X001": "P1 control (10X001)"}
SLIDE_COL = {"ACH0013": "#D55E00", "ACH0012": "#E69F00", "10X001": "#0072B2"}
NICHE_COL = {f"ctniche_{i}": c for i, c in enumerate(
    ["#D51F26", "#272E6A", "#208A42", "#89288F", "#F47D2B", "#FEE500", "#8A9FD1", "#C06CAB", "#D8A767"], start=1)}
PAL = [c for i, c in enumerate(mpl.colormaps["tab20"].colors) if i not in (14, 15)] + list(mpl.colormaps["Dark2"].colors)

d = pd.read_csv(RES / "mm_proto_sweep.csv")
ok = d[d.reached_k]
K = d.groupby("slide").target_k.first().to_dict()

fig, axes = plt.subplots(1, 4, figsize=(183 * MM, 50 * MM))
for ax, col, ylab in [(axes[0], "coherence_gap", "Coherence above shuffled"),
                      (axes[1], "ARI_to_rna_space", "ARI to RNA + space\n(TF weight 0)"),
                      (axes[2], "seed_ARI_mean", "Seed stability (lower bound)"),
                      (axes[3], "ARI_celltype_niche", "ARI to authors' cell-type niches")]:
    for s in SLIDES:
        for w_sp, ls in [(0.75, "-"), (0.5, "--")]:
            g = ok[(ok.slide == s) & (ok.spatial_weight == w_sp)].groupby("tf_weight")
            m, n = g[col].mean(), g.size()
            ax.plot(m.index, m.values, ls, marker="o", ms=2.2, lw=1, color=SLIDE_COL[s],
                    label=f"{SLIDES[s].split(' (')[0]}, space {w_sp}")
            if col == "coherence_gap":
                for x, y, k in zip(m.index, m.values, n.values):
                    if k < 5:
                        ax.annotate(f"n={k}", (x, y), textcoords="offset points", xytext=(2, -6), fontsize=4.5,
                                    color=SLIDE_COL[s])
    ax.set_xlabel("TF-layer weight (RNA weight = 1)")
    ax.set_ylabel(ylab)
    if col != "coherence_gap":
        ax.set_ylim(0, 1.02)
axes[0].legend(frameon=False, fontsize=4.8, loc="lower left")
fig.text(0.5, -0.06, "SpatialLeiden multimodal: RNA graph + TF-activity graph (decoupler ULM, CollecTRI) + spatial grid. "
         f"K = authors' cell-type niche count ({', '.join(f'{SLIDES[s].split(' (')[0]} {K[s]}' for s in SLIDES)}). "
         "Means over seeds that reached K (n shown when < 5): at high TF weight the spatial-resolution search alone "
         "cannot restore K.", ha="center", fontsize=5.3, wrap=True)
fig.tight_layout()
fig.savefig(FIG / "fig_mm_curves.png", dpi=300, bbox_inches="tight")
fig.savefig(FIG / "fig_mm_curves.pdf", bbox_inches="tight")
plt.close(fig)


def draw(ax, xy, colours, title):
    pitch = np.median(cKDTree(xy).query(xy, k=2)[0][:, 1])
    dd = np.full(len(xy), 1.05 * pitch)
    coll = EllipseCollection(dd, dd, np.zeros(len(xy)), units="xy", offsets=xy, offset_transform=ax.transData,
                             facecolors=colours, linewidths=0)
    coll.set_rasterized(True)
    ax.add_collection(coll)
    ax.set_xlim(xy[:, 0].min() - pitch, xy[:, 0].max() + pitch)
    ax.set_ylim(xy[:, 1].max() + pitch, xy[:, 1].min() - pitch)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(title)


def match_to(ref, lab):
    tab = pd.crosstab(pd.Series(lab), pd.Series(ref))
    r, c = linear_sum_assignment(-tab.to_numpy())
    mapping = {tab.index[i]: tab.columns[j] for i, j in zip(r, c)}
    for n, k in enumerate(k for k in tab.index if k not in mapping):
        mapping[k] = f"x{n}"
    return pd.Series(lab).map(mapping).to_numpy()


TFW = [0, 0.25, 0.5, 1.0, 2.0]
for s in SLIDES:
    lab = pd.read_csv(RES / f"labels_{s}.csv.gz", index_col=0, dtype=str)
    sp = pd.read_csv(DATA / "visium_tables" / f"{s}_spots.csv").set_index("barcode").loc[lab.index]
    xy = sp[["pxl_col_fullres", "pxl_row_fullres"]].to_numpy(float)
    fig, axes = plt.subplots(1, len(TFW) + 1, figsize=(183 * MM, 40 * MM), gridspec_kw={"wspace": 0.04})
    ref = lab["sp0.75_tf0_s0"].to_numpy()
    cats = sorted(set(ref), key=int)
    col = {c: PAL[i % len(PAL)] for i, c in enumerate(cats)}
    for ax, w in zip(axes, TFW):
        key = f"sp0.75_tf{w}_s0" if f"sp0.75_tf{w}_s0" in lab else f"sp0.75_tf{int(w)}_s0"
        m = ref if w == 0 else match_to(ref, lab[key].to_numpy())
        nk = len(set(lab[key]))
        flag = "" if nk == K[s] else f"\n{nk} domains (missed K)"
        draw(ax, xy, [col.get(v, "#777777") for v in m], f"TF weight {w:g}{flag}")
    draw(axes[-1], xy, [NICHE_COL[v] for v in sp.celltype_niche], "Authors' cell-type niches")
    fig.suptitle(f"{SLIDES[s]}: RNA (weight 1) + TF activity + space (weight 0.75), K = {K[s]}, seed 0; "
                 "colours matched to TF weight 0 (= RNA + space)", fontsize=6.5, y=1.02)
    fig.savefig(FIG / f"fig_mm_maps_{s}.png", dpi=250, bbox_inches="tight")
    fig.savefig(FIG / f"fig_mm_maps_{s}.pdf", dpi=300, bbox_inches="tight")
    plt.close(fig)
print("figures written to", FIG)
