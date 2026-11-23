"""Multipanel: SpatialLeiden on RNA alone, TF activity alone, and RNA + TF, on P1 and P9.

Columns: RNA only and TF activity only (spatial weight 0, the experiment-0 no-space
baselines), then three runs
with the spatial graph. Settings are identical so only the input layers differ:
    RNA only           experiment 0, K = authors' cell-type niches, layer_ratio 0, seed 0 (label r0_s0)
    RNA + space        experiment 0, K = authors' cell-type niches, layer_ratio 0.75, seed 0
                       (results/kuppe_2022/exp0_k_celltype_niche, label r0.75_s0)
    TF + space         same, TF graph instead of RNA (results/kuppe_2022/exp0tf_k_celltype_niche)
    RNA + TF + space   multimodal prototype: RNA 1, TF 0.5, spatial 0.75, seed 0
                       (results/kuppe_2022/mm_proto_rna_tf, label sp0.75_tf0.5_s0)
Colours: each map's domains are matched (Hungarian, by overlap) to the RNA + space domains,
so equal colours mean best-matching domains. The authors' niches keep their own colours.

Usage (from Derived_Modalities/):
    python 04_multimodal_weight_grid/proto_rna_tf/fig_rna_tf_combo.py data/kuppe_2022 results/kuppe_2022
Output: results/kuppe_2022/mm_proto_rna_tf/figures/fig_rna_tf_combo.png/.pdf and
        fig_rna_tf_combo_settings.csv (the exact runs shown, with K reached and coherence)
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
OUT = RES / "mm_proto_rna_tf" / "figures"
OUT.mkdir(parents=True, exist_ok=True)
MM = 1 / 25.4
mpl.rcParams.update({"font.family": "Arial", "font.size": 6, "axes.titlesize": 6.5, "pdf.fonttype": 42})
SLIDES = {"10X001": "P1 control\n(10X001)", "ACH0012": "P9 IZ\n(ACH0012)", "ACH0013": "P9 IZ rep2\n(ACH0013)"}
SLIDE_COL = {"10X001": "#0072B2", "ACH0013": "#D55E00", "ACH0012": "#E69F00"}
NICHE_COL = {f"ctniche_{i}": c for i, c in enumerate(
    ["#D51F26", "#272E6A", "#208A42", "#89288F", "#F47D2B", "#FEE500", "#8A9FD1", "#C06CAB", "#D8A767"], start=1)}
PAL = [c for i, c in enumerate(mpl.colormaps["tab20"].colors) if i not in (14, 15)] + list(mpl.colormaps["Dark2"].colors)
# Marker panel: P9 is an ischaemic-zone section (GT/IZ, 2 days after onset); the paper reports
# "an increased abundance of a macrophage population defined by SPP1 expression in the
# ischaemic sample group" (Kuppe et al. 2022, Fig. 6j-k, Extended Data Fig. 11f). SPP1 is shown
# as log-normalised expression (normalize_total 1e4, log1p) on one colour scale for all slides.
# NPPB: the paper's stressed-cardiomyocyte marker ("NPPB showed a distinct upregulation and increased
# chromatin accessibility in vCM3"; reported in the border zone after infarction).
MARKERS = ["SPP1", "NPPB"]
H5AD = {"10X001": "Visium_control_P1.h5ad", "ACH0013": "Visium_GT_IZ_P9_rep2.h5ad", "ACH0012": "Visium_GT_IZ_P9.h5ad"}
RUNS = [("RNA only", "exp0_k_celltype_niche", "r0_s0"),
        ("RNA + space", "exp0_k_celltype_niche", "r0.75_s0"),
        ("TF activity only", "exp0tf_k_celltype_niche", "r0_s0"),
        ("TF activity + space", "exp0tf_k_celltype_niche", "r0.75_s0"),
        ("RNA + TF + space", "mm_proto_rna_tf", "sp0.75_tf0.5_s0")]


def draw(ax, xy, colours):
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


def match_to(ref, lab):
    tab = pd.crosstab(pd.Series(lab), pd.Series(ref))
    r, c = linear_sum_assignment(-tab.to_numpy())
    mapping = {tab.index[i]: tab.columns[j] for i, j in zip(r, c)}
    for n, k in enumerate(k for k in tab.index if k not in mapping):
        mapping[k] = f"x{n}"
    return pd.Series(lab).map(mapping).to_numpy()


def coherence(lab, i, j):
    lab = np.asarray(lab)
    return float((lab[i] == lab[j]).mean())


NM = len(MARKERS)
WRAP = {"TF activity only": "TF activity\nonly", "RNA + space": "RNA\n+ space", "TF activity + space": "TF activity\n+ space", "RNA + TF + space": "RNA + TF\n+ space"}
cols = (["H&E"] + [WRAP.get(r[0], r[0]) for r in RUNS] + ["Authors'\nniches"]
        + [f"{m}\nexpression" for m in MARKERS])
J_RUN0 = 1                      # first clustering column
J_AUTH = 1 + len(RUNS)          # authors' niches
J_MK0 = J_AUTH + 1              # first marker column
fig, axes = plt.subplots(len(SLIDES), len(cols), figsize=(183 * MM, 32 * MM * len(SLIDES)),
                         gridspec_kw={"wspace": 0.04, "hspace": 0.08})
rows = []
expr = {m: {} for m in MARKERS}
for run in SLIDES:  # marker expression, log-normalised, in the h5ad barcode order
    import anndata as ad
    h5 = ad.read_h5ad(DATA / "processed1" / H5AD[run])
    counts = h5.raw.X
    tot = np.asarray(counts.sum(1)).ravel()
    for mk in MARKERS:
        g = list(h5.raw.var_names).index(mk)
        expr[mk][run] = pd.Series(np.log1p(np.asarray(counts[:, g].todense()).ravel() / tot * 1e4),
                                  index=h5.obs_names)
vmax = {mk: float(np.quantile(np.concatenate([v.to_numpy() for v in expr[mk].values()]), 0.99)) for mk in MARKERS}
for i, (run, label) in enumerate(SLIDES.items()):
    spots = pd.read_csv(DATA / "visium_tables" / f"{run}_spots.csv").set_index("barcode")
    labs = {}
    for name, folder, key in RUNS:
        t = pd.read_csv(RES / folder / f"labels_{run}.csv.gz", index_col=0, dtype=str)
        labs[name] = t[key]
    idx = labs["RNA + space"].index  # every run is checked against this order below
    for name in labs:
        assert list(labs[name].index) == list(idx), "spot order differs between runs"
    sp = spots.loc[idx]
    xy = sp[["pxl_col_fullres", "pxl_row_fullres"]].to_numpy(float)
    # neighbour pairs on the Visium grid (array coordinates), for the coherence numbers in the CSV
    pos = {(r, c): n for n, (r, c) in enumerate(zip(sp.array_row, sp.array_col))}
    ei, ej = zip(*[(n, pos[(r + dr, c + dc)]) for (r, c), n in pos.items()
                   for dr, dc in [(0, 2), (0, -2), (1, 1), (1, -1), (-1, 1), (-1, -1)] if (r + dr, c + dc) in pos])
    ei, ej = np.array(ei), np.array(ej)
    # H&E
    ax = axes[i, 0]
    img = np.asarray(Image.open(DATA / "visium_tables" / f"{run}_lowres.png").convert("RGB"))
    f = json.loads((DATA / "visium_tables" / f"{run}_scalefactors.json").read_text())["tissue_lowres_scalef"]
    ax.imshow(img)
    ax.set_xlim(xy[:, 0].min() * f - 8, xy[:, 0].max() * f + 8)
    ax.set_ylim(xy[:, 1].max() * f + 8, xy[:, 1].min() * f - 8)
    ax.axis("off")
    ax.text(-0.06, 0.5, label, transform=ax.transAxes, rotation=90, ha="right", va="center",
            fontweight="bold", color=SLIDE_COL[run], fontsize=7)
    ref = labs["RNA + space"].to_numpy()
    cmap = {c: PAL[k % len(PAL)] for k, c in enumerate(sorted(set(ref), key=int))}
    for jm, mk in enumerate(MARKERS, start=J_MK0):
        norm = mpl.colors.Normalize(0, vmax[mk], clip=True)
        draw(axes[i, jm], xy, mpl.colormaps["magma"](norm(expr[mk][run].loc[idx].to_numpy())))
        rows.append({"slide": run, "run": f"{mk} expression", "key": "log1p(CP10k)",
                     "n_domains": np.nan, "coherence_edge_weighted": np.nan,
                     "median": round(float(expr[mk][run].median()), 3),
                     "share_spots_detected": round(float((expr[mk][run] > 0).mean()), 3)})
    for j, (name, _, key) in enumerate(RUNS, start=J_RUN0):
        lab = labs[name].to_numpy()
        m = ref if name == "RNA + space" else match_to(ref, lab)
        draw(axes[i, j], xy, [cmap.get(v, "#777777") for v in m])
        rows.append({"slide": run, "run": name, "key": key, "n_domains": len(set(lab)),
                     "coherence_edge_weighted": round(coherence(lab, ei, ej), 3)})
    draw(axes[i, J_AUTH], xy, [NICHE_COL[v] for v in sp.celltype_niche])
for j, c in enumerate(cols):
    axes[0, j].set_title(c, fontweight="bold")
# colour scale for each marker column, below its last panel
for jm, mk in enumerate(MARKERS, start=J_MK0):
    cax = axes[-1, jm].inset_axes([0.1, -0.12, 0.8, 0.05])
    cb = fig.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(0, vmax[mk]), cmap="magma"), cax=cax,
                      orientation="horizontal")
    cb.ax.tick_params(labelsize=5, length=1.5, width=0.4)
    cb.set_label("log1p(CP10k)", fontsize=5.5)
    cb.outline.set_linewidth(0.4)
# settings are in the module docstring and in fig_rna_tf_combo_settings.csv (no footnote on the figure)
fig.savefig(OUT / "fig_rna_tf_combo.png", dpi=250, bbox_inches="tight")
fig.savefig(OUT / "fig_rna_tf_combo.pdf", dpi=300, bbox_inches="tight")
pd.DataFrame(rows).to_csv(OUT / "fig_rna_tf_combo_settings.csv", index=False)
print(pd.DataFrame(rows).to_string(index=False))
