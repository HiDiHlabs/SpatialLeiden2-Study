"""One multipanel comparing three joint clusterings (all 27 slides, seed 0, spatial weight 0.75):
    RNA + space          results/kuppe_2022/exp2_joint_harmony        r0.75_s0
    TF + space           results/kuppe_2022/exp2_joint_harmony_tf     r0.75_s0 (TF_FOLDER below; while the
                         sweep was finishing, the same run recomputed from its logged resolutions)
    RNA + TF + space     results/kuppe_2022/exp3_joint_rna_tf         r0.75_s0 (TF weight 0.5)
Domains are numbered D1.. by size within each clustering (as in domain_enrichment_by_zone.py).
Colours: each RNA domain has its own colour; each RNA + TF domain takes the colour of the RNA
domain most of its spots come from (same rule for TF + space). If a larger domain of the same
clustering already took that colour, the domain is NEW (marked "new" / "*") and gets its own
colour (RNA + TF + space is coloured first, so its new domain is pink).
Labels: the two cell types with the highest median cell2location proportion in the domain.
Top: maps on four slides (P1 control, BZ P3, P9 IZ, P9 IZ rep2): H&E | RNA + space | TF + space |
RNA + TF + space | NPPB (stressed cardiomyocytes) with the new RNA + TF domain(s) outlined.
Bottom: a where each RNA + TF domain comes from (share of its spots per RNA domain);
b-d mean share of each domain per zone (slide = unit; from domain_enrichment_by_zone.py).

Mode "ccc" (third argument): the same figure with CCC (our LIANA+ layer) instead of TF, at the CCC
spatial weight 0.25: CCC + space (exp2_joint_harmony_ccc_ours, r0.25_s0) and RNA + CCC + space
(exp3_joint_rna_ccc_ours, r0.25_s0; recomputed from its logged resolutions while the sweep was
finishing). RNA + space stays at 0.75.

Usage (from Derived_Modalities/; run domain_enrichment_by_zone.py for all three clusterings first):
    python 03_single_modality_clustering/exp2_joint_harmony/fig_joint_rna_vs_rnatf.py data/kuppe_2022 results/kuppe_2022 [tf|ccc]
Output: results/kuppe_2022/joint_domain_enrichment/fig_joint_rna_vs_rna{tf,ccc}.png/.pdf
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
OUT = RES / "joint_domain_enrichment"
COL = "r0.75_s0"
MM = 1 / 25.4
mpl.rcParams.update({"font.family": "Arial", "font.size": 6, "axes.titlesize": 6.5, "pdf.fonttype": 42})
SLIDES = {"10X001": "P1 control", "10X0026": "BZ (P3)", "ACH0012": "P9 IZ", "ACH0013": "P9 IZ rep2"}
ZONES = ["CTRL", "RZ", "BZ", "IZ", "FZ"]
RNA_COLS = ["#1f4e9c", "#8fb8e0", "#c49c6b", "#d62728", "#ff9f40", "#8a5fbf", "#7f7f7f", "#d4c62a", "#cfcfcf"]
NEW_COLS = ["#e43fb1", "#17becf", "#2ca02c", "#9c6b30", "#4b0082", "#00b894", "#ff6f91", "#6b8e23"]
MODE = sys.argv[3] if len(sys.argv) > 3 else "tf"


def first(*folders):
    return next((f for f in folders if (RES / f / "exp2_labels.csv.gz").exists()), folders[-1])


if MODE == "tf":
    L, WL = "TF", "0.75"
    SINGLE = (first("exp2_joint_harmony_tf", "exp2_joint_harmony_tf_r075_fromlog"), "r0.75_s0")
    COMBO = ("exp3_joint_rna_tf", "r0.75_s0")
else:
    L, WL = "CCC", "0.25"
    SINGLE = ("exp2_joint_harmony_ccc_ours", "r0.25_s0")
    COMBO = (first("exp3_joint_rna_ccc_ours", "exp3_joint_rna_ccc_ours_r025_fromlog"), "r0.25_s0")
TF_FOLDER = SINGLE[0]
FIG = f"fig_joint_rna_vs_rna{MODE}"


def load(folder, col=COL):
    lab = pd.read_csv(RES / folder / "exp2_labels.csv.gz", index_col=0, dtype=str)[col]
    order = lab.value_counts().index
    return lab.map({d: f"D{i + 1}" for i, d in enumerate(order)})


rna, tf, tfo = load("exp2_joint_harmony"), load(*COMBO), load(*SINGLE)
tf, tfo = tf.loc[rna.index], tfo.loc[rna.index]
spots = []
for r in pd.unique(rna.index.str.split("_", n=1).str[0]):
    t = pd.read_csv(DATA / "visium_tables" / f"{r}_spots.csv").set_index("barcode")
    t.index = r + "_" + t.index
    spots.append(t)
spots = pd.concat(spots).loc[rna.index]
props = [c for c in spots.columns if c.startswith("prop_")]


def describe(lab):
    med = spots.groupby(lab.to_numpy())[props].median()
    return {d: "+".join(c.replace("prop_", "") for c in med.loc[d].sort_values(ascending=False).index[:2]) for d in med.index}


rna_doms = [f"D{i}" for i in range(1, rna.nunique() + 1)]
tf_doms = [f"D{i}" for i in range(1, tf.nunique() + 1)]
rna_col = dict(zip(rna_doms, RNA_COLS))
new_pool = {}


def colour(lab):
    doms = [f"D{i}" for i in range(1, lab.nunique() + 1)]
    src_ = pd.crosstab(lab, rna, normalize="index").reindex(index=doms, columns=rna_doms, fill_value=0)
    col, used, new_ = {}, set(), []
    for d in doms:  # D1 is the largest, so larger domains claim a colour first
        major = src_.loc[d].idxmax()
        if major not in used:
            col[d] = rna_col[major]
            used.add(major)
        else:
            col[d] = NEW_COLS[len(new_pool) % len(NEW_COLS)]  # every new domain its own colour
            new_pool[(id(lab), d)] = col[d]
            new_.append(d)
    return doms, src_, col, new_


tf_doms, src, tf_col, new = colour(tf)
tfo_doms, _src_o, tfo_col, new_o = colour(tfo)
desc_r, desc_t, desc_o = describe(rna), describe(tf), describe(tfo)

# expression of NPPB per slide shown
nppb = {}
tables = {r: set(spots.index[spots.index.str.startswith(r + "_")].str.split("_", n=1).str[1]) for r in SLIDES}
for f in sorted((DATA / "processed1").glob("Visium_*.h5ad")):
    names = set(ad.read_h5ad(f, backed="r").obs_names)
    run = [r for r, b in tables.items() if b == names]
    if not run:
        continue
    h = ad.read_h5ad(f)
    tot = np.asarray(h.raw.X.sum(1)).ravel()
    g = list(h.raw.var_names).index("NPPB")
    nppb[run[0]] = pd.Series(np.log1p(np.asarray(h.raw.X[:, g].todense()).ravel() / tot * 1e4), index=run[0] + "_" + h.obs_names)
vmax = float(np.quantile(np.concatenate([v.to_numpy() for v in nppb.values()]), 0.99))


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


def outline(sp, mask):
    """Hexagon edges between a masked spot and any grid neighbour outside the mask."""
    xy = sp[["pxl_col_fullres", "pxl_row_fullres"]].to_numpy(float)
    grid = np.c_[sp.array_col, sp.array_row].astype(float)
    coef, *_ = np.linalg.lstsq(np.c_[grid, np.ones(len(grid))], xy, rcond=None)
    pos = {(int(r), int(c)): n for n, (r, c) in enumerate(zip(sp.array_row, sp.array_col))}
    segs = []
    for n in np.flatnonzero(mask):
        r, c = int(sp.array_row.iloc[n]), int(sp.array_col.iloc[n])
        for dr, dc in [(0, 2), (0, -2), (1, 1), (1, -1), (-1, 1), (-1, -1)]:
            m = pos.get((r + dr, c + dc))
            if m is not None and mask[m]:
                continue
            v = np.array([dc, dr]) @ coef[:2]
            mid = xy[n] + v / 2
            perp = np.array([-v[1], v[0]]) / (2 * np.sqrt(3))
            segs.append([mid - perp, mid + perp])
    return np.array(segs) if segs else np.zeros((0, 2, 2))


fig = plt.figure(figsize=(183 * MM, 245 * MM))
gs = fig.add_gridspec(3, 1, height_ratios=[2.6, 0.95, 0.95], hspace=0.3)
top = gs[0].subgridspec(len(SLIDES), 5, wspace=0.04, hspace=0.06)
for i, (run, label) in enumerate(SLIDES.items()):
    idx = spots.index[spots.index.str.startswith(run + "_")]
    sp = spots.loc[idx]
    xy = sp[["pxl_col_fullres", "pxl_row_fullres"]].to_numpy(float)
    axs = [fig.add_subplot(top[i, j]) for j in range(5)]
    img = np.asarray(Image.open(DATA / "visium_tables" / f"{run}_lowres.png").convert("RGB"))
    f = json.loads((DATA / "visium_tables" / f"{run}_scalefactors.json").read_text())["tissue_lowres_scalef"]
    axs[0].imshow(img)
    axs[0].set_xlim(xy[:, 0].min() * f - 8, xy[:, 0].max() * f + 8)
    axs[0].set_ylim(xy[:, 1].max() * f + 8, xy[:, 1].min() * f - 8)
    axs[0].axis("off")
    axs[0].text(-0.06, 0.5, label, transform=axs[0].transAxes, rotation=90, ha="right", va="center",
                fontweight="bold", fontsize=7)
    draw(axs[1], xy, [rna_col[d] for d in rna.loc[idx]])
    draw(axs[2], xy, [tfo_col[d] for d in tfo.loc[idx]])
    draw(axs[3], xy, [tf_col[d] for d in tf.loc[idx]])
    draw(axs[4], xy, mpl.colormaps["magma"](mpl.colors.Normalize(0, vmax, clip=True)(nppb[run].loc[idx].to_numpy())))
    mask = tf.loc[idx].isin(new).to_numpy()
    for a, c in [(axs[3], "black"), (axs[4], "white")]:
        a.add_collection(LineCollection(outline(sp, mask), colors=c, linewidths=0.5))
    if i == 0:
        for a, t in zip(axs, ["H&E", "RNA + space 0.75\n(joint, 27 slides)", f"{L} + space {WL}\n(joint, 27 slides)",
                              f"RNA + {L} + space {WL}\n(joint, 27 slides)",
                              "NPPB expression\n(outline: new domain)"]):
            a.set_title(t, fontweight="bold")
    if i == len(SLIDES) - 1:
        cax = axs[4].inset_axes([0.1, -0.08, 0.8, 0.04])
        cb = fig.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(0, vmax), cmap="magma"), cax=cax,
                          orientation="horizontal")
        cb.ax.tick_params(labelsize=5, length=1.5, width=0.4)
        cb.set_label("log1p(CP10k)", fontsize=5)
mid = gs[1].subgridspec(1, 2, width_ratios=[1, 1.6], wspace=0.25)
bot = gs[2].subgridspec(1, 3, wspace=0.45)
ax = fig.add_subplot(mid[0])
ax.imshow(src.to_numpy(), cmap="Greys", vmin=0, vmax=1, aspect="auto")
for (r_, c_), v in np.ndenumerate(src.to_numpy()):
    if v >= 0.05:
        ax.text(c_, r_, f"{v:.2f}", ha="center", va="center", fontsize=4.5, color="white" if v > 0.5 else "black")
ax.set_xticks(range(len(rna_doms)), rna_doms)
ax.set_yticks(range(len(tf_doms)), [d + (" (new)" if d in new else "") for d in tf_doms])
for t, d in zip(ax.get_xticklabels(), rna_doms):
    t.set_color(rna_col[d])
for t, d in zip(ax.get_yticklabels(), tf_doms):
    t.set_color(tf_col[d])
ax.set_xlabel("RNA + space domain")
ax.set_ylabel(f"RNA + {L} + space domain")
ax.set_title(f"a  Where each RNA + {L} domain\ncomes from (share of its spots)", fontweight="bold", loc="left")
for j, (tag, col_, doms, cols, title) in enumerate([("exp2_joint_harmony", COL, rna_doms, rna_col, "b  RNA + space"),
                                                    (SINGLE[0], SINGLE[1], tfo_doms, tfo_col, f"c  {L} + space"),
                                                    (COMBO[0], COMBO[1], tf_doms, tf_col, f"d  RNA + {L} + space")], start=1):
    t = pd.read_csv(OUT / f"{tag}_{col_}_tests.csv").set_index("domain").loc[doms]
    M = t[[f"mean_share_{z}" for z in ZONES]].to_numpy()
    a = fig.add_subplot(bot[j - 1])
    a.imshow(M, cmap="Blues", vmin=0, vmax=0.75, aspect="auto")
    for (r_, c_), v in np.ndenumerate(M):
        a.text(c_, r_, f"{v:.2f}", ha="center", va="center", fontsize=4.5, color="white" if v > 0.45 else "black")
    sig = []
    for d in doms:
        s = ""
        if t.loc[d, "padj_CTRL_gt_patients"] < 0.1:
            s = " ▲ctrl"
        elif min(t.loc[d, f"padj_{z}_gt_CTRL"] for z in ZONES[1:]) < 0.1:
            s = " ▲pat"
        sig.append(s)
    a.set_xticks(range(len(ZONES)), ZONES)
    newd = new_o if j == 2 else (new if j == 3 else [])
    a.set_yticks(range(len(doms)), [f"{d}{s}" + ("*" if d in newd else "") for d, s in zip(doms, sig)])
    for tk, d in zip(a.get_yticklabels(), doms):
        tk.set_color(cols[d])
    a.set_title(f"{title}\nmean share per zone", fontweight="bold", loc="left")
axl = fig.add_subplot(mid[1])
axl.axis("off")
h1 = [Patch(color=rna_col[d], label=f"{d} {desc_r[d]}") for d in rna_doms]
h2 = [Patch(color=tf_col[d], label=f"{d} {desc_t[d]}" + (" (new)" if d in new else "")) for d in tf_doms]
h3 = [Patch(color=tfo_col[d], label=f"{d} {desc_o[d]}" + (" (new)" if d in new_o else "")) for d in tfo_doms]
for k_, (h, t_) in enumerate([(h1, "RNA + space"), (h3, f"{L} + space"), (h2, f"RNA + {L} + space")]):
    lg = axl.legend(handles=h, title=t_, loc="upper left", bbox_to_anchor=(0.0 + 0.34 * k_, 1.0),
                    frameon=False, fontsize=5, title_fontsize=5.5, handlelength=0.9, labelspacing=0.3)
    axl.add_artist(lg)
for a in fig.axes:
    a.tick_params(labelsize=5, length=1.5, width=0.4)
fig.text(0.5, 0.075, "▲ctrl / ▲pat: enriched in control / in a patient zone vs control (BH-adjusted p < 0.1; slide = unit, "
         "4 control slides); * new domain (split out of an RNA domain)", ha="center", fontsize=5)
fig.savefig(OUT / f"{FIG}.png", dpi=250, bbox_inches="tight")
fig.savefig(OUT / f"{FIG}.pdf", dpi=300, bbox_inches="tight")
print(f"new domains in RNA + {L}:", new, {d: desc_t[d] for d in new}, f"| in {L} + space:", new_o, {d: desc_o[d] for d in new_o})
print("labels from", SINGLE, COMBO)
