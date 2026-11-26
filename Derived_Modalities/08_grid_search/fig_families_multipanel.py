"""Summary multipanel of the grid-search region families (Nature format, 183 x 170 mm, Arial 5/6/8 pt).

Families (consensus spots from candidate_families.py), grouped by theme:
  3 myocardium IZ, 7 myocardium BZ, 33 stressed CM, 0 inflammatory macrophage, 13 lipid niche,
  2 fibrotic matrix (THBS4/CD36), 6 fibrotic complement stroma
a  TFs: union of each family's top 3 consistent TFs (fig_region_tf_ccc.py rankings); dot colour = mean
   within-slide effect (inside vs outside, SD units), dot size = share of slides with a positive effect
b  the same for ligand -> receptor pairs (our LIANA+ local scores)
c  mean share of the slide per zone (slide = unit)
d  maps on each family's representative slide (largest share among slides of its main zone):
   H&E with the region outlined | region | top TF activity | top LR score
   (TF and LR colour scales: 1st-99th percentile of that slide)
e  stressed-CM domain across the grid: share of runs (all K, 5 seeds) with a domain at Jaccard >= 0.5 to the
   stressed-CM domain of 24 Nov (exp3_joint_rna_tf r0.75_s0, D5), per layer configuration and spatial weight

Usage (from Derived_Modalities/, env spatialleiden_joint):
    python 08_grid_search/fig_families_multipanel.py data/kuppe_2022 results/kuppe_2022
Output: results/kuppe_2022/grid_search/fig_families_multipanel.png/.pdf
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
G = RES / "grid_search"
MM = 1 / 25.4
FS, FH, FL = 5, 6, 8
mpl.rcParams.update({"font.family": "Arial", "font.size": FS, "axes.titlesize": FH, "pdf.fonttype": 42,
                     "axes.linewidth": 0.4})
FAMS = [(3, "Myocardium (IZ)"), (7, "Myocardium (BZ)"), (33, "Stressed CM"), (0, "Inflammatory Mφ"),
        (13, "Lipid niche"), (2, "Fibrotic matrix"), (6, "Fibrotic complement")]
ZONES = ["CTRL", "RZ", "BZ", "IZ", "FZ"]
PINK = "#e43fb1"

spots = pd.read_csv(G / "spots.csv.gz", index_col=0)
cons = pd.read_csv(G / "family_consensus_spots.csv.gz", index_col=0).loc[spots.index]
run = spots["sample"].to_numpy()
zone_of = spots.groupby("sample").major_labl.first()


def top(f, kind, n=3):
    t = pd.read_csv(G / "regions" / f"family_{f}" / f"{kind}_ranking.csv", index_col=0)
    c = t[t.frac_positive >= 0.75]
    return (c if len(c) else t).head(n), t


tf_rank = {f: top(f, "tf")[1] for f, _ in FAMS}
lr_rank = {f: top(f, "lr")[1] for f, _ in FAMS}
tf_list = list(dict.fromkeys(x for f, _ in FAMS for x in top(f, "tf")[0].index))
lr_list = list(dict.fromkeys(x for f, _ in FAMS for x in top(f, "lr", 2)[0].index))


def dotplot(ax, ranks, feats, title, cmap):
    eff = np.array([[ranks[f].mean_effect.get(x, np.nan) for x in feats] for f, _ in FAMS])
    frac = np.array([[ranks[f].frac_positive.get(x, np.nan) for x in feats] for f, _ in FAMS])
    lim = np.nanquantile(np.abs(eff), 0.98)
    yy, xx = np.meshgrid(range(len(FAMS)), range(len(feats)), indexing="ij")
    sca = ax.scatter(xx.ravel(), yy.ravel(), s=(np.nan_to_num(frac.ravel()) ** 2) * 14 + 0.5,
                     c=eff.ravel(), cmap=cmap, vmin=-lim, vmax=lim, linewidths=0)
    ax.set_xticks(range(len(feats)), [x.replace("^", "→").replace("_", "/").replace("ITGA7/ITGB1", "ITGA7/B1") for x in feats], rotation=90, fontsize=FS)
    ax.set_yticks(range(len(FAMS)), [f"{n} ({f})" for f, n in FAMS], fontsize=FS)
    ax.set_xlim(-0.6, len(feats) - 0.4)
    ax.set_ylim(len(FAMS) - 0.5, -0.5)
    ax.set_title(title, fontweight="bold", fontsize=FH, loc="left")
    ax.tick_params(length=1.5, width=0.4, pad=1)
    ax.grid(color="#eeeeee", lw=0.3)
    ax.set_axisbelow(True)
    return sca


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
    xy = sp[["pxl_col_fullres", "pxl_row_fullres"]].to_numpy(float)
    grid = np.c_[sp.array_col, sp.array_row].astype(float)
    coef, *_ = np.linalg.lstsq(np.c_[grid, np.ones(len(grid))], xy, rcond=None)
    pos = {(int(a), int(b)): n for n, (a, b) in enumerate(zip(sp.array_row, sp.array_col))}
    segs = []
    for n in np.flatnonzero(mask):
        a, b = int(sp.array_row.iloc[n]), int(sp.array_col.iloc[n])
        for da, db in [(0, 2), (0, -2), (1, 1), (1, -1), (-1, 1), (-1, -1)]:
            q = pos.get((a + da, b + db))
            if q is not None and mask[q]:
                continue
            v = np.array([db, da]) @ coef[:2]
            mid = xy[n] + v / 2
            perp = np.array([-v[1], v[0]]) / (2 * np.sqrt(3))
            segs.append([mid - perp, mid + perp])
    return np.array(segs) if segs else np.zeros((0, 2, 2))


fig = plt.figure(figsize=(183 * MM, 170 * MM))
outer = fig.add_gridspec(2, 1, height_ratios=[1, 2.05], hspace=0.5)
top_ = outer[0].subgridspec(1, 4, width_ratios=[len(tf_list), len(lr_list), 4.6, 0.6], wspace=0.12)
ax_a = fig.add_subplot(top_[0])
sa = dotplot(ax_a, tf_rank, tf_list, "a  Characteristic TFs", "PuOr_r")
ax_b = fig.add_subplot(top_[1])
sb = dotplot(ax_b, lr_rank, lr_list, "b  Characteristic ligand → receptor pairs", "PuOr_r")
ax_b.set_yticklabels([])
ax_a.text(-0.02, 1.0, "", transform=ax_a.transAxes)
# c zone shares
ax_c = fig.add_subplot(top_[2])
shares = np.zeros((len(FAMS), len(ZONES)))
for i, (f, _) in enumerate(FAMS):
    per = pd.Series(cons[f"family_{f}"].to_numpy()).groupby(run).mean()
    shares[i] = [per[zone_of.index[zone_of == z]].mean() for z in ZONES]
ax_c.imshow(shares, cmap="Blues", vmin=0, vmax=max(0.3, shares.max()), aspect="auto")
for (a_, b_), v in np.ndenumerate(shares):
    ax_c.text(b_, a_, f"{100 * v:.0f}" if v >= 0.005 else "", ha="center", va="center", fontsize=FS - 0.5,
              color="white" if v > 0.18 else "black")
ax_c.set_xticks(range(len(ZONES)), ["ctrl" if z == "CTRL" else z for z in ZONES], fontsize=FS)
ax_c.set_yticks([])
ax_c.set_title("c  % of slide", fontweight="bold", fontsize=FH, loc="left")
ax_c.tick_params(length=1.5, width=0.4, pad=1)
# legend for dot plots
ax_l = fig.add_subplot(top_[3])
ax_l.axis("off")
cax = ax_l.inset_axes([0.1, 0.45, 0.35, 0.5])
cb = fig.colorbar(sb, cax=cax)
cb.ax.tick_params(labelsize=FS - 0.5, length=1.5, width=0.4)
cb.set_label("effect (SD)", fontsize=FS)
cb.outline.set_linewidth(0.4)
for k_, fr in enumerate([0.5, 0.75, 1.0]):
    ax_l.scatter(0.3, 0.3 - k_ * 0.12, s=fr ** 2 * 14 + 0.5, c="grey", transform=ax_l.transAxes, clip_on=False)
    ax_l.text(0.55, 0.3 - k_ * 0.12, f"{fr:.0%}", transform=ax_l.transAxes, va="center", fontsize=FS - 0.5)
ax_l.text(0.0, 0.4, "slides +", transform=ax_l.transAxes, fontsize=FS - 0.5)

# d maps + e grid check
bot = outer[1].subgridspec(1, 2, width_ratios=[2.2, 1], wspace=0.38)
gm = bot[0].subgridspec(len(FAMS), 4, wspace=0.08, hspace=0.38)
tf_cache, lr_cache = {}, {}
for i, (f, name) in enumerate(FAMS):
    mask_all = cons[f"family_{f}"].to_numpy()
    per = pd.Series(mask_all).groupby(run).mean()
    main_zone = pd.Series(spots.major_labl.to_numpy()[mask_all]).mode().iloc[0]
    r = per[zone_of.index[zone_of == main_zone]].idxmax()
    idx = spots.index[run == r]
    sp = pd.read_csv(DATA / "visium_tables" / f"{r}_spots.csv").set_index("barcode").loc[idx.str.split("_", n=1).str[1]]
    xy = sp[["pxl_col_fullres", "pxl_row_fullres"]].to_numpy(float)
    m = cons.loc[idx, f"family_{f}"].to_numpy()
    segs = outline(sp, m)
    tfn = top(f, "tf", 1)[0].index[0]
    lrn = top(f, "lr", 1)[0].index[0]
    t = ad.read_h5ad(DATA / "modalities" / "collectri_all27" / f"collectri_ulm_{r}.h5ad")
    tv = pd.Series(np.asarray(t[:, tfn].X).ravel(), index=r + "_" + t.obs_names).loc[idx].to_numpy()
    c = pd.read_csv(DATA / "modalities" / "ccc_all27_ours" / f"ccc_{r}.csv", index_col=0)
    lv = c.set_axis(r + "_" + c.index, axis=0).loc[idx, lrn].to_numpy()
    ax = fig.add_subplot(gm[i, 0])
    img = np.asarray(Image.open(DATA / "visium_tables" / f"{r}_lowres.png").convert("RGB"))
    sf = json.loads((DATA / "visium_tables" / f"{r}_scalefactors.json").read_text())["tissue_lowres_scalef"]
    ax.imshow(img)
    ax.set_xlim(xy[:, 0].min() * sf - 8, xy[:, 0].max() * sf + 8)
    ax.set_ylim(xy[:, 1].max() * sf + 8, xy[:, 1].min() * sf - 8)
    ax.axis("off")
    ax.add_collection(LineCollection(segs * sf, colors="black", linewidths=0.35))
    pat = sp.patient.iloc[0]
    ax.text(-0.08, 0.5, f"{name}\n{zone_of[r]} {pat}, {m.mean():.0%}", transform=ax.transAxes, rotation=0,
            ha="right", va="center", fontsize=FS)
    ax = fig.add_subplot(gm[i, 1])
    draw(ax, xy, [PINK if v else "#E6E6E6" for v in m])
    ax.add_collection(LineCollection(segs, colors="black", linewidths=0.35))
    for j, (vals, lab_, cm) in enumerate([(tv, f"{tfn}", "viridis"), (lv, lrn.replace("^", "→").replace("_", "/"), "magma")], start=2):
        ax = fig.add_subplot(gm[i, j])
        nrm = mpl.colors.Normalize(np.quantile(vals, 0.01), np.quantile(vals, 0.99), clip=True)
        draw(ax, xy, mpl.colormaps[cm](nrm(vals)))
        ax.add_collection(LineCollection(segs, colors="white", linewidths=0.35))
        ax.text(0.5, -0.04, lab_.replace("ITGA7/ITGB1", "ITGA7/B1"), transform=ax.transAxes, ha="center", va="top",
                fontsize=FS - 0.5, fontstyle="italic" if j == 2 else "normal")
    if i == 0:
        for j, t_ in enumerate(["H&E", "Region", "Top TF\n(activity)", "Top LR pair\n(local score)"]):
            fig.axes[-4 + j].set_title(t_, fontweight="bold", fontsize=FH)
        fig.axes[-4].text(-1.9, 1.35, "d", transform=fig.axes[-4].transAxes, fontweight="bold", fontsize=FL)
# e
d = pd.read_csv(G / "stressed_cm_in_grid.csv")
CFG = ["rna", "rna+tf_0.25", "rna+tf_0.5", "rna+tf_1.0", "rna+ccc_0.25", "rna+ccc_0.5", "rna+ccc_1.0",
       "rna+tf+ccc_0.25", "rna+tf+ccc_0.5", "rna+tf+ccc_1.0"]
W = [0.0, 0.25, 0.75, 1.2, 2.0]
mat = d.assign(hit=d.best_j >= 0.5).groupby(["config", "w"]).hit.mean().unstack().reindex(index=CFG, columns=W)
ge = bot[1].subgridspec(2, 1, height_ratios=[1.6, 1], hspace=0.6)
ax_e = fig.add_subplot(ge[0])
ax_e.imshow(mat.to_numpy(), cmap="Purples", vmin=0, vmax=0.5, aspect="auto")
for (a_, b_), v in np.ndenumerate(mat.to_numpy()):
    ax_e.text(b_, a_, f"{v:.0%}", ha="center", va="center", fontsize=FS - 0.5, color="white" if v > 0.3 else "black")
lab_cfg = [c.replace("rna+", "RNA+").replace("rna", "RNA").replace("tf", "TF").replace("ccc", "CCC").replace("_", " w ")
           for c in CFG]
ax_e.set_yticks(range(len(CFG)), lab_cfg, fontsize=FS)
ax_e.set_xticks(range(len(W)), [str(w) for w in W], fontsize=FS)
ax_e.set_xlabel("spatial weight", fontsize=FS)
ax_e.tick_params(length=1.5, width=0.4, pad=1)
ax_e.set_title("e  Stressed-CM domain recovered\n(share of runs, Jaccard ≥ 0.5)", fontweight="bold", fontsize=FH, loc="left")
ax_k = fig.add_subplot(ge[1])
for c_, col, lab_ in [("rna", "#555555", "RNA"), ("rna+tf_0.5", "#2a9d8f", "RNA+TF w 0.5"), ("rna+ccc_0.5", "#9c3d54", "RNA+CCC w 0.5")]:
    med = d[d.config == c_].groupby("K").best_j.median()
    ax_k.plot(med.index, med.values, marker="o", ms=2, lw=0.8, color=col, label=lab_)
ax_k.axhline(0.5, color="black", lw=0.4, ls="--")
ax_k.set_xticks([9, 15, 25, 40])
ax_k.set_xlabel("target K", fontsize=FS)
ax_k.set_ylabel("median best Jaccard", fontsize=FS)
ax_k.legend(frameon=False, fontsize=FS - 0.5, loc="lower left", bbox_to_anchor=(0.0, 1.02), ncol=3, handlelength=1.0, columnspacing=0.6)
ax_k.tick_params(labelsize=FS, length=1.5, width=0.4)
for s_ in ("top", "right"):
    ax_k.spines[s_].set_visible(False)
fig.savefig(G / "fig_families_multipanel.png", dpi=300, bbox_inches="tight")
fig.savefig(G / "fig_families_multipanel.pdf", dpi=300, bbox_inches="tight")
print("TFs", tf_list, "\nLR", lr_list)
