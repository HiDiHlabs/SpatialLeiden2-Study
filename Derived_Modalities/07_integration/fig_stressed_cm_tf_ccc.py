"""Characteristic TFs and ligand-receptor (LR) pairs of the stressed-cardiomyocyte domain (joint RNA + TF +
space, exp3_joint_rna_tf r0.75_s0, domain D5 by size; NPPB / ANKRD1 high), and maps on Visium slides.

Ranking, within slide (so slide / patient effects cannot drive it): on every slide where the domain covers
>= 2 % of the spots, effect = (mean inside - mean outside) / SD of the slide, per feature; then the mean
effect over those slides and the number of slides with a positive effect. Features:
  TF   CollecTRI ULM score (decoupler), 622 TFs scored on all 27 slides (data/.../collectri_all27)
  LR   our LIANA+ local cosine score (data/.../ccc_all27_ours), the 286 pairs present on all 27 slides
Top features = highest mean effect among those positive on >= 75 % of the slides.
Caveat: the domain was found with TF activity as one of its layers, so TF enrichment is partly built in;
the LR scores were not used for this clustering (RNA + CCC finds an overlapping domain independently).

Figure: a bar charts of the top 12 TFs and LR pairs (mean effect; text = slides positive / slides);
b maps on control P7 (ACH003), BZ P3 (10X0026), IZ P3 (10X0017), IZ P9 (ACH0012): H&E | domain |
top 2 TFs | top 2 LR pairs | NPPB; the domain outlined on every map. Colour scales 1st-99th percentile
over the shown slides.

Usage (from Derived_Modalities/, env spatialleiden_joint):
    python 07_integration/fig_stressed_cm_tf_ccc.py data/kuppe_2022 results/kuppe_2022
Output: results/kuppe_2022/stressed_cm/fig_stressed_cm_tf_ccc.png/.pdf, tf_ranking.csv, lr_ranking.csv
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
OUT = RES / "stressed_cm"
OUT.mkdir(parents=True, exist_ok=True)
MM = 1 / 25.4
FS, FH, FL = 5, 6, 8
mpl.rcParams.update({"font.family": "Arial", "font.size": FS, "axes.titlesize": FH, "pdf.fonttype": 42})
SHOW = {"ACH003": "Control P7", "10X0026": "BZ P3", "10X0017": "IZ P3", "ACH0012": "IZ P9"}
PINK = "#e43fb1"

lab = pd.read_csv(RES / "exp3_joint_rna_tf" / "exp2_labels.csv.gz", index_col=0, dtype=str)["r0.75_s0"]
order = lab.value_counts().index
lab = lab.map({d: f"D{i + 1}" for i, d in enumerate(order)})
dom = lab == "D5"
run = lab.index.str.split("_", n=1).str[0]
runs = pd.unique(run)

pairs = (RES / "ccc_harmony_ours" / "ccc_pairs.txt").read_text().split()
tf_sets, tfs = {}, None
for r in runs:
    t = ad.read_h5ad(DATA / "modalities" / "collectri_all27" / f"collectri_ulm_{r}.h5ad")
    tf_sets[r] = t
    tfs = set(t.var_names) if tfs is None else tfs & set(t.var_names)
tfs = sorted(tfs)


def per_slide_effect(get):
    rows = []
    for r in runs:
        idx = lab.index[run == r]
        m = dom.loc[idx].to_numpy()
        if m.mean() < 0.02:
            continue
        X = get(r).loc[idx]
        sd = X.std().replace(0, np.nan)
        rows.append(((X[m].mean() - X[~m].mean()) / sd).rename(r))
    E = pd.concat(rows, axis=1)
    out = pd.DataFrame({"mean_effect": E.mean(1), "n_positive": (E > 0).sum(1), "n_slides": E.notna().sum(1)})
    out["frac_positive"] = out.n_positive / out.n_slides
    return out.sort_values("mean_effect", ascending=False), E


def tf_get(r):
    t = tf_sets[r][:, tfs]
    return pd.DataFrame(np.asarray(t.X), index=r + "_" + t.obs_names, columns=tfs)


ccc = {r: pd.read_csv(DATA / "modalities" / "ccc_all27_ours" / f"ccc_{r}.csv", index_col=0)[pairs] for r in runs}


def lr_get(r):
    c = ccc[r]
    return c.set_axis(r + "_" + c.index, axis=0)


tf_rank, _ = per_slide_effect(tf_get)
lr_rank, _ = per_slide_effect(lr_get)
tf_rank.to_csv(OUT / "tf_ranking.csv")
lr_rank.to_csv(OUT / "lr_ranking.csv")
tf_top = tf_rank[tf_rank.frac_positive >= 0.75].head(12)
lr_top = lr_rank[lr_rank.frac_positive >= 0.75].head(12)
print(f"domain on {int(tf_rank.n_slides.max())} slides (>= 2 %)")
print("TFs:\n", tf_top.round(2).to_string())
print("LR pairs:\n", lr_top.round(2).to_string())


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


# features on the shown slides
nppb = {}
for f in sorted((DATA / "processed1").glob("Visium_*.h5ad")):
    nm = set(ad.read_h5ad(f, backed="r").obs_names)
    for r in SHOW:
        if set(lab.index[run == r].str.split("_", n=1).str[1]) == nm:
            h = ad.read_h5ad(f)
            tot = np.asarray(h.raw.X.sum(1)).ravel()
            g = list(h.raw.var_names).index("NPPB")
            nppb[r] = pd.Series(np.log1p(np.asarray(h.raw.X[:, g].todense()).ravel() / tot * 1e4), index=r + "_" + h.obs_names)
feat = [(f"{t}\nTF activity", "tf", t, "viridis") for t in tf_top.index[:2]] + \
       [(f"{p.replace('^', ' → ')}\nLR score", "lr", p, "magma") for p in lr_top.index[:2]] + \
       [("NPPB\nexpression", "g", "NPPB", "magma")]


def values(r, kind, key):
    idx = lab.index[run == r]
    if kind == "tf":
        return tf_get(r).loc[idx, key].to_numpy()
    if kind == "lr":
        return lr_get(r).loc[idx, key].to_numpy()
    return nppb[r].loc[idx].to_numpy()


lims = {}
for title, kind, key, cm in feat:
    allv = np.concatenate([values(r, kind, key) for r in SHOW])
    lims[title] = (np.quantile(allv, 0.01), np.quantile(allv, 0.99))

fig = plt.figure(figsize=(183 * MM, 182 * MM))
gs = fig.add_gridspec(2, 1, height_ratios=[0.85, 2.3], hspace=0.3)
ga = gs[0].subgridspec(1, 2, wspace=0.6)
for j, (tab, title, colr) in enumerate([(tf_top, "Top TFs (CollecTRI activity)", "#2a9d8f"),
                                        (lr_top, "Top ligand → receptor pairs (local score)", "#9c3d54")]):
    ax = fig.add_subplot(ga[j])
    y = np.arange(len(tab))[::-1]
    ax.barh(y, tab.mean_effect, color=colr, height=0.7)
    ax.set_yticks(y, [i.replace("^", " → ") for i in tab.index], fontsize=FS)
    for yy, (_, rr) in zip(y, tab.iterrows()):
        ax.text(rr.mean_effect + 0.02, yy, f"{int(rr.n_positive)}/{int(rr.n_slides)}", va="center", fontsize=FS - 0.5)
    ax.set_xlabel("effect inside vs outside the domain\n(SD units, mean over slides)", fontsize=FS)
    ax.set_title(title, fontweight="bold", fontsize=FH, loc="left")
    ax.tick_params(labelsize=FS, length=1.5, width=0.4)
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)
    ax.set_xlim(0, tab.mean_effect.max() * 1.25)
    if j == 0:
        ax.text(-0.45, 1.05, "a", transform=ax.transAxes, fontweight="bold", fontsize=FL)
gb = gs[1].subgridspec(len(SHOW), 2 + len(feat), wspace=0.04, hspace=0.06)
for i, (r, label) in enumerate(SHOW.items()):
    idx = lab.index[run == r]
    sp = pd.read_csv(DATA / "visium_tables" / f"{r}_spots.csv").set_index("barcode").loc[idx.str.split("_", n=1).str[1]]
    xy = sp[["pxl_col_fullres", "pxl_row_fullres"]].to_numpy(float)
    m = dom.loc[idx].to_numpy()
    segs = outline(sp, m)
    ax = fig.add_subplot(gb[i, 0])
    img = np.asarray(Image.open(DATA / "visium_tables" / f"{r}_lowres.png").convert("RGB"))
    f = json.loads((DATA / "visium_tables" / f"{r}_scalefactors.json").read_text())["tissue_lowres_scalef"]
    ax.imshow(img)
    ax.set_xlim(xy[:, 0].min() * f - 8, xy[:, 0].max() * f + 8)
    ax.set_ylim(xy[:, 1].max() * f + 8, xy[:, 1].min() * f - 8)
    ax.axis("off")
    ax.add_collection(LineCollection(segs * f, colors="black", linewidths=0.4))
    ax.text(-0.06, 0.5, f"{label}\n({m.mean():.0%} domain)", transform=ax.transAxes, rotation=90, ha="right",
            va="center", fontweight="bold", fontsize=FS)
    if i == 0:
        ax.set_title("H&E", fontweight="bold", fontsize=FH)
        ax.text(-0.5, 1.12, "b", transform=ax.transAxes, fontweight="bold", fontsize=FL)
    ax = fig.add_subplot(gb[i, 1])
    draw(ax, xy, [PINK if v else "#E6E6E6" for v in m])
    ax.add_collection(LineCollection(segs, colors="black", linewidths=0.4))
    if i == 0:
        ax.set_title("Stressed CM\ndomain (spots)", fontweight="bold", fontsize=FH)
    for j, (title, kind, key, cm) in enumerate(feat, start=2):
        ax = fig.add_subplot(gb[i, j])
        nrm = mpl.colors.Normalize(*lims[title], clip=True)
        draw(ax, xy, mpl.colormaps[cm](nrm(values(r, kind, key))))
        ax.add_collection(LineCollection(segs, colors="white", linewidths=0.4))
        if i == 0:
            ax.set_title(title, fontweight="bold", fontsize=FH)
        if i == len(SHOW) - 1:
            cax = ax.inset_axes([0.1, -0.1, 0.8, 0.05])
            cb = fig.colorbar(mpl.cm.ScalarMappable(norm=nrm, cmap=cm), cax=cax, orientation="horizontal")
            cb.ax.tick_params(labelsize=FS, length=1.5, width=0.4)
            cb.outline.set_linewidth(0.4)
fig.savefig(OUT / "fig_stressed_cm_tf_ccc.png", dpi=250, bbox_inches="tight")
fig.savefig(OUT / "fig_stressed_cm_tf_ccc.pdf", dpi=300, bbox_inches="tight")
