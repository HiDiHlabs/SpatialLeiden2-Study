"""Multipanel of the stressed-cardiomyocyte region (Nature format, 183 x 170 mm, Arial 5/6/8 pt).

Region = union of the consensus spots of the four grid-search families with the stressed-CM signature
(NPPB and ANKRD1 higher inside on every slide; HSF / UPR TFs): families 33, 37, 41, 42
(08_grid_search/candidate_families.py). Within-slide effects as in fig_region_tf_ccc.py: on every slide where
the region covers >= 2 % of the spots, (mean inside - mean outside) / SD of the slide, averaged over slides.
a  maps on four slides (two border zone, one ischaemic, one control): H&E with the region outlined, region,
   NPPB, ANKRD1 (the authors' vCM3 markers), top TF, top LR pair
b  the authors' cardiomyocyte-state markers (NPPB, ANKRD1, MYH7) and other cell-type markers, effect per slide
c  share of the slide per zone (one dot per slide)
d  top TFs, e top ligand -> receptor pairs (bars: mean effect; text: slides positive)
f  how the four families were found: their instances per layer configuration (grid search)
g  the region split across the authors' cell-type and molecular niches (27 slides)

Usage (from Derived_Modalities/, env spatialleiden_joint):
    python 08_grid_search/fig_stressed_cm_multipanel.py data/kuppe_2022 results/kuppe_2022
Output: results/kuppe_2022/grid_search/fig_stressed_cm_multipanel.png/.pdf, stressed_cm_region_effects.csv
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

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))
from plot_patient_region_umap import PATIENT_COLOURS, REGION_COLOURS  # noqa: E402

DATA, RES = Path(sys.argv[1]), Path(sys.argv[2])
G = RES / "grid_search"
MM = 1 / 25.4
FS, FH, FL = 5, 6, 8
mpl.rcParams.update({"font.family": "Arial", "font.size": FS, "axes.titlesize": FH, "pdf.fonttype": 42,
                     "axes.linewidth": 0.4})
FAMS = [33, 37, 41, 42]
PINK = "#e43fb1"
ZONES = ["CTRL", "RZ", "BZ", "IZ", "FZ"]
ZCOL = {"CTRL": REGION_COLOURS["control"], **{z: REGION_COLOURS[z] for z in ZONES[1:]}}
GENES = ["NPPB", "ANKRD1", "MYH7", "TTN", "POSTN", "COL1A1", "SPP1", "CD163", "PECAM1"]

spots = pd.read_csv(G / "spots.csv.gz", index_col=0)
cons = pd.read_csv(G / "family_consensus_spots.csv.gz", index_col=0).loc[spots.index]
region = cons[[f"family_{f}" for f in FAMS]].any(axis=1).to_numpy()
run = spots["sample"].to_numpy()
runs = pd.unique(run)
zone_of = spots.groupby("sample").major_labl.first()
pat_of = spots.groupby("sample").patient.first()

# features: genes (log1p CP10k), TF (ULM), LR (local score)
tables = {}
for t in (DATA / "visium_tables").glob("*_spots.csv"):
    tables[t.name[:-10]] = set(pd.read_csv(t, usecols=["barcode"]).barcode)
E = []
for f in sorted((DATA / "processed1").glob("Visium_*.h5ad")):
    h = ad.read_h5ad(f)
    r = [k for k, b in tables.items() if b == set(h.obs_names)][0]
    tot = np.asarray(h.raw.X.sum(1)).ravel()
    nm = list(h.raw.var_names)
    E.append(pd.DataFrame({g: (np.log1p(np.asarray(h.raw.X[:, nm.index(g)].todense()).ravel() / tot * 1e4)
                               if g in nm else 0.0) for g in GENES}, index=r + "_" + h.obs_names))
E = pd.concat(E).loc[spots.index]
pairs = (RES / "ccc_harmony_ours" / "ccc_pairs.txt").read_text().split()
TF, LR, tfs = {}, {}, None
for r in runs:
    t = ad.read_h5ad(DATA / "modalities" / "collectri_all27" / f"collectri_ulm_{r}.h5ad")
    TF[r] = pd.DataFrame(np.asarray(t.X), index=r + "_" + t.obs_names, columns=t.var_names)
    tfs = set(t.var_names) if tfs is None else tfs & set(t.var_names)
    c = pd.read_csv(DATA / "modalities" / "ccc_all27_ours" / f"ccc_{r}.csv", index_col=0)[pairs]
    LR[r] = c.set_axis(r + "_" + c.index, axis=0)
tfs = sorted(tfs)


def effects(get):
    cols = []
    for r in runs:
        s = run == r
        m = region[s]
        if m.mean() < 0.02:
            continue
        X = get(r).loc[spots.index[s]]
        sd = X.std().replace(0, np.nan)
        cols.append(((X[m].mean() - X[~m].mean()) / sd).rename(r))
    Ef = pd.concat(cols, axis=1)
    return pd.DataFrame({"mean_effect": Ef.mean(1), "n_positive": (Ef > 0).sum(1), "n_slides": Ef.notna().sum(1)}), Ef


g_eff, g_mat = effects(lambda r: E.loc[spots.index[run == r]])
t_eff, _ = effects(lambda r: TF[r][tfs])
l_eff, _ = effects(lambda r: LR[r])
pd.concat([g_eff.assign(kind="gene"), t_eff.assign(kind="TF"), l_eff.assign(kind="LR")]).to_csv(
    G / "stressed_cm_region_effects.csv")
nsl = int(g_eff.n_slides.max())
t_top = t_eff[t_eff.n_positive >= 0.75 * t_eff.n_slides].sort_values("mean_effect", ascending=False).head(10)
l_top = l_eff[l_eff.n_positive >= 0.75 * l_eff.n_slides].sort_values("mean_effect", ascending=False).head(10)
print(f"region: {region.sum()} spots, on {nsl} slides (>= 2 %)")
print(g_eff.round(2).to_string())
print(t_top.round(2).to_string())
print(l_top.round(2).to_string())


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


share = pd.Series(region).groupby(run).mean()
bz = share[zone_of[share.index] == "BZ"].sort_values(ascending=False).index[:2].tolist()
iz = share[zone_of[share.index] == "IZ"].sort_values(ascending=False).index[:1].tolist()
ct = share[zone_of[share.index] == "CTRL"].sort_values(ascending=False).index[:1].tolist()
SHOW = bz + iz + ct
tf1, lr1 = t_top.index[0], l_top.index[0]
cols_a = [("H&E", None), ("Region", None), ("NPPB", "g"), ("ANKRD1", "g"), (tf1, "tf"), (lr1, "lr")]

fig = plt.figure(figsize=(183 * MM, 170 * MM))
outer = fig.add_gridspec(2, 1, height_ratios=[1.25, 1], hspace=0.32)
ga = outer[0].subgridspec(len(SHOW), len(cols_a), wspace=0.04, hspace=0.1)
lims = {}
for name, kind in cols_a[2:]:
    vals = np.concatenate([(E.loc[spots.index[run == r], name] if kind == "g" else
                            (TF[r].loc[spots.index[run == r], name] if kind == "tf" else
                             LR[r].loc[spots.index[run == r], name])).to_numpy() for r in SHOW])
    lims[name] = (np.quantile(vals, 0.01), np.quantile(vals, 0.99))
for i, r in enumerate(SHOW):
    idx = spots.index[run == r]
    sp = pd.read_csv(DATA / "visium_tables" / f"{r}_spots.csv").set_index("barcode").loc[idx.str.split("_", n=1).str[1]]
    xy = sp[["pxl_col_fullres", "pxl_row_fullres"]].to_numpy(float)
    m = region[run == r]
    segs = outline(sp, m)
    for j, (name, kind) in enumerate(cols_a):
        ax = fig.add_subplot(ga[i, j])
        if name == "H&E":
            img = np.asarray(Image.open(DATA / "visium_tables" / f"{r}_lowres.png").convert("RGB"))
            sf = json.loads((DATA / "visium_tables" / f"{r}_scalefactors.json").read_text())["tissue_lowres_scalef"]
            ax.imshow(img)
            ax.set_xlim(xy[:, 0].min() * sf - 8, xy[:, 0].max() * sf + 8)
            ax.set_ylim(xy[:, 1].max() * sf + 8, xy[:, 1].min() * sf - 8)
            ax.axis("off")
            ax.add_collection(LineCollection(segs * sf, colors="black", linewidths=0.35))
            ax.text(-0.06, 0.5, f"{zone_of[r]} {pat_of[r]}\n{m.mean():.0%} region", transform=ax.transAxes,
                    ha="right", va="center", fontsize=FS)
            if i == 0:
                ax.set_title("H&E", fontweight="bold", fontsize=FH)
            continue
        if name == "Region":
            draw(ax, xy, [PINK if v else "#E6E6E6" for v in m])
            ax.add_collection(LineCollection(segs, colors="black", linewidths=0.35))
        else:
            v = (E.loc[idx, name] if kind == "g" else TF[r].loc[idx, name] if kind == "tf" else LR[r].loc[idx, name]).to_numpy()
            cm = "magma" if kind != "tf" else "viridis"
            nrm = mpl.colors.Normalize(*lims[name], clip=True)
            draw(ax, xy, mpl.colormaps[cm](nrm(v)))
            ax.add_collection(LineCollection(segs, colors="white", linewidths=0.35))
            if i == len(SHOW) - 1:
                cax = ax.inset_axes([0.1, -0.1, 0.8, 0.05])
                cb = fig.colorbar(mpl.cm.ScalarMappable(norm=nrm, cmap=cm), cax=cax, orientation="horizontal")
                cb.ax.tick_params(labelsize=FS - 0.5, length=1.5, width=0.4)
                cb.outline.set_linewidth(0.4)
        if i == 0:
            lab = {"g": f"{name}\n(expression)", "tf": f"{name}\n(TF activity)",
                   "lr": f"{name.replace('^', ' → ').replace('_', '/')}\n(LR score)"}.get(kind, name)
            ax.set_title(lab, fontweight="bold", fontsize=FH, fontstyle="italic" if kind in ("g", "tf") else "normal")
    if i == 0:
        fig.axes[-len(cols_a)].text(-0.75, 1.25, "a", transform=fig.axes[-len(cols_a)].transAxes, fontweight="bold", fontsize=FL)

bot = outer[1].subgridspec(1, 5, width_ratios=[1.15, 0.75, 1.0, 1.15, 1.2], wspace=0.95)
# b genes
ax = fig.add_subplot(bot[0])
gg = g_eff.loc[GENES]
y = np.arange(len(GENES))[::-1]
ax.barh(y, gg.mean_effect, color=[PINK if g in ("NPPB", "ANKRD1", "MYH7") else "#9e9e9e" for g in GENES], height=0.7)
for yy, (g, rr) in zip(y, gg.iterrows()):
    ax.text(max(rr.mean_effect, 0) + 0.03, yy, f"{int(rr.n_positive)}/{int(rr.n_slides)}", va="center", fontsize=FS - 0.5)
ax.set_yticks(y, GENES, fontstyle="italic")
ax.axvline(0, color="black", lw=0.4)
ax.set_xlabel("effect inside vs outside (SD)")
ax.set_title("b  Marker genes", fontweight="bold", loc="left")
# c zone dots
ax = fig.add_subplot(bot[1])
rng = np.random.default_rng(0)
for k, z in enumerate(ZONES):
    v = share[zone_of[share.index] == z].to_numpy()
    ax.scatter(k + rng.uniform(-0.15, 0.15, len(v)), v * 100, s=4, color=ZCOL[z], linewidths=0)
    ax.plot([k - 0.25, k + 0.25], [v.mean() * 100] * 2, color="black", lw=0.6)
ax.set_xticks(range(len(ZONES)), ["ctrl" if z == "CTRL" else z for z in ZONES], rotation=90)
ax.set_ylabel("% of slide")
ax.set_title("c  Zones", fontweight="bold", loc="left")
# d TFs, e LR
for pos, tab, title, colr, kk in [(2, t_top, "d  Top TFs", "#2a9d8f", "tf"), (3, l_top, "e  Top LR pairs", "#9c3d54", "lr")]:
    ax = fig.add_subplot(bot[pos])
    y = np.arange(len(tab))[::-1]
    ax.barh(y, tab.mean_effect, color=colr, height=0.7)
    for yy, (_, rr) in zip(y, tab.iterrows()):
        ax.text(rr.mean_effect - 0.02, yy, f"{int(rr.n_positive)}/{int(rr.n_slides)}", va="center", ha="right",
                fontsize=FS - 0.5, color="white")
    labs = [x.replace("^", "→").replace("_", "/").replace("ITGA5/ITGB1", "ITGA5/B1").replace("ITGA7/ITGB1", "ITGA7/B1").replace("ITGAV/ITGB1", "ITGAV/B1") for x in tab.index]
    ax.set_yticks(y, labs, fontstyle="italic" if kk == "tf" else "normal")
    ax.set_xlim(0, tab.mean_effect.max() * 1.08)
    ax.set_xlabel("effect (SD)")
    ax.set_title(title, fontweight="bold", loc="left")
# f layers + g authors' niches
gfg = bot[4].subgridspec(2, 1, hspace=0.9, height_ratios=[1, 1.3])
ax = fig.add_subplot(gfg[0])
cand = pd.read_csv(G / "candidates_with_family.csv")
cc = cand[cand.family.isin(FAMS)]
kind = np.select([cc.config.str.contains("tf") & cc.config.str.contains("ccc"), cc.config.str.contains("ccc"),
                  cc.config.str.contains("tf")], ["RNA+TF+CCC", "RNA+CCC", "RNA+TF"], "RNA")
tab = pd.crosstab(cc.family, kind).reindex(index=FAMS, columns=["RNA+TF", "RNA+CCC", "RNA+TF+CCC"], fill_value=0)
left = np.zeros(len(FAMS))
for c_, col in zip(tab.columns, ["#2a9d8f", "#9c3d54", "#6d5d8a"]):
    ax.barh(range(len(FAMS)), tab[c_], left=left, color=col, height=0.7, label=c_)
    left += tab[c_].to_numpy()
ax.set_yticks(range(len(FAMS)), [f"fam. {f}" for f in FAMS])
ax.set_ylim(len(FAMS) - 0.5, -0.5)
ax.set_xlabel("grid instances")
ax.legend(frameon=False, fontsize=FS - 0.5, loc="lower left", bbox_to_anchor=(0, 1.02), ncol=3, handlelength=0.8,
          columnspacing=0.5)
ax.set_title("f  Found with", fontweight="bold", loc="left", pad=10)
ax = fig.add_subplot(gfg[1])
STALLION = list(PATIENT_COLOURS.values())
for k_, (col, name) in enumerate([("celltype_niche", "cell-type"), ("molecular_niche", "molecular")]):
    order = sorted(spots[col].unique(), key=lambda x: int(x.split("_")[-1]))
    shr = pd.Series(spots[col].to_numpy()[region]).value_counts(normalize=True).reindex(order, fill_value=0)
    left = 0
    for o, v in shr.items():
        ax.barh(k_, v, left=left, color=dict(zip(order, STALLION))[o], height=0.6)
        if v > 0.12:
            ax.text(left + v / 2, k_, o.split("_")[-1], ha="center", va="center", fontsize=FS - 0.5, color="white")
        left += v
ax.set_yticks([0, 1], ["cell-type", "molecular"])
ax.set_ylim(1.5, -0.5)
ax.set_xlim(0, 1)
ax.set_xlabel("share of region's spots")
ax.set_title("g  Authors' niches", fontweight="bold", loc="left")
for a in fig.axes:
    a.tick_params(labelsize=FS, length=1.5, width=0.4, pad=1)
    if a.axison:
        for s_ in ("top", "right"):
            a.spines[s_].set_visible(False)
fig.savefig(G / "fig_stressed_cm_multipanel.png", dpi=300, bbox_inches="tight")
fig.savefig(G / "fig_stressed_cm_multipanel.pdf", dpi=300, bbox_inches="tight")
