"""How one region family is formed from grid-search domain instances (candidate_families.py).

a  every instance of the family on one slide (the slide where the consensus region is largest), labelled with
   its run: layer configuration, target K, spatial weight, seed
b  how often each spot is part of an instance (share of instances), and the consensus region (>= 50 %)
c  pairwise Jaccard of the instances over ALL spots of all 27 slides; instances are linked when >= 0.5
   (single linkage), which is how they were grouped into one family
d  the authors' annotations on the same slide (cell-type niches, molecular niches; authors' palette) and,
   over all 27 slides, how the consensus spots split across the authors' niches
Families with many instances: panel a shows at most 18, sampled evenly across configurations (seed 0 of
the sampler); the matrix in c covers all instances (numbers printed only up to 15).

Usage (from Derived_Modalities/, env spatialleiden_joint):
    python 08_grid_search/fig_family_grouping.py data/kuppe_2022 results/kuppe_2022 7
Output: results/kuppe_2022/grid_search/regions/family_<n>/fig_family_<n>_grouping.png/.pdf
"""

import json
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import EllipseCollection, LineCollection
from PIL import Image
from scipy.spatial import cKDTree

DATA, RES, FAM = Path(sys.argv[1]), Path(sys.argv[2]), int(sys.argv[3])
G = RES / "grid_search"
OUT = G / "regions" / f"family_{FAM}"
OUT.mkdir(parents=True, exist_ok=True)
MM = 1 / 25.4
FS, FH, FL = 5, 6, 8
mpl.rcParams.update({"font.family": "Arial", "font.size": FS, "axes.titlesize": FH, "pdf.fonttype": 42})
PINK = "#e43fb1"

spots = pd.read_csv(G / "spots.csv.gz", index_col=0)
run = spots["sample"].to_numpy()
cand = pd.read_csv(G / "candidates.csv").reset_index(drop=True)
# re-derive the family membership exactly as candidate_families.py
import scipy.sparse as sp_
from scipy.sparse.csgraph import connected_components
rows, cols = [], []
masks = []
for i, r in cand.iterrows():
    lab = np.load(G / "runs" / f"{r.run_id}.npz")["labels"]
    idx = np.flatnonzero(lab == r.domain)
    rows += [i] * len(idx)
    cols += idx.tolist()
M = sp_.csr_matrix((np.ones(len(rows), np.float32), (rows, cols)), shape=(len(cand), len(spots)))
inter = (M @ M.T).toarray()
size = np.asarray(M.sum(1)).ravel()
jac = inter / (size[:, None] + size[None, :] - inter)
_, fam = connected_components(sp_.csr_matrix(jac >= 0.5), directed=False)
cand["family"] = fam
members = cand[cand.family == FAM].sort_values(["config", "target_k", "spatial_weight", "seed"])
idx_m = members.index.to_numpy()
inst = np.asarray(M[idx_m].todense()).astype(bool)
freq = inst.mean(0)
consensus = freq >= 0.5
r_show = pd.Series(consensus).groupby(run).mean().idxmax()
sel = run == r_show
sp = pd.read_csv(DATA / "visium_tables" / f"{r_show}_spots.csv").set_index("barcode").loc[
    spots.index[sel].str.split("_", n=1).str[1]]
xy = sp[["pxl_col_fullres", "pxl_row_fullres"]].to_numpy(float)
zone = spots.major_labl[sel].iloc[0]
pat = spots.patient[sel].iloc[0]


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


def outline(mask):
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


def cfg_label(c):
    return c.replace("rna+", "RNA+").replace("rna", "RNA").replace("tf", "TF").replace("ccc", "CCC").replace("_", " w ")


n = len(members)
SHOW_A = members.groupby("config", group_keys=False).apply(
    lambda g: g.iloc[np.linspace(0, len(g) - 1, min(len(g), max(1, round(18 * len(g) / n)))).astype(int)])
SHOW_A = SHOW_A.head(18) if len(SHOW_A) > 18 else SHOW_A
pos_in_members = {ix: k for k, ix in enumerate(members.index)}
ncol = 6
nrow_a = int(np.ceil(n / ncol))
nrow_a = int(np.ceil(len(SHOW_A) / ncol))
fig = plt.figure(figsize=(183 * MM, (30 * nrow_a + 62 + 52) * MM))
gs = fig.add_gridspec(3, 1, height_ratios=[nrow_a, 2.0, 1.7], hspace=0.3)
ga = gs[0].subgridspec(nrow_a, ncol, wspace=0.08, hspace=0.45)
for kk, (i, r) in enumerate(SHOW_A.iterrows()):
    k = pos_in_members[i]
    ax = fig.add_subplot(ga[kk // ncol, kk % ncol])
    m = inst[k][sel]
    draw(ax, [PINK if v else "#E6E6E6" for v in m])
    ax.set_title(f"#{k + 1} {cfg_label(r.config)}\nK {r.target_k}, space {r.spatial_weight}, seed {r.seed}", fontsize=FS)
    if kk == 0:
        ax.text(-0.25, 1.45, "a", transform=ax.transAxes, fontweight="bold", fontsize=FL)
        ax.text(0.0, 1.45, f"Family {FAM}: {n} domain instances" + ("" if len(SHOW_A) == n else
                f" ({len(SHOW_A)} shown, sampled across configurations)"), transform=ax.transAxes,
                fontweight="bold", fontsize=FH, va="bottom")
gb = gs[1].subgridspec(1, 4, width_ratios=[1, 1, 1, 1.5], wspace=0.25)
ax = fig.add_subplot(gb[0])
img = np.asarray(Image.open(DATA / "visium_tables" / f"{r_show}_lowres.png").convert("RGB"))
sf = json.loads((DATA / "visium_tables" / f"{r_show}_scalefactors.json").read_text())["tissue_lowres_scalef"]
ax.imshow(img)
ax.set_xlim(xy[:, 0].min() * sf - 8, xy[:, 0].max() * sf + 8)
ax.set_ylim(xy[:, 1].max() * sf + 8, xy[:, 1].min() * sf - 8)
ax.axis("off")
ax.add_collection(LineCollection(outline(consensus[sel]) * sf, colors="black", linewidths=0.4))
ax.set_title(f"H&E, {zone} {pat}\n(consensus outlined)", fontweight="bold", fontsize=FH)
ax.text(-0.15, 1.12, "b", transform=ax.transAxes, fontweight="bold", fontsize=FL)
ax = fig.add_subplot(gb[1])
nrm = mpl.colors.Normalize(0, 1)
draw(ax, mpl.colormaps["Purples"](nrm(freq[sel])))
ax.set_title("Share of instances\ncontaining each spot", fontweight="bold", fontsize=FH)
cax = ax.inset_axes([0.1, -0.08, 0.8, 0.04])
cb = fig.colorbar(mpl.cm.ScalarMappable(norm=nrm, cmap="Purples"), cax=cax, orientation="horizontal")
cb.ax.tick_params(labelsize=FS, length=1.5, width=0.4)
ax = fig.add_subplot(gb[2])
draw(ax, [PINK if v else "#E6E6E6" for v in consensus[sel]])
ax.set_title(f"Consensus region\n(spots in ≥ 50 % of instances)", fontweight="bold", fontsize=FH)
ax = fig.add_subplot(gb[3])
J = jac[np.ix_(idx_m, idx_m)]
ax.imshow(J, cmap="Greys", vmin=0, vmax=1)
for (a_, b_), v in np.ndenumerate(J):
    if a_ != b_ and n <= 15:
        ax.text(b_, a_, f"{v:.2f}", ha="center", va="center", fontsize=FS - 1.5, color="white" if v > 0.6 else "black")
labs = [f"{k + 1}" for k in range(n)]
step = 1 if n <= 15 else int(np.ceil(n / 15))
ax.set_xticks(range(0, n, step), labs[::step], fontsize=FS, rotation=90)
ax.set_yticks(range(0, n, step), labs[::step], fontsize=FS)
ax.set_title("Jaccard between instances\n(all 27 slides; ≥ 0.5 = linked)", fontweight="bold", fontsize=FH)
ax.tick_params(length=1.5, width=0.4)
ax.text(-0.18, 1.12, "c", transform=ax.transAxes, fontweight="bold", fontsize=FL)
for k in range(n):  # number the instance titles to match the matrix
    pass
# d: the authors' annotations
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))
from plot_patient_region_umap import PATIENT_COLOURS  # noqa: E402
STALLION = list(PATIENT_COLOURS.values())
CT = sorted(spots.celltype_niche.unique(), key=lambda x: int(x.split("_")[-1]))
MOL = sorted(spots.molecular_niche.unique(), key=lambda x: int(x.split("_")[-1]))
ct_col, mol_col = dict(zip(CT, STALLION)), dict(zip(MOL, STALLION))
gd = gs[2].subgridspec(1, 4, width_ratios=[1, 1, 1.25, 1.25], wspace=0.35)
segs = outline(consensus[sel])
for j, (col, cmap, title) in enumerate([("celltype_niche", ct_col, "Authors' cell-type\nniches"),
                                         ("molecular_niche", mol_col, "Authors' molecular\nniches")]):
    ax = fig.add_subplot(gd[j])
    draw(ax, [cmap[v] for v in spots[col].to_numpy()[sel]])
    ax.add_collection(LineCollection(segs, colors="black", linewidths=0.4))
    ax.set_title(title + " (consensus outlined)", fontweight="bold", fontsize=FH)
    if j == 0:
        ax.text(-0.15, 1.12, "d", transform=ax.transAxes, fontweight="bold", fontsize=FL)
for j, (col, cmap, order, title) in enumerate([("celltype_niche", ct_col, CT, "cell-type niche"),
                                                ("molecular_niche", mol_col, MOL, "molecular niche")], start=2):
    ax = fig.add_subplot(gd[j])
    lab_all = spots[col].to_numpy()
    share_of_region = pd.Series(lab_all[consensus]).value_counts(normalize=True).reindex(order, fill_value=0)
    niche_covered = pd.Series(consensus).groupby(lab_all).mean().reindex(order, fill_value=0)
    y = np.arange(len(order))
    ax.barh(y, share_of_region.values, color=[cmap[o] for o in order], height=0.7)
    for yy, (a1, a2) in enumerate(zip(share_of_region.values, niche_covered.values)):
        if a1 >= 0.02:
            ax.text(a1 + 0.01, yy, f"{a1:.0%} (= {a2:.0%} of niche)", va="center", fontsize=FS - 0.5)
    ax.set_yticks(y, [o.split("_")[-1] for o in order], fontsize=FS)
    ax.set_ylim(len(order) - 0.5, -0.5)
    ax.set_xlim(0, max(0.6, share_of_region.max() * 1.6))
    ax.set_xlabel("share of the region's spots (27 slides)", fontsize=FS)
    ax.set_title(f"Region split across the authors'\n{title}s", fontweight="bold", fontsize=FH)
    ax.tick_params(labelsize=FS, length=1.5, width=0.4)
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)
fig.savefig(OUT / f"fig_family_{FAM}_grouping.png", dpi=250, bbox_inches="tight")
fig.savefig(OUT / f"fig_family_{FAM}_grouping.pdf", dpi=300, bbox_inches="tight")
print(members[["config", "target_k", "spatial_weight", "seed", "domain", "n_spots", "fold_vs_ctrl", "j_rna"]].to_string())
print("consensus spots", int(consensus.sum()), "shown on", r_show)
for col in ["celltype_niche", "molecular_niche"]:
    v = pd.Series(spots[col].to_numpy()[consensus]).value_counts(normalize=True).head(3)
    print(col, ", ".join(f"{k} {x:.0%}" for k, x in v.items()))
