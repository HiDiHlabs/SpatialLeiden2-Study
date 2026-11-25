"""Kuppe et al. 2022 Figure 3 (molecular niches) redrawn for the JOINT analysis of all 27 slides, with the
authors' molecular niches next to our joint SpatialLeiden domains of RNA + TF + space and RNA + CCC + space.
Nothing from the published figure is reused; everything is recomputed from the spot tables and our runs.

Domains (seed 0, K = 9): RNA + TF + space (exp3_joint_rna_tf, r0.75_s0), RNA + CCC + space
(exp3_joint_rna_ccc_ours, r0.25_s0). Colours: molecular niches in the authors' discrete palette (ArchR
"stallion", niche order); each of our domains takes the colour of the molecular niche it overlaps most
(Hungarian matching over all 88,704 spots), so equal colours = matching groups; domains without a match are
grey; the verified stressed-cardiomyocyte domain (RNA + TF D5, RNA + CCC D6; D numbered by size) is pink.
Panels (as Fig. 3a-f of the paper):
  a  UMAP of the RNA Harmony embedding (07_integration/fig_joint_umaps.py cache), coloured by molecular niche
     and by each clustering
  b  spatial maps on three slides: control P7 (ACH003), border zone P3 (10X0026), ischaemic P9 rep2 (ACH0013;
     arrow: molecular niche 11 = the artery, Fig. 3b of the paper)
  c  scaled median cell2location composition per niche / domain (z across rows); * = higher than in the
     other niches / domains (one-sided Mann-Whitney on spots, BH adj. P < 0.05, as in the paper)
  d  mean share per patient group (myogenic, ischaemic, fibrotic; slide = unit)
  e  control (ACH003) and border zone (10X0026): myogenic molecular niches 1, 2, 4 vs our matched domains,
     MYBPC3 and ANKRD2 (log1p CP10k)
Also writes agreement (ARI) and the niche-to-domain matching.

Usage (from Derived_Modalities/, env spatialleiden_joint):
    python 07_integration/fig_joint_kuppe_fig3.py data/kuppe_2022 results/kuppe_2022
Output: results/kuppe_2022/joint_kuppe_fig3/fig_joint_kuppe_fig3.png/.pdf, niche_domain_matching.csv
"""

import json
import sys
from pathlib import Path

import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import EllipseCollection
from matplotlib.patches import Patch
from PIL import Image
from scipy.optimize import linear_sum_assignment
from scipy.spatial import cKDTree
from scipy.stats import mannwhitneyu
from sklearn.metrics import adjusted_rand_score
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))
from plot_patient_region_umap import PATIENT_COLOURS  # noqa: E402

DATA, RES = Path(sys.argv[1]), Path(sys.argv[2])
OUT = RES / "joint_kuppe_fig3"
OUT.mkdir(parents=True, exist_ok=True)
MM = 1 / 25.4
FS, FH, FL = 5, 6, 8
mpl.rcParams.update({"font.family": "Arial", "font.size": FS, "axes.titlesize": FH, "pdf.fonttype": 42})
STALLION = list(PATIENT_COLOURS.values())
PINK, GREYS = "#e43fb1", ["#9E9E9E", "#BDBDBD", "#757575", "#D5D5D5"]
CLUST = {"RNA + TF + space": ("exp3_joint_rna_tf", "r0.75_s0", "D5"),
         "RNA + CCC + space": ("exp3_joint_rna_ccc_ours", "r0.25_s0", "D6")}
SLIDES = {"ACH003": "Control P7", "10X0026": "BZ P3", "ACH0013": "IZ P9"}
H5AD = {"ACH003": None, "10X0026": None, "ACH0013": None}
GROUPS = {"group_1": "myogenic", "group_2": "ischaemic", "group_3": "fibrotic"}

e = np.load(RES / "joint_umaps" / "umaps.npz", allow_pickle=True)
names = pd.Index(e["obs_names"].astype(str))
umap = e["RNA"]
run = names.str.split("_", n=1).str[0]
sp = []
for r in pd.unique(run):
    t = pd.read_csv(DATA / "visium_tables" / f"{r}_spots.csv").set_index("barcode")
    t.index = r + "_" + t.index
    sp.append(t)
spots = pd.concat(sp).loc[names]
MOL = sorted(spots.molecular_niche.unique(), key=lambda x: int(x.split("_")[-1]))
mol_col = dict(zip(MOL, STALLION))
mol = spots.molecular_niche.to_numpy()


def load(folder, col):
    lab = pd.read_csv(RES / folder / "exp2_labels.csv.gz", index_col=0, dtype=str)[col].loc[names]
    order = lab.value_counts().index
    return lab.map({d: f"D{i + 1}" for i, d in enumerate(order)}).to_numpy()


labs, cols, match, rows = {}, {}, {}, []
for k, (f, c, stressed) in CLUST.items():
    lab = load(f, c)
    tab = pd.crosstab(pd.Series(lab, name="domain"), pd.Series(mol, name="niche"))
    ri, ci = linear_sum_assignment(-tab.to_numpy())
    m = {tab.index[i]: tab.columns[j] for i, j in zip(ri, ci)}
    col, g = {}, 0
    for d in sorted(tab.index, key=lambda x: int(x[1:])):
        if d == stressed:
            col[d] = PINK
        elif d in m:
            col[d] = mol_col[m[d]]
        else:
            col[d], g = GREYS[g % len(GREYS)], g + 1
        inter = tab.loc[d, m[d]] if d in m else 0
        rows.append({"clustering": k, "domain": d, "n_spots": int(tab.loc[d].sum()), "matched_niche": m.get(d, ""),
                     "overlap_spots": int(inter), "stressed_CM": d == stressed,
                     "jaccard_with_match": inter / (tab.loc[d].sum() + tab[m[d]].sum() - inter) if d in m else np.nan})
    labs[k], cols[k], match[k] = lab, col, m
ari = {k: adjusted_rand_score(mol, labs[k]) for k in labs}
ari_ct = {k: adjusted_rand_score(spots.celltype_niche, labs[k]) for k in labs}
pd.DataFrame(rows).to_csv(OUT / "niche_domain_matching.csv", index=False)
print({k: round(v, 3) for k, v in ari.items()}, "(ARI vs molecular niches)",
      {k: round(v, 3) for k, v in ari_ct.items()}, "(vs cell-type niches)")

props = [c for c in spots.columns if c.startswith("prop_")]
ctn = [p.replace("prop_", "") for p in props]


def comp_table(groups, order):
    P = spots[props]
    med = P.groupby(groups).median().reindex(order)
    z = (med - med.mean()) / med.std(ddof=0).replace(0, 1)
    pv = np.ones(med.shape)
    for i, gname in enumerate(order):
        inn = np.asarray(groups) == gname
        for j, c in enumerate(props):
            pv[i, j] = mannwhitneyu(P[c].to_numpy()[inn], P[c].to_numpy()[~inn], alternative="greater").pvalue
    adj = multipletests(pv.ravel(), method="fdr_bh")[1].reshape(pv.shape)
    return z, adj


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


def umap_arrows(ax):
    for s_ in ax.spines.values():
        s_.set_visible(False)
    ax.set(xticks=[], yticks=[])
    arrow = dict(arrowstyle="-|>", lw=0.6, color="black", mutation_scale=5, shrinkA=0, shrinkB=0)
    ax.annotate("", xy=(0.25, 0.0), xytext=(0.0, 0.0), xycoords="axes fraction", arrowprops=arrow)
    ax.annotate("", xy=(0.0, 0.25), xytext=(0.0, 0.0), xycoords="axes fraction", arrowprops=arrow)
    ax.text(0.0, -0.025, "UMAP 1", transform=ax.transAxes, ha="left", va="top", fontsize=FS)
    ax.text(-0.025, 0.0, "UMAP 2", transform=ax.transAxes, ha="left", va="bottom", rotation=90, rotation_mode="anchor",
            fontsize=FS)


def letter(ax, s, x=-0.15, y=1.05):
    ax.text(x, y, s, transform=ax.transAxes, fontweight="bold", fontsize=FL, va="bottom", ha="left")


fig = plt.figure(figsize=(183 * MM, 245 * MM))
gs = fig.add_gridspec(4, 1, height_ratios=[1.0, 1.55, 1.15, 1.2], hspace=0.3)
COLS3 = [("Molecular niches\n(authors)", mol, mol_col)] + [(k, labs[k], cols[k]) for k in labs]
# a: UMAPs
ga = gs[0].subgridspec(1, 3, wspace=0.15)
rng = np.random.default_rng(0)
order = rng.permutation(len(names))
for j, (title, lab, cmap) in enumerate(COLS3):
    ax = fig.add_subplot(ga[j])
    ax.scatter(umap[order, 0], umap[order, 1], c=pd.Series(lab).map(cmap).to_numpy()[order], s=0.1, linewidths=0,
               rasterized=True)
    if j == 0:
        umap_arrows(ax)
        letter(ax, "a")
    else:
        ax.axis("off")
        ax.text(0.98, 0.02, f"ARI vs niches {ari[title]:.2f}", transform=ax.transAxes, ha="right", fontsize=FS)
    ax.set_title(title, fontweight="bold", fontsize=FH)
# b: spatial maps (rows = slides, columns = niches / clusterings)
gbm = gs[1].subgridspec(len(SLIDES), 4, wspace=0.03, hspace=0.06)
gbb = [[gbm[i, j] for i in range(len(SLIDES))] for j in range(4)]
for i, (r, label) in enumerate(SLIDES.items()):
    m = (run == r)
    xy = spots.loc[m, ["pxl_col_fullres", "pxl_row_fullres"]].to_numpy(float)
    ax = fig.add_subplot(gbb[0][i])
    img = np.asarray(Image.open(DATA / "visium_tables" / f"{r}_lowres.png").convert("RGB"))
    f = json.loads((DATA / "visium_tables" / f"{r}_scalefactors.json").read_text())["tissue_lowres_scalef"]
    ax.imshow(img)
    ax.set_xlim(xy[:, 0].min() * f - 8, xy[:, 0].max() * f + 8)
    ax.set_ylim(xy[:, 1].max() * f + 8, xy[:, 1].min() * f - 8)
    ax.axis("off")
    ax.text(-0.06, 0.5, label, transform=ax.transAxes, rotation=90, ha="right", va="center", fontsize=FS, fontweight="bold")
    if i == 0:
        ax.set_title("H&E", fontweight="bold", fontsize=FH)
        letter(ax, "b", x=-0.35)
    for j, (title, lab, cmap) in enumerate(COLS3, start=1):
        ax = fig.add_subplot(gbb[j][i])
        draw(ax, xy, pd.Series(np.asarray(lab)[m]).map(cmap).to_numpy())
        if i == 0:
            ax.set_title(title.replace("\n", " "), fontweight="bold", fontsize=FH)
        if r == "ACH0013":  # arrow to the artery (molecular niche 11)
            m11 = (mol[m] == "molniche_11")
            if m11.any():
                tip = np.median(xy[m11], axis=0)
                ax.annotate("", xy=tip, xytext=tip + np.array([-2500, -1800]),
                            arrowprops=dict(arrowstyle="-|>", color="black", lw=0.6, mutation_scale=5))
# c: composition heatmaps
gc = gs[2].subgridspec(1, 4, width_ratios=[1, 1, 1, 0.08], wspace=0.55)
norm = mpl.colors.TwoSlopeNorm(vmin=-2.5, vcenter=0, vmax=2.5)
for j, (title, lab, cmap) in enumerate(COLS3):
    order_g = MOL if j == 0 else sorted(set(lab), key=lambda x: int(x[1:]))
    order_g = [g for g in order_g if (np.asarray(lab) == g).sum() >= 50]
    z, adj = comp_table(np.asarray(lab), order_g)
    ax = fig.add_subplot(gc[j])
    ax.imshow(z.to_numpy().clip(-2.5, 2.5), cmap="RdBu_r", norm=norm, aspect="auto")
    for (a_, b_), v in np.ndenumerate(adj):
        if v < 0.05 and z.to_numpy()[a_, b_] > 0:
            ax.text(b_, a_, "*", ha="center", va="center", fontsize=FS)
    ax.set_xticks(range(len(ctn)), ctn, rotation=90, fontsize=FS)
    yl = [g.split("_")[-1] if j == 0 else g for g in order_g]
    ax.set_yticks(range(len(order_g)), yl, fontsize=FS)
    for t_, g in zip(ax.get_yticklabels(), order_g):
        t_.set_color(cmap[g] if cmap[g] not in GREYS else "#555555")
        if cmap[g] == PINK:
            t_.set_fontweight("bold")
    ax.tick_params(length=1.5, width=0.4, pad=1)
    ax.set_title(title.replace("\n", " "), fontweight="bold", fontsize=FH)
    if j == 0:
        letter(ax, "c", x=-0.3)
cax = fig.add_subplot(gc[3])
cb = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap="RdBu_r"), cax=cax)
cb.ax.tick_params(labelsize=FS, length=1.5, width=0.4)
cb.set_label("scaled median composition", fontsize=FS)
cb.outline.set_linewidth(0.4)
# d + e
gde = gs[3].subgridspec(1, 2, width_ratios=[0.9, 2.1], wspace=0.15)
ax = fig.add_subplot(gde[0])
grp = spots.patient_group.map(GROUPS).to_numpy()
x, ticks, tl = 0, [], []
for title, lab, cmap in COLS3:
    df = pd.DataFrame({"run": run, "grp": grp, "lab": np.asarray(lab)})
    share = pd.crosstab([df.run, df.grp], df.lab, normalize="index")
    mean = share.groupby(level="grp").mean().reindex(["myogenic", "ischaemic", "fibrotic"])
    keys = MOL if title.startswith("Molecular") else sorted(mean.columns, key=lambda q: int(q[1:]))
    for gname in mean.index:
        bottom = 0
        for kk in keys:
            v = mean.loc[gname, kk] if kk in mean.columns else 0
            ax.bar(x, v, bottom=bottom, color=cmap[kk], width=0.8, linewidth=0)
            bottom += v
        ticks.append(x)
        tl.append(gname)
        x += 1
    x += 0.6
ax.set_xticks(ticks, tl, rotation=90, fontsize=FS)
ax.set_ylim(0, 1)
ax.set_ylabel("mean share per slide", fontsize=FS)
ax.tick_params(labelsize=FS, length=1.5, width=0.4)
for s_ in ("top", "right"):
    ax.spines[s_].set_visible(False)
for xx, t in zip([1, 4.6, 8.2], ["niches", "RNA + TF", "RNA + CCC"]):
    ax.text(xx, 1.03, t, ha="center", fontsize=FS, fontweight="bold", transform=ax.get_xaxis_transform())
letter(ax, "d", x=-0.35, y=1.1)
ge = gde[1].subgridspec(2, 6, wspace=0.05, hspace=0.08)
MYO = ["molniche_1", "molniche_2", "molniche_4"]
genes = {}
for f in sorted((DATA / "processed1").glob("Visium_*.h5ad")):
    nm = set(ad.read_h5ad(f, backed="r").obs_names)
    for r in ["ACH003", "10X0026"]:
        if set(spots.index[run == r].str.split("_", n=1).str[1]) == nm:
            h = ad.read_h5ad(f)
            tot = np.asarray(h.raw.X.sum(1)).ravel()
            vn = list(h.raw.var_names)
            genes[r] = {g: pd.Series(np.log1p(np.asarray(h.raw.X[:, vn.index(g)].todense()).ravel() / tot * 1e4),
                                     index=r + "_" + h.obs_names) for g in ["MYBPC3", "ANKRD2"]}
gmax = {g: float(np.quantile(np.concatenate([genes[r][g].to_numpy() for r in genes]), 0.99)) for g in ["MYBPC3", "ANKRD2"]}
for i, r in enumerate(["ACH003", "10X0026"]):
    m = (run == r)
    idx = spots.index[m]
    xy = spots.loc[m, ["pxl_col_fullres", "pxl_row_fullres"]].to_numpy(float)
    ax = fig.add_subplot(ge[i, 0])
    img = np.asarray(Image.open(DATA / "visium_tables" / f"{r}_lowres.png").convert("RGB"))
    f = json.loads((DATA / "visium_tables" / f"{r}_scalefactors.json").read_text())["tissue_lowres_scalef"]
    ax.imshow(img)
    ax.set_xlim(xy[:, 0].min() * f - 8, xy[:, 0].max() * f + 8)
    ax.set_ylim(xy[:, 1].max() * f + 8, xy[:, 1].min() * f - 8)
    ax.axis("off")
    ax.text(-0.06, 0.5, SLIDES[r], transform=ax.transAxes, rotation=90, ha="right", va="center", fontsize=FS,
            fontweight="bold")
    if i == 0:
        letter(ax, "e", x=-0.45, y=1.12)
    panels = [("Niches 1, 2, 4", [mol_col[q] if q in MYO else "#EEEEEE" for q in mol[m]])]
    for k in labs:
        mm = {d for d, q in match[k].items() if q in MYO}
        panels.append((f"{k.replace(' + space', '').replace(' + ', '+')}\nmatched", [cols[k][d] if d in mm else "#EEEEEE" for d in labs[k][m]]))
    for j, (t, c) in enumerate(panels, start=1):
        ax = fig.add_subplot(ge[i, j])
        draw(ax, xy, c)
        if i == 0:
            ax.set_title(t, fontweight="bold", fontsize=FH)
    for j, g in enumerate(["MYBPC3", "ANKRD2"], start=4):
        ax = fig.add_subplot(ge[i, j])
        nrm = mpl.colors.Normalize(0, gmax[g], clip=True)
        draw(ax, xy, mpl.colormaps["magma"](nrm(genes[r][g].loc[idx].to_numpy())))
        if i == 0:
            ax.set_title(g, fontweight="bold", fontsize=FH, fontstyle="italic")
        if i == 1:
            cax = ax.inset_axes([0.1, -0.1, 0.8, 0.05])
            cb = fig.colorbar(mpl.cm.ScalarMappable(norm=nrm, cmap="magma"), cax=cax, orientation="horizontal")
            cb.ax.tick_params(labelsize=FS, length=1.5, width=0.4)
            cb.outline.set_linewidth(0.4)
# legend: niches + rule
lg = [Patch(color=mol_col[q], label=q.split("_")[-1]) for q in MOL] + [Patch(color=PINK, label="stressed CM"),
                                                                       Patch(color=GREYS[0], label="no niche match")]
fig.legend(handles=lg, loc="lower center", ncol=len(lg), frameon=False, fontsize=FS, handlelength=0.9,
           columnspacing=0.8, bbox_to_anchor=(0.5, 0.045),
           title="Molecular niche (our domains: colour of the best-matching niche)", title_fontsize=FH)
fig.savefig(OUT / "fig_joint_kuppe_fig3.png", dpi=250, bbox_inches="tight")
fig.savefig(OUT / "fig_joint_kuppe_fig3.pdf", dpi=300, bbox_inches="tight")
print(pd.DataFrame(rows).round(2).to_string(index=False))
