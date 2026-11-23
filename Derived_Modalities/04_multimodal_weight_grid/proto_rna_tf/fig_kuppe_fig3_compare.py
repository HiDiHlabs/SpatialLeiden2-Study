"""Kuppe et al. 2022 Figure 3 (molecular niches), redrawn on our three slides with the authors'
molecular niches next to our SpatialLeiden domains. Nothing from the published figure is
reused; everything is recomputed from the deposited spot tables and our runs.

Ours = RNA + TF + space (mm_proto_rna_tf, sp0.75_tf0.5_s0: RNA 1, TF 0.5, spatial 0.75,
K = authors' cell-type niche count, seed 0). Each of our domains is coloured by the authors'
molecular niche it overlaps most (Hungarian matching per slide), so equal colours = matching
niche; domains without a match are grey shades.

Figure A (fig3_compare_maps): per slide, H&E | authors' molecular niches | our domains |
    UMAP of spots (RNA PCA of that slide) coloured by authors' niches | by our domains.
    (Kuppe Fig. 3a-b; arrows mark molecular niche 11 = MYH11+ vSMCs, the artery on P9 rep2.)
Figure B (fig3_compare_composition): Kuppe Fig. 3c-f analogues.
    c  scaled median cell2location composition per niche/domain (z across rows); * = higher in
       that niche/domain than in the others (one-sided Mann-Whitney, BH adj. P < 0.05)
    d  niche/domain shares per slide
    e,f  P1 control and P9 IZ: H&E, myogenic molecular niches 1/2/4 vs matched domains of ours,
       MYBPC3 and ANKRD2 expression (log1p CP10k). (The paper shows control and border zone; we
       have no border-zone slide here.)

Usage (env spatialleiden, from Derived_Modalities/):
    python 04_multimodal_weight_grid/proto_rna_tf/fig_kuppe_fig3_compare.py data/kuppe_2022 results/kuppe_2022
"""

import importlib.util
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve()
spec = importlib.util.spec_from_file_location(
    "exp0_sweep", HERE.parents[2] / "03_single_modality_clustering" / "exp0_layer_ratio_sweep" / "exp0_sweep.py")
e0 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(e0)

import anndata as ad  # noqa: E402
import matplotlib as mpl  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scanpy as sc  # noqa: E402
from matplotlib.collections import EllipseCollection  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from PIL import Image  # noqa: E402
from scipy.optimize import linear_sum_assignment  # noqa: E402
from scipy.spatial import cKDTree  # noqa: E402
from scipy.stats import mannwhitneyu  # noqa: E402
from statsmodels.stats.multitest import multipletests  # noqa: E402

DATA, RES = Path(sys.argv[1]), Path(sys.argv[2])
OUT = RES / "kuppe_fig3_compare"
OUT.mkdir(parents=True, exist_ok=True)
MM = 1 / 25.4
mpl.rcParams.update({"font.family": "Arial", "font.size": 6, "axes.titlesize": 6.5, "pdf.fonttype": 42})
SLIDES = {"10X001": "P1 control\n(10X001)", "ACH0012": "P9 IZ\n(ACH0012)", "ACH0013": "P9 IZ rep2\n(ACH0013)"}
SLIDE_COL = {"10X001": "#0072B2", "ACH0013": "#D55E00", "ACH0012": "#E69F00"}
H5AD = {"10X001": "Visium_control_P1.h5ad", "ACH0013": "Visium_GT_IZ_P9_rep2.h5ad", "ACH0012": "Visium_GT_IZ_P9.h5ad"}
MOL = [f"molniche_{i}" for i in range(1, 13)]
MOL_COL = dict(zip(MOL, ["#D51F26", "#272E6A", "#208A42", "#89288F", "#F47D2B", "#FEE500", "#8A9FD1",
                         "#C06CAB", "#D8A767", "#90D5E4", "#89C75F", "#F37B7D"]))
GREYS = ["#9E9E9E", "#BDBDBD", "#757575", "#E0E0E0"]
OURS = ("mm_proto_rna_tf", "sp0.75_tf0.5_s0")


def short(n):
    return n.replace("molniche_", "")


def draw(ax, xy, colours, s=1.05):
    pitch = np.median(cKDTree(xy).query(xy, k=2)[0][:, 1])
    d = np.full(len(xy), s * pitch)
    coll = EllipseCollection(d, d, np.zeros(len(xy)), units="xy", offsets=xy, offset_transform=ax.transData,
                             facecolors=colours, linewidths=0)
    coll.set_rasterized(True)
    ax.add_collection(coll)
    ax.set_xlim(xy[:, 0].min() - pitch, xy[:, 0].max() + pitch)
    ax.set_ylim(xy[:, 1].max() + pitch, xy[:, 1].min() - pitch)
    ax.set_aspect("equal")
    ax.axis("off")


def he(ax, run, xy):
    img = np.asarray(Image.open(DATA / "visium_tables" / f"{run}_lowres.png").convert("RGB"))
    f = json.loads((DATA / "visium_tables" / f"{run}_scalefactors.json").read_text())["tissue_lowres_scalef"]
    ax.imshow(img)
    ax.set_xlim(xy[:, 0].min() * f - 8, xy[:, 0].max() * f + 8)
    ax.set_ylim(xy[:, 1].max() * f + 8, xy[:, 1].min() * f - 8)
    ax.axis("off")


# ---- load per slide ----
S = {}
for run in SLIDES:
    rna = e0.load_slide(DATA, run)
    e0.preprocess(rna)
    sc.tl.umap(rna, random_state=0)
    lab = pd.read_csv(RES / OURS[0] / f"labels_{run}.csv.gz", index_col=0, dtype=str)[OURS[1]]
    assert list(lab.index) == list(rna.obs_names)
    sp = pd.read_csv(DATA / "visium_tables" / f"{run}_spots.csv").set_index("barcode").loc[lab.index]
    tab = pd.crosstab(lab, sp.molecular_niche)
    r, c = linear_sum_assignment(-tab.to_numpy())
    match = {tab.index[i]: tab.columns[j] for i, j in zip(r, c)}
    g = 0
    col = {}
    for d in sorted(set(lab), key=int):
        if d in match:
            col[d] = MOL_COL[match[d]]
        else:
            col[d] = GREYS[g % len(GREYS)]
            g += 1
    h = ad.read_h5ad(DATA / "processed1" / H5AD[run])
    tot = np.asarray(h.raw.X.sum(1)).ravel()
    names = list(h.raw.var_names)
    genes = {gn: pd.Series(np.log1p(np.asarray(h.raw.X[:, names.index(gn)].todense()).ravel() / tot * 1e4),
                           index=h.obs_names).loc[lab.index].to_numpy() for gn in ["MYBPC3", "ANKRD2"]}
    S[run] = dict(sp=sp, lab=lab, umap=rna.obsm["X_umap"], match=match, col=col, genes=genes,
                  xy=sp[["pxl_col_fullres", "pxl_row_fullres"]].to_numpy(float))
    print(run, "our domain -> molecular niche:", {d: short(m) for d, m in sorted(match.items(), key=lambda x: int(x[0]))})

# ---- Figure A: maps + UMAPs ----
fig, axes = plt.subplots(3, 5, figsize=(183 * MM, 118 * MM), gridspec_kw={"wspace": 0.06, "hspace": 0.1})
for i, (run, label) in enumerate(SLIDES.items()):
    s = S[run]
    he(axes[i, 0], run, s["xy"])
    draw(axes[i, 1], s["xy"], [MOL_COL[m] for m in s["sp"].molecular_niche])
    draw(axes[i, 2], s["xy"], [s["col"][d] for d in s["lab"]])
    for j, colours in [(3, [MOL_COL[m] for m in s["sp"].molecular_niche]), (4, [s["col"][d] for d in s["lab"]])]:
        ax = axes[i, j]
        ax.scatter(*s["umap"].T, c=colours, s=0.6, linewidths=0, rasterized=True)
        ax.set_xticks([]); ax.set_yticks([])
        for sp_ in ax.spines.values():
            sp_.set_linewidth(0.4)
    axes[i, 0].text(-0.08, 0.5, label, transform=axes[i, 0].transAxes, rotation=90, ha="right", va="center",
                    fontweight="bold", color=SLIDE_COL[run], fontsize=7)
    if run == "ACH0013":  # arrow at the artery (molecular niche 11) on both maps
        m11 = (s["sp"].molecular_niche == "molniche_11").to_numpy()
        tip = np.median(s["xy"][m11], axis=0)
        for j in (1, 2):
            axes[i, j].annotate("", xy=tip, xytext=tip + np.array([-2500, -1800]),
                                arrowprops=dict(arrowstyle="-|>", color="black", lw=0.8, mutation_scale=6))
for j, t in enumerate(["H&E", "Authors'\nmolecular niches", "Ours: RNA + TF\n+ space", "UMAP, authors'\nniches",
                       "UMAP, our\ndomains"]):
    axes[0, j].set_title(t, fontweight="bold")
present = sorted({m for s in S.values() for m in s["sp"].molecular_niche}, key=lambda m: int(short(m)))
fig.legend(handles=[Patch(color=MOL_COL[m], label=short(m)) for m in present] + [Patch(color=GREYS[0], label="ours, no match")],
           loc="lower center", ncol=len(present) + 1, frameon=False, fontsize=5.5, handlelength=1, columnspacing=0.8,
           bbox_to_anchor=(0.5, 0.03), title="Molecular niche (colour of our domains = best-matching niche)", title_fontsize=5.5)
fig.savefig(OUT / "fig3_compare_maps.png", dpi=250, bbox_inches="tight")
fig.savefig(OUT / "fig3_compare_maps.pdf", dpi=300, bbox_inches="tight")
plt.close(fig)


# ---- Figure B: composition heatmaps, shares, myogenic niches + markers ----
def comp_table(groups, P):
    """Median composition per group, z-scored across groups per cell type, with one-sided tests."""
    med = P.groupby(groups).median()
    z = (med - med.mean()) / med.std(ddof=0).replace(0, 1)
    pv = pd.DataFrame(1.0, index=med.index, columns=med.columns)
    for gname in med.index:
        inn = (groups == gname).to_numpy()
        for c in med.columns:
            pv.loc[gname, c] = mannwhitneyu(P[c][inn], P[c][~inn], alternative="greater").pvalue
    adj = multipletests(pv.to_numpy().ravel(), method="fdr_bh")[1].reshape(pv.shape)
    return z, pd.DataFrame(adj, index=pv.index, columns=pv.columns)


props = [c for c in S["10X001"]["sp"].columns if c.startswith("prop_")]
pooled = pd.concat([S[r]["sp"] for r in SLIDES])
za, aa = comp_table(pooled.molecular_niche, pooled[props])
za = za.loc[sorted(za.index, key=lambda m: int(short(m)))]
aa = aa.loc[za.index]
zo, ao, rowcol = [], [], []
for run in SLIDES:
    s = S[run]
    g = s["lab"].map(lambda d: f"{run}:{d}")
    z, a = comp_table(g, s["sp"][props])
    order = sorted(z.index, key=lambda x: int(x.split(":")[1]))
    zo.append(z.loc[order]); ao.append(a.loc[order])
    rowcol += [s["col"][x.split(":")[1]] for x in order]
zo, ao = pd.concat(zo), pd.concat(ao)
ct = [p.replace("prop_", "") for p in props]

fig = plt.figure(figsize=(183 * MM, 170 * MM))
gs = fig.add_gridspec(2, 1, height_ratios=[1.25, 1], hspace=0.36)
top = gs[0].subgridspec(1, 3, width_ratios=[1, 1.15, 0.7], wspace=0.75)
norm = mpl.colors.TwoSlopeNorm(vmin=-2.5, vcenter=0, vmax=2.5)


def heat(ax, z, a, ylabels, title, rowcolours=None):
    ax.imshow(z.to_numpy().clip(-2.5, 2.5), cmap="RdBu_r", norm=norm, aspect="auto")
    for (i, j), v in np.ndenumerate(a.to_numpy()):
        if v < 0.05 and z.to_numpy()[i, j] > 0:
            ax.text(j, i, "*", ha="center", va="center", fontsize=6, color="black")
    ax.set_xticks(range(len(ct)), ct, rotation=90, fontsize=5.5)
    ax.set_yticks(range(len(ylabels)), ylabels, fontsize=5)
    ax.tick_params(length=1.5, width=0.4)
    ax.tick_params(axis="y", pad=7 if rowcolours is not None else 2)
    ax.set_title(title, fontweight="bold")
    if rowcolours is not None:
        for i, c in enumerate(rowcolours):
            ax.add_patch(mpl.patches.Rectangle((-1.15, i - 0.5), 0.55, 1, color=c, clip_on=False))


ax_c1 = fig.add_subplot(top[0])
heat(ax_c1, za, aa, [short(m) for m in za.index], "c  Authors' molecular niches\n(3 slides pooled)",
     [MOL_COL[m] for m in za.index])
ax_c2 = fig.add_subplot(top[1])
ylab = [f"{'P1' if r == '10X001' else ('P9' if r == 'ACH0012' else 'P9r2')} d{d}" for r, d in (x.split(":") for x in zo.index)]
heat(ax_c2, zo, ao, ylab, "c  Our domains\n(per slide)", rowcol)
cax = ax_c2.inset_axes([1.05, 0.0, 0.05, 0.22])
cb = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap="RdBu_r"), cax=cax)
cb.set_label("scaled median\ncomposition", fontsize=5)
cb.ax.tick_params(labelsize=5, length=1.5, width=0.4)
cb.outline.set_linewidth(0.4)
# d: shares per slide, authors (A) vs ours (O)
ax_d = fig.add_subplot(top[2])
x = 0
ticks, tl = [], []
for run in SLIDES:
    s = S[run]
    for kind in ("A", "O"):
        if kind == "A":
            sh = s["sp"].molecular_niche.value_counts(normalize=True)
            cols = [MOL_COL[m] for m in sh.index]
        else:
            sh = s["lab"].value_counts(normalize=True)
            cols = [s["col"][d] for d in sh.index]
        bottom = 0
        for v, c in zip(sh.to_numpy(), cols):
            ax_d.bar(x, v, bottom=bottom, color=c, width=0.8, linewidth=0)
            bottom += v
        ticks.append(x); tl.append(f"{'P1' if run == '10X001' else ('P9' if run == 'ACH0012' else 'P9 r2')} {'authors' if kind == 'A' else 'ours'}")
        x += 1
    x += 0.5
ax_d.set_xticks(ticks, tl, rotation=90, fontsize=5.5)
ax_d.set_ylabel("share of spots", fontsize=6, labelpad=1)
ax_d.set_ylim(0, 1)
ax_d.tick_params(labelsize=5, length=1.5, width=0.4)
for k in ("top", "right"):
    ax_d.spines[k].set_visible(False)
ax_d.set_title("d  Niche shares", fontweight="bold")
# e,f: myogenic niches and markers on P1 control and P9 IZ
bot = gs[1].subgridspec(2, 5, wspace=0.06, hspace=0.12)
MYO = ["molniche_1", "molniche_2", "molniche_4"]
gvmax = {gn: float(np.quantile(np.concatenate([S[r]["genes"][gn] for r in ("10X001", "ACH0012")]), 0.99))
         for gn in ["MYBPC3", "ANKRD2"]}
for i, run in enumerate(["10X001", "ACH0012"]):
    s = S[run]
    axs = [fig.add_subplot(bot[i, j]) for j in range(5)]
    he(axs[0], run, s["xy"])
    draw(axs[1], s["xy"], [MOL_COL[m] if m in MYO else "#EEEEEE" for m in s["sp"].molecular_niche])
    draw(axs[2], s["xy"], [s["col"][d] if s["match"].get(d) in MYO else "#EEEEEE" for d in s["lab"]])
    for j, gn in [(3, "MYBPC3"), (4, "ANKRD2")]:
        nrm = mpl.colors.Normalize(0, gvmax[gn], clip=True)
        draw(axs[j], s["xy"], mpl.colormaps["magma"](nrm(s["genes"][gn])))
        if i == 1:
            cax = axs[j].inset_axes([0.1, -0.1, 0.8, 0.05])
            cb = fig.colorbar(mpl.cm.ScalarMappable(norm=nrm, cmap="magma"), cax=cax, orientation="horizontal")
            cb.ax.tick_params(labelsize=5, length=1.5, width=0.4)
            cb.set_label("log1p(CP10k)", fontsize=5)
            cb.outline.set_linewidth(0.4)
    axs[0].text(-0.08, 0.5, SLIDES[run], transform=axs[0].transAxes, rotation=90, ha="right", va="center",
                fontweight="bold", color=SLIDE_COL[run], fontsize=6.5)
    if i == 0:
        for ax, t in zip(axs, ["e,f  H&E", "Authors' niches\n1, 2, 4 (myogenic)", "Ours, matched\nto 1, 2, 4",
                               "MYBPC3", "ANKRD2"]):
            ax.set_title(t, fontweight="bold")
fig.savefig(OUT / "fig3_compare_composition.png", dpi=250, bbox_inches="tight")
fig.savefig(OUT / "fig3_compare_composition.pdf", dpi=300, bbox_inches="tight")
pd.DataFrame([{"slide": r, "our_domain": d, "matched_molecular_niche": m, "n_spots": int((S[r]["lab"] == d).sum()),
               "overlap": int(((S[r]["lab"] == d) & (S[r]["sp"].molecular_niche == m)).sum())}
              for r in SLIDES for d, m in S[r]["match"].items()]).to_csv(OUT / "domain_to_molecular_niche.csv", index=False)
