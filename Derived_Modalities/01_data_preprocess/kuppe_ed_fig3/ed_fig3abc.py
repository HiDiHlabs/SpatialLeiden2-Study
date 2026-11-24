"""Reproduce Kuppe et al. 2022, Extended Data Fig. 3a-c (cell-type niches).

a, b: UMAP of the Visium spots of all 27 sections built from their cell-type compositions,
      coloured by patient (a) and region (b).
c:    cell-type niche proportions per sample (the authors' celltype_niche labels).

Paper Methods ("Niche definitions from spatial transcriptomics data"): the cell2location
cell-type proportions of each spot are transformed into isometric log ratios (ILR) and spots
are clustered on them. No Harmony / batch correction is involved: compositions are compared
directly across sections. The UMAP parameters are not stated; here: ILR (11 cell types -> 10
coordinates, Helmert basis), 15-NN euclidean graph, scanpy UMAP, seed 0.

Input: data/kuppe_2022/visium_tables/<run>_spots.csv (prop_* = cell2location proportions,
celltype_niche, patient, major_labl, patient_region_id), runs = the 27 of the paper (no ACH0021).

Usage (from Derived_Modalities/):
    python 01_data_preprocess/kuppe_ed_fig3/ed_fig3abc.py data/kuppe_2022 results/kuppe_2022/kuppe_ed_fig3
"""
import sys
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.linalg import helmert

SEED = 0
# x-axis order of the paper's Extended Data Fig. 3c (patient_region_id, as printed)
PAPER_ORDER = ["control_P1", "control_P17", "control_P7", "control_P8", "FZ_P14", "FZ_P18", "FZ_P20",
               "FZ/GT_P19", "FZ/GT_P4", "GT/IZ_P13", "GT/IZ_P15", "GT/IZ_P9", "GT/IZ_P9_rep2", "IZ_P10",
               "IZ_P15", "IZ_P16", "IZ_P3", "IZ/BZ_P2", "RZ_P11", "RZ_P3", "RZ_P6", "RZ_P9", "RZ/BZ_P12",
               "RZ/BZ_P2", "RZ/BZ_P3", "RZ/FZ_P5", "RZ/GT_P2"]
data, out = Path(sys.argv[1]), Path(sys.argv[2])
out.mkdir(parents=True, exist_ok=True)

meta = pd.read_csv(data / "visium_raw" / "metadata-Visium.csv")
runs = [r for r in meta.hca_sample_id if r != "ACH0021"]
s = pd.concat([pd.read_csv(data / "visium_tables" / f"{r}_spots.csv") for r in runs], ignore_index=True)
P = [c for c in s.columns if c.startswith("prop_")]
assert len(runs) == 27 and (s[P] > 0).all().all()

# ILR: clr, then project on an orthonormal basis of the clr hyperplane (Helmert)
x = np.log(s[P].to_numpy())
clr = x - x.mean(1, keepdims=True)
ilr = clr @ helmert(len(P)).T  # (n, D-1)

a = ad.AnnData(obs=s[["patient", "major_labl", "patient_region_id", "celltype_niche", "hca_sample_id"]].astype(str))
a.obsm["X_ilr"] = ilr
sc.pp.neighbors(a, n_neighbors=15, use_rep="X_ilr", random_state=SEED)
sc.tl.umap(a, random_state=SEED)
u = a.obsm["X_umap"]
pd.DataFrame(u, columns=["UMAP1", "UMAP2"], index=s.hca_sample_id + "_" + s.barcode).join(
    a.obs.set_index(s.hca_sample_id + "_" + s.barcode)).to_csv(out / "ed_fig3_umap.csv.gz")

# ---- figure -------------------------------------------------------------------------------

def umap_arrows(ax, title):
    """Paper style: no frame or ticks, two short arrows labelled UMAP 1 / UMAP 2 at the bottom left."""
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set(xticks=[], yticks=[], xlabel="", ylabel="")
    ax.set_title(title, fontsize=12)
    kw = dict(xycoords="axes fraction", arrowprops=dict(arrowstyle="-|>", lw=1.3, color="black", mutation_scale=12, shrinkA=0, shrinkB=0, joinstyle="miter", capstyle="projecting"))
    ax.annotate("", xy=(0.15, 0.0), xytext=(0.0, 0.0), **kw)
    ax.annotate("", xy=(0.0, 0.15), xytext=(0.0, 0.0), **kw)
    ax.text(0.075, -0.035, "UMAP 1", transform=ax.transAxes, ha="center", va="top", fontsize=10)
    ax.text(-0.035, 0.075, "UMAP 2", transform=ax.transAxes, ha="right", va="center", rotation=90, fontsize=10)


fig = plt.figure(figsize=(17, 5.6))
gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.7])
order = np.random.default_rng(SEED).permutation(len(s))
pats = sorted(a.obs.patient.unique(), key=lambda p: int(p[1:]))
# patient colours: ArchR::paletteDiscrete default ("stallion"), used for all discrete UMAPs in the
# authors' plot_figures notebooks; order as in the paper's Extended Data Fig. 3a legend
STALLION = ["#d51f26", "#272e6a", "#208a42", "#89288f", "#f47d2b", "#fee500", "#8a9fd1", "#c06cab",
            "#d8a767", "#90d5e4", "#89c75f", "#f37b7d", "#9983bd", "#d24b27", "#3bbca8", "#6e4b9e",
            "#0c727c", "#7e1416", "#e6c2dc", "#3d3d3d"]
assert pats == [f"P{i}" for i in range(1, 21)]
pcol = dict(zip(pats, STALLION))
# region colours: cols.region in the authors' plot_figures/fig4.ipynb (saezlab/visium_heart);
# niche colours: niche_cols in st_snRNAseq/05.1_nicheanalysis/plot_compniche_bars.R
zcol = {"control": "#4CAF50", "FZ": "#009DB2", "IZ": "#B50084", "RZ": "#FFC107", "BZ": "#F57F17"}
zlab = a.obs.major_labl.replace({"CTRL": "control"})

ax = fig.add_subplot(gs[0])
ax.scatter(*u[order].T, c=[pcol[p] for p in a.obs.patient.values[order]], s=0.3, rasterized=True)
ax.legend(handles=[plt.Line2D([], [], ls="", marker="o", markersize=6, color=pcol[p], label=p) for p in pats],
          fontsize=9, markerscale=1, handletextpad=0.3, labelspacing=0.35, loc="center left", bbox_to_anchor=(1, 0.5), frameon=False)
umap_arrows(ax, "patients")
ax.text(-0.08, 1.04, "a", transform=ax.transAxes, fontsize=15, fontweight="bold")

ax = fig.add_subplot(gs[1])
ax.scatter(*u[order].T, c=[zcol[z] for z in zlab.values[order]], s=0.3, rasterized=True)
ax.legend(handles=[plt.Line2D([], [], ls="", marker="o", markersize=6, color=c, label=z) for z, c in zcol.items()],
          fontsize=10, markerscale=1, handletextpad=0.3, labelspacing=0.35, loc="lower left", bbox_to_anchor=(0.92, 0.02), frameon=False)
umap_arrows(ax, "regions")
ax.text(-0.08, 1.04, "b", transform=ax.transAxes, fontsize=15, fontweight="bold")

ax = fig.add_subplot(gs[2])
niches = sorted(a.obs.celltype_niche.unique(), key=lambda n: int(n.split("_")[-1]))
ct = pd.crosstab(s.patient_region_id, a.obs.celltype_niche.values, normalize="index")[niches]
assert set(ct.index) == set(PAPER_ORDER), set(ct.index) ^ set(PAPER_ORDER)
ct = ct.loc[PAPER_ORDER]
ncol = ["#d51f26", "#272e6a", "#208a42", "#89288f", "#f47d2b", "#fee500", "#8a9fd1", "#c06cab", "#d8a767"]
bottom = np.zeros(len(ct))
for n, c in reversed(list(zip(niches, ncol))):  # niche_9 at the bottom, niche_1 on top, as in the paper
    ax.bar(range(len(ct)), ct[n], bottom=bottom, color=c, width=0.9, edgecolor="white", linewidth=0.3, label=n.replace("ctniche", "niche"))
    bottom += ct[n].to_numpy()
h, l = ax.get_legend_handles_labels()
ax.legend(h[::-1], l[::-1], fontsize=10, handlelength=1.6, handleheight=1.6, labelspacing=0.15, loc="center left", bbox_to_anchor=(1, 0.5), frameon=False)
ax.set_xticks(range(len(ct)), ct.index, rotation=45, ha="right", rotation_mode="anchor", fontsize=9)
ax.set(ylabel="cell type proportion", ylim=(0, 1), xlim=(-0.6, len(ct) - 0.4), yticks=[0, 0.25, 0.5, 0.75, 1])
ax.spines[["top", "right"]].set_visible(False)  # theme_classic, as in the paper
ax.spines[["left", "bottom"]].set_linewidth(2)
ax.tick_params(width=1.5, length=5, labelsize=9)
ax.text(-0.06, 1.04, "c", transform=ax.transAxes, fontsize=15, fontweight="bold")
fig.tight_layout()
fig.savefig(out / "ed_fig3abc.png", dpi=150)
fig.savefig(out / "ed_fig3abc.pdf")
ct.to_csv(out / "ed_fig3c_niche_proportions.csv")
print("done:", out)
