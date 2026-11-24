"""Reproduce Kuppe et al. 2022, Extended Data Fig. 3d (feature plots on the composition UMAP).

Top row: cell2location proportions of cardiomyocytes, fibroblasts, cycling cells, pericytes,
adipocytes, endothelial cells, myeloid cells and vSMCs on the composition UMAP of ed_fig3abc.py.
Bottom row: marker genes RYR2, PDGFRA, MKI67, ABCC9, FASN, VWF, IL7R, MYH11 as a gene-weighted
kernel density on the UMAP. The paper uses the R package Nebulosa; here an approximation:
log-normalised expression (CP10k, log1p, per section from the raw counts) as weights of a 2D
histogram on a 400 x 400 grid, Gaussian-smoothed (sigma = 4 grid cells), read back at each spot.
A gene absent from a section's h5ad (MKI67 in 5 sections, IL7R in 3) is taken as zero there:
the processed files leave out genes without counts.

Usage (from Derived_Modalities/, after ed_fig3abc.py):
    python 01_data_preprocess/kuppe_ed_fig3/ed_fig3d.py data/kuppe_2022 results/kuppe_2022/kuppe_ed_fig3
"""
import sys
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter

data, out = Path(sys.argv[1]), Path(sys.argv[2])
CT = [("cardiomyocytes", "prop_CM"), ("fibroblasts", "prop_Fib"), ("cycling cells", "prop_prolif"),
      ("pericytes", "prop_PC"), ("adipocytes", "prop_Adipo"), ("endothelial cells", "prop_Endo"),
      ("myeloid cells", "prop_Myeloid"), ("vSMCs", "prop_vSMCs")]
GENES = ["RYR2", "PDGFRA", "MKI67", "ABCC9", "FASN", "VWF", "IL7R", "MYH11"]

um = pd.read_csv(out / "ed_fig3_umap.csv.gz", index_col=0)
meta = pd.read_csv(data / "visium_raw" / "metadata-Visium.csv")
runs = [r for r in meta.hca_sample_id if r != "ACH0021"]
spots = pd.concat([pd.read_csv(data / "visium_tables" / f"{r}_spots.csv").assign(run=r) for r in runs], ignore_index=True)
spots.index = spots.run + "_" + spots.barcode
spots = spots.loc[um.index]

# log-normalised expression of the 8 genes, per section from the raw counts
expr = pd.DataFrame(0.0, index=um.index, columns=GENES)
for f in sorted((data / "processed1").glob("Visium_*.h5ad")):
    a = ad.read_h5ad(f)
    raw = a.raw.to_adata()
    run = [r for r in runs if (data / "visium_tables" / f"{r}_spots.csv").exists()
           and set(pd.read_csv(data / "visium_tables" / f"{r}_spots.csv", usecols=["barcode"]).barcode) == set(a.obs_names)]
    assert len(run) == 1, f.name
    X = raw.X.tocsr() if hasattr(raw.X, "tocsr") else raw.X
    tot = np.asarray(X.sum(1)).ravel()
    idx = [f"{run[0]}_{b}" for b in raw.obs_names]
    for g in GENES:
        if g in raw.var_names:
            c = np.asarray(X[:, raw.var_names.get_loc(g)].todense()).ravel()
            expr.loc[idx, g] = np.log1p(c / tot * 1e4)

u = um[["UMAP1", "UMAP2"]].to_numpy()
NB = 400
lo, hi = u.min(0), u.max(0)
ij = np.clip(((u - lo) / (hi - lo) * (NB - 1)).round().astype(int), 0, NB - 1)


def density(w):
    h = np.zeros((NB, NB))
    np.add.at(h, (ij[:, 0], ij[:, 1]), w)
    h = gaussian_filter(h, 4)
    return h[ij[:, 0], ij[:, 1]]


order = np.random.default_rng(0).permutation(len(u))
fig, ax = plt.subplots(2, 8, figsize=(24, 6.4))
for j, (name, col) in enumerate(CT):
    v = spots[col].to_numpy()
    sc_ = ax[0, j].scatter(*u[order].T, c=v[order], s=0.2, cmap="viridis",
                           vmax=np.quantile(v, 0.99), rasterized=True)
    ax[0, j].set_title(name, fontsize=10)
for j, g in enumerate(GENES):
    d = density(expr[g].to_numpy())
    o = np.argsort(d)
    ax[1, j].scatter(*u[o].T, c=d[o], s=0.2, cmap="plasma", rasterized=True)
    ax[1, j].set_title(g, fontsize=10)
for a_ in ax.ravel():
    a_.set(xticks=[], yticks=[])
ax[0, 0].set(ylabel="UMAP 2")
ax[1, 0].set(ylabel="UMAP 2", xlabel="UMAP 1")
fig.colorbar(sc_, ax=ax[0, :], shrink=0.6, pad=0.01, label="proportion")
fig.text(0.92, 0.25, "density\nmin → max", fontsize=8)
fig.suptitle("d  composition UMAP: cell-type proportions (top) and marker-gene density (bottom)")
fig.savefig(out / "ed_fig3d.png", dpi=130, bbox_inches="tight")
fig.savefig(out / "ed_fig3d.pdf", bbox_inches="tight")
print("done")
