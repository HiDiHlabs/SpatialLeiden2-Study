"""TF activity for TFs whose behaviour the biology of P1 and P9 predicts.

P9 (GT/IZ, 2 days after symptom onset): infarct core with necrotic cardiomyocytes,
neutrophils and oedema -> hypoxia (HIF1A) and NF-kB inflammation (RELA) should be high.
P1 (non-transplanted donor heart): healthy myocardium -> cardiomyocyte identity (MEF2C)
should be high. MEF2C is the reverse-direction check: ULM scores grow with the number of
detected targets, and P9 is deeper, so depth alone would raise every TF on P9.

Raw ULM scores on one shared colour scale per TF across slides (1st-99th percentile of all
spots), plus the distribution per slide. From collectri_ulm_<run>.h5ad.

Usage (from Derived_Modalities/): python 02_modalities/collectri/collectri_expected_tfs.py data/kuppe_2022 results/kuppe_2022/collectri
"""

import sys
from pathlib import Path

import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import EllipseCollection
from scipy.spatial import cKDTree

DATA, RES = Path(sys.argv[1]), Path(sys.argv[2])
SLIDES = {"10X001": "P1 control", "ACH0013": "P9 IZ rep2", "ACH0012": "P9 IZ"}
SLIDE_COL = {"10X001": "#0072B2", "ACH0013": "#D55E00", "ACH0012": "#E69F00"}
TFS = [("HIF1A", "hypoxia\nexpected P9 > P1"), ("RELA", "NF-κB inflammation\nexpected P9 > P1"),
       ("MEF2C", "cardiomyocyte identity\nexpected P1 > P9")]
mpl.rcParams.update({"font.family": "Arial", "font.size": 6, "axes.titlesize": 6.5, "axes.labelsize": 6,
                     "xtick.labelsize": 5.5, "ytick.labelsize": 5.5, "axes.linewidth": 0.5,
                     "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42})
objs = {r: ad.read_h5ad(DATA / "modalities" / f"collectri_ulm_{r}.h5ad") for r in SLIDES}

fig, axes = plt.subplots(len(TFS), 6, figsize=(183 / 25.4, 1.9 * len(TFS) + 0.3),
                         gridspec_kw={"width_ratios": [1, 1, 1, 0.05, 0.25, 1.1], "wspace": 0.08, "hspace": 0.45})
for i in range(len(TFS)):
    axes[i, 4].axis("off")  # spacer between colour bar and violins
for i, (tf, why) in enumerate(TFS):
    allv = np.concatenate([objs[r][:, tf].X.ravel() for r in SLIDES])
    norm = mpl.colors.Normalize(*np.quantile(allv, [0.01, 0.99]), clip=True)
    for j, (run, label) in enumerate(SLIDES.items()):
        ax = axes[i, j]
        o = objs[run]
        x, y = o.obsm["spatial"].T
        pitch = np.median(cKDTree(np.c_[x, y]).query(np.c_[x, y], k=2)[0][:, 1])
        d = np.full(len(x), 1.05 * pitch)
        coll = EllipseCollection(d, d, np.zeros(len(x)), units="xy", offsets=np.c_[x, y],
                                 offset_transform=ax.transData,
                                 facecolors=mpl.colormaps["magma"](norm(o[:, tf].X.ravel())), linewidths=0)
        coll.set_rasterized(True)
        ax.add_collection(coll)
        ax.set_xlim(x.min() - pitch, x.max() + pitch)
        ax.set_ylim(y.max() + pitch, y.min() - pitch)
        ax.set_aspect("equal")
        ax.axis("off")
        ax.set_title(f"{label} ({run})", color=SLIDE_COL[run], fontweight="bold", fontsize=6)
        if j == 0:  # row heading above the row: TF and the expectation
            ax.text(0.0, 1.22, rf"$\bf{{{tf}}}$   " + why.replace("\n", ", "), transform=ax.transAxes,
                    ha="left", va="bottom", fontsize=7)
    ax = axes[i, 5]
    data = [objs[r][:, tf].X.ravel() for r in SLIDES]
    vp = ax.violinplot(data, showmedians=True, showextrema=False)
    for b, r in zip(vp["bodies"], SLIDES):
        b.set_facecolor(SLIDE_COL[r])
        b.set_alpha(0.75)
    vp["cmedians"].set_color("black")
    ax.set_xticks([1, 2, 3], [SLIDES[r] for r in SLIDES])
    ax.set_title("Distribution per slide (ULM score)", fontsize=6)
    cb = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap="magma"), cax=axes[i, 3])
    cb.ax.tick_params(labelsize=5)
    cb.ax.yaxis.set_ticks_position("right")
fig.suptitle("TF activity (decoupler ULM, CollecTRI) for TFs predicted by the biology; one colour scale per TF, "
             "shared across slides. Raw scores compared across slides also carry the depth difference "
             "(P9 rep2 is ~2.6× deeper than P1); MEF2C, which goes the other way, shows the contrast is "
             "not depth alone.", fontsize=5.5, y=1.0, wrap=True)
out = RES / "figures"
out.mkdir(parents=True, exist_ok=True)
fig.savefig(out / "collectri_expected_tfs_p1_p9.png", dpi=250, bbox_inches="tight")
fig.savefig(out / "collectri_expected_tfs_p1_p9.pdf", dpi=300, bbox_inches="tight")
print("written", out / "collectri_expected_tfs_p1_p9.png")
