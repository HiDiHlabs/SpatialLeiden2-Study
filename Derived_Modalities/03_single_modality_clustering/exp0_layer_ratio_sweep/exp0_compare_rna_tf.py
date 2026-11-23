"""Compare SpatialLeiden on RNA vs on TF activity (CollecTRI) for the same slides.

Reads two finished sweeps (same K setting, same seeds and layer_ratio grid):
    results/kuppe_2022/<rna_dir>/  and  results/kuppe_2022/<tf_dir>/
Only runs that reached K are averaged. Outputs, in <tf_dir>/compare_rna_tf/:
    fig_rna_vs_tf_curves.png/.pdf   coherence gap, ARI to ratio 0, seed stability; RNA solid, TF dashed
    fig_rna_vs_tf_maps_<run>.png    per slide: domains at ratios 0 / 0.5 / 1 / 2 from RNA (top row)
                                    and from TF activity (bottom row), seed 0, plus the authors' niches
    rna_vs_tf_ari.csv              ARI between the RNA and the TF partition, same slide, ratio and seed
    depth_association.csv          how much of each partition is explained by spot depth (UMIs):
                                    eta^2 of log10 UMIs across clusters, at every ratio, both layers

Usage (from Derived_Modalities/):
    python 03_single_modality_clustering/exp0_layer_ratio_sweep/exp0_compare_rna_tf.py \
        data/kuppe_2022 results/kuppe_2022 exp0_k_celltype_niche exp0tf_k_celltype_niche
"""

import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import EllipseCollection
from scipy.optimize import linear_sum_assignment
from scipy.spatial import cKDTree
from sklearn.metrics import adjusted_rand_score

DATA, RES, RNA_DIR, TF_DIR = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3], sys.argv[4]
OUT = RES / TF_DIR / "compare_rna_tf"
OUT.mkdir(parents=True, exist_ok=True)
MM = 1 / 25.4
mpl.rcParams.update({"font.family": "Arial", "font.size": 6, "axes.titlesize": 6.5,
                     "axes.labelsize": 6, "xtick.labelsize": 5.5, "ytick.labelsize": 5.5,
                     "legend.fontsize": 5.5, "axes.linewidth": 0.5, "axes.spines.top": False,
                     "axes.spines.right": False, "pdf.fonttype": 42})
SLIDES = {"ACH0013": "P9 IZ rep2 (ACH0013)", "ACH0012": "P9 IZ (ACH0012)", "10X001": "P1 control (10X001)"}
SLIDE_COL = {"ACH0013": "#D55E00", "ACH0012": "#E69F00", "10X001": "#0072B2"}
NICHE_COL = {f"ctniche_{i}": c for i, c in enumerate(
    ["#D51F26", "#272E6A", "#208A42", "#89288F", "#F47D2B", "#FEE500", "#8A9FD1", "#C06CAB", "#D8A767"], start=1)}
PAL = [c for i, c in enumerate(mpl.colormaps["tab20"].colors) if i not in (14, 15)] + list(mpl.colormaps["Dark2"].colors)


def load(d):
    df = pd.read_csv(RES / d / "exp0_sweep.csv")
    return df[df.n_clusters == df.target_k]


rna, tf = load(RNA_DIR), load(TF_DIR)
K = rna.groupby("slide").target_k.first().to_dict()

# ---------------- curves: RNA solid, TF dashed
fig, axes = plt.subplots(1, 3, figsize=(183 * MM, 55 * MM))
for ax, col, ylab in [(axes[0], "coherence_gap", "Coherence above shuffled"),
                      (axes[1], "ARI_to_ratio0", "ARI to ratio-0 partition"),
                      (axes[2], "seed_ARI_mean", "Seed stability (lower bound)")]:
    for s in SLIDES:
        for d, ls, lab in [(rna, "-", "RNA"), (tf, "--", "TF activity")]:
            m = d[d.slide == s].groupby("layer_ratio")[col].mean()
            ax.plot(m.index, m.values, ls, marker="o", ms=2.2, lw=1, color=SLIDE_COL[s],
                    label=f"{SLIDES[s]}, {lab}")
    ax.set_xlabel("layer_ratio")
    ax.set_ylabel(ylab)
    if col != "coherence_gap":
        ax.set_ylim(0, 1.02)
axes[0].legend(frameon=False, fontsize=5, loc="upper left")
fig.text(0.5, -0.03, "Solid: expression graph from RNA. Dashed: from TF activity (decoupler ULM, CollecTRI). "
         f"Spatial graph, K ({', '.join(f'{SLIDES[s].split(' (')[0]} {K[s]}' for s in SLIDES)}), seeds "
         "and layer_ratio grid identical. Means over seeds that reached K.", ha="center", fontsize=5.5)
fig.tight_layout()
fig.savefig(OUT / "fig_rna_vs_tf_curves.png", dpi=300, bbox_inches="tight")
fig.savefig(OUT / "fig_rna_vs_tf_curves.pdf", bbox_inches="tight")
plt.close(fig)


# ---------------- maps and ARI between layers
def draw(ax, xy, colours, title):
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
    ax.set_title(title)


def match_to(ref, lab):
    tab = pd.crosstab(pd.Series(lab), pd.Series(ref))
    r, c = linear_sum_assignment(-tab.to_numpy())
    mapping = {tab.index[i]: tab.columns[j] for i, j in zip(r, c)}
    for n, k in enumerate(k for k in tab.index if k not in mapping):  # unmatched clusters
        mapping[k] = str(100 + n)
    return pd.Series(lab).map(mapping).to_numpy()


def eta2(values, labels):
    """Share of the variance of `values` explained by cluster membership (0-1)."""
    v = pd.Series(values)
    g = v.groupby(pd.Series(labels).to_numpy())
    return float((g.size() * (g.mean() - v.mean()) ** 2).sum() / ((v - v.mean()) ** 2).sum())


ari_rows, depth_rows = [], []
SHOW = [0, 0.5, 1.0, 2.0]
for s in SLIDES:
    lr = pd.read_csv(RES / RNA_DIR / f"labels_{s}.csv.gz", index_col=0, dtype=str)
    lt = pd.read_csv(RES / TF_DIR / f"labels_{s}.csv.gz", index_col=0, dtype=str)
    assert list(lr.index) == list(lt.index), "spot order differs between layers"
    umi = pd.read_csv(DATA / "visium_tables" / f"{s}_spots.csv").set_index("barcode").loc[lr.index, "nCount_Spatial"]
    log_umi = np.log10(umi.to_numpy())
    ratios = sorted({float(c.split("_s")[0][1:]) for c in lr.columns if c.startswith("r")})
    seeds = sorted({int(c.split("_s")[1]) for c in lr.columns if c.startswith("r")})
    for ratio in ratios:
        for seed in seeds:
            key = f"r{ratio:g}_s{seed}" if f"r{ratio:g}_s{seed}" in lr else f"r{ratio}_s{seed}"
            if key not in lr or key not in lt:
                key = f"r{int(ratio)}_s{seed}" if ratio == int(ratio) else key
            ari_rows.append({"slide": s, "layer_ratio": ratio, "seed": seed,
                             "ARI_rna_vs_tf": adjusted_rand_score(lr[key], lt[key])})
            for layer, lab in [("RNA", lr[key]), ("TF", lt[key])]:
                depth_rows.append({"slide": s, "layer": layer, "layer_ratio": ratio, "seed": seed,
                                   "eta2_log10_umi": eta2(log_umi, lab.to_numpy())})
    xy = lr[["x", "y"]].astype(float).to_numpy()
    fig, axes = plt.subplots(2, len(SHOW) + 1, figsize=(183 * MM, 78 * MM),
                             gridspec_kw={"wspace": 0.04, "hspace": 0.18})
    for row, (name, lab_df) in enumerate([("RNA", lr), ("TF activity", lt)]):
        ref = lab_df["r0_s0"].to_numpy()
        cats = sorted(set(ref), key=int)
        col = {c: PAL[k % len(PAL)] for k, c in enumerate(cats)}
        for j, ratio in enumerate(SHOW):
            key = "r0_s0" if ratio == 0 else f"r{ratio}_s0"
            m = ref if ratio == 0 else match_to(ref, lab_df[key].to_numpy())
            draw(axes[row, j], xy, [col.get(v, "#777777") for v in m], f"{name}, layer_ratio {ratio:g}")
    draw(axes[0, -1], xy, [NICHE_COL[v] for v in lr.celltype_niche], "Authors' cell-type niches")
    axes[1, -1].axis("off")  # ARI between layers is kept in rna_vs_tf_ari.csv, not on the maps
    fig.suptitle(f"{SLIDES[s]}: SpatialLeiden on RNA (top) vs TF activity (bottom), K = {K[s]}, seed 0; "
                 "colours matched to the ratio-0 partition of each row", fontsize=6.5, y=0.99)
    fig.savefig(OUT / f"fig_rna_vs_tf_maps_{s}.png", dpi=250, bbox_inches="tight")
    fig.savefig(OUT / f"fig_rna_vs_tf_maps_{s}.pdf", dpi=300, bbox_inches="tight")
    plt.close(fig)

ari = pd.DataFrame(ari_rows)
ari.to_csv(OUT / "rna_vs_tf_ari.csv", index=False)
dep = pd.DataFrame(depth_rows)
dep.to_csv(OUT / "depth_association.csv", index=False)
print(ari.groupby(["slide", "layer_ratio"]).ARI_rna_vs_tf.mean().unstack(0).round(2).to_string())
print("\nShare of log10-UMI variance explained by cluster (eta^2), mean of seeds:")
print(dep.groupby(["slide", "layer", "layer_ratio"]).eta2_log10_umi.mean().unstack([0, 1]).round(2).to_string())
