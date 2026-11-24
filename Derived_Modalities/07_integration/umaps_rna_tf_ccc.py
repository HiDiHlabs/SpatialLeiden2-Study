"""Draft in the style of Kuppe et al. 2022 Extended Data Fig. 1a-b: UMAPs of the
Harmony-integrated Visium spots (all 27 slides, 88,704 spots), coloured by patient (a) and by
zone (b), for the three layers used in the SpatialLeiden panels, plus the joint SpatialLeiden
domains on the same UMAPs: expression only (c, spatial weight 0, K = 9) and with the spatial
graph at weight 0.5 (d). Domains from experiment 4 (exp4_unfixed_k.py: expression resolution
set to give K = 9 at weight 0, spatial resolution 1.0, K free at weight 0.5), seed 0;
results/kuppe_2022/exp4_unfixed_k_{rna,tf,ccc}/exp2_labels.csv.gz. Colours: D1.. by size at
weight 0; weight-0.5 domains take the number of the weight-0 domain they overlap most.

UMAPs (no recomputation):
  RNA  results/kuppe_2022/exp1_harmony_covariates/umaps.npz, key "authors"
       (Harmony on patient + sample + batch; 15-NN; seed 0)
  TF   results/kuppe_2022/tf_harmony/tf_embeddings.npz, key "umap_harmony" (same covariates)
  CCC  colleague's CCC_all27_TFmatched_Harmony.h5ad, obsm "X_umap_harmony_TFmatched"
       (Harmony on patient + sample), mapped to Visium runs by barcode sets (exp4_loader.py)
Patient and zone from data/kuppe_2022/visium_tables/<run>_spots.csv; colours from
common/plot_patient_region_umap.py (team palette).

Usage (from Derived_Modalities/, spatialleiden env):
    python 07_integration/umaps_rna_tf_ccc.py data/kuppe_2022 /path/CCC_all27_TFmatched_Harmony.h5ad \
        figures/kuppe_ed_fig1ab
"""
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE / "common"))
sys.path.insert(0, str(HERE / "03_single_modality_clustering" / "exp2_joint_harmony"))
import exp4_loader  # noqa: E402
from plot_patient_region_umap import PATIENT_COLOURS, REGION_COLOURS, normalise_region  # noqa: E402

DATA, CCC_H5AD, OUT = Path(sys.argv[1]), sys.argv[2], Path(sys.argv[3])
OUT.mkdir(parents=True, exist_ok=True)
RES = HERE / "results" / "kuppe_2022"
mpl.rcParams.update({"font.family": "Arial", "font.size": 6, "pdf.fonttype": 42, "ps.fonttype": 42})

r = np.load(RES / "exp1_harmony_covariates" / "umaps.npz", allow_pickle=True)
names = pd.Index(np.load(RES / "exp1_harmony_covariates" / "embeddings.npz", allow_pickle=True)["obs_names"].astype(str))
t = np.load(RES / "tf_harmony" / "tf_embeddings.npz", allow_pickle=True)
U = {
    "RNA": pd.DataFrame(r["authors"], index=names),
    "TF activity (CollecTRI)": pd.DataFrame(t["umap_harmony"], index=t["obs_names"].astype(str)).loc[names],
    "CCC (ligand-receptor)": exp4_loader.load(f"h5ad:{CCC_H5AD}:X_umap_harmony_TFmatched", DATA).loc[names],
}
meta = []
for s in pd.unique(names.str.split("_", n=1).str[0]):
    d = pd.read_csv(DATA / "visium_tables" / f"{s}_spots.csv", usecols=["barcode", "patient", "major_labl"])
    d.index = s + "_" + d.barcode
    meta.append(d)
meta = pd.concat(meta).loc[names]
zone = meta.major_labl.map(lambda z: normalise_region("control" if z == "CTRL" else z)).to_numpy()
patient = meta.patient.astype(str).to_numpy()
pats = sorted(set(patient), key=lambda p: int(p[1:]))
zones = [z for z in ["control", "RZ", "BZ", "IZ", "FZ"] if z in set(zone)]
pcol = {p: PATIENT_COLOURS.get(p, "#999999") for p in pats}
order = np.random.default_rng(0).permutation(len(names))  # no slide drawn on top of the others

from scipy.optimize import linear_sum_assignment  # noqa: E402
DOM_PAL = ["#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000",
           "#8C6BB1", "#A6761D", "#66C2A5", "#E7298A", "#B3DE69", "#FB8072"]
LABDIR = {"RNA": "exp4_unfixed_k_rna", "TF activity (CollecTRI)": "exp4_unfixed_k_tf",
          "CCC (ligand-receptor)": "exp4_unfixed_k_ccc"}
W_SP = "0.5"


def domains(layer):
    lab = pd.read_csv(RES / LABDIR[layer] / "exp2_labels.csv.gz", index_col=0, dtype=str).loc[names]
    ref = lab["r0_s0"]
    num = {c: i for i, c in enumerate(ref.value_counts().index)}
    d0 = ref.map(num)
    tab = pd.crosstab(lab[f"r{W_SP}_s0"], d0)
    a_, b_ = linear_sum_assignment(-tab.to_numpy())
    m = {tab.index[i]: tab.columns[j] for i, j in zip(a_, b_)}
    nxt = len(num)
    for extra in lab[f"r{W_SP}_s0"][~lab[f"r{W_SP}_s0"].isin(m)].value_counts().index:
        m[extra] = nxt
        nxt += 1
    return d0.to_numpy(), lab[f"r{W_SP}_s0"].map(m).to_numpy()


DOM = {layer: domains(layer) for layer in U}

fig, axes = plt.subplots(4, 3, figsize=(183 / 25.4, 250 / 25.4),
                         gridspec_kw={"wspace": 0.05, "hspace": 0.12, "left": 0.02, "right": 0.86,
                                      "top": 0.96, "bottom": 0.015})
for c, (title, xy) in enumerate(U.items()):
    xy = xy.to_numpy()[order]
    for rrow, (lab, cmap) in enumerate([(patient, pcol), (zone, REGION_COLOURS)]):
        ax = axes[rrow, c]
        ax.scatter(xy[:, 0], xy[:, 1], s=0.15, c=[cmap[v] for v in lab[order]], linewidths=0, rasterized=True)
        ax.set_xticks([])
        ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)
        ax.set_xlabel("UMAP1", fontsize=5, loc="left")
        ax.set_ylabel("UMAP2", fontsize=5, loc="bottom")
        if rrow == 0:
            ax.set_title(f"{title}\n(Harmony, 27 slides, n = {len(names):,} spots)", fontsize=6)
    for rrow, d in [(2, DOM[title][0]), (3, DOM[title][1])]:
        ax = axes[rrow, c]
        dd = d[order]
        ax.scatter(xy[:, 0], xy[:, 1], s=0.15, c=[DOM_PAL[i % len(DOM_PAL)] for i in dd], linewidths=0,
                   rasterized=True)
        ax.set_xticks([])
        ax.set_yticks([])
        for s_ in ax.spines.values():
            s_.set_visible(False)
        k = len(set(d))
        ax.text(0.02, 0.98, f"K = {k}", transform=ax.transAxes, va="top", fontsize=6)
        ax.set_xlabel("UMAP1", fontsize=5, loc="left")
        ax.set_ylabel("UMAP2", fontsize=5, loc="bottom")
h_p = [plt.Line2D([], [], ls="", marker="o", ms=3, mfc=pcol[p], mec="none", label=p) for p in pats]
h_z = [plt.Line2D([], [], ls="", marker="o", ms=3.5, mfc=REGION_COLOURS[z], mec="none", label=z) for z in zones]
fig.legend(handles=h_p, title="Patient", loc="center left", bbox_to_anchor=(0.865, 0.86), ncol=2,
           frameon=False, fontsize=5, title_fontsize=6, handletextpad=0.1, columnspacing=0.6)
fig.legend(handles=h_z, title="Zone", loc="center left", bbox_to_anchor=(0.865, 0.62), frameon=False,
           fontsize=5.5, title_fontsize=6, handletextpad=0.1)
n_dom = max(max(len(set(DOM[k][0])), len(set(DOM[k][1]))) for k in DOM)
h_d = [plt.Line2D([], [], ls="", marker="o", ms=3.5, mfc=DOM_PAL[i % len(DOM_PAL)], mec="none", label=f"D{i + 1}")
       for i in range(n_dom)]
fig.legend(handles=h_d, title="SpatialLeiden\ndomain", loc="center left", bbox_to_anchor=(0.865, 0.25),
           frameon=False, fontsize=5.5, title_fontsize=6, handletextpad=0.1, ncol=2, columnspacing=0.6)
for y, letter, txt in [(0.985, "a", "Patient"), (0.745, "b", "Zone"),
                       (0.505, "c", "SpatialLeiden domains, expression only (spatial weight 0; K = 9)"),
                       (0.262, "d", f"SpatialLeiden domains, spatial weight {W_SP} (K free)")]:
    fig.text(0.0, y, letter, fontsize=8, fontweight="bold", va="top")
    fig.text(0.02, y - 0.002, txt, fontsize=6, va="top")
fig.savefig(OUT / "kuppe_ed_fig1ab_umaps_rna_tf_ccc.pdf", dpi=450)
fig.savefig(OUT / "kuppe_ed_fig1ab_umaps_rna_tf_ccc.png", dpi=300)
print("written", OUT / "kuppe_ed_fig1ab_umaps_rna_tf_ccc.png")
