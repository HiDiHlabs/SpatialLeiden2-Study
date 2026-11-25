"""Yesterday's multipanels (P1 control, P9 IZ, P9 IZ rep2), redrawn with domains from the JOINT
clustering of all 27 slides on Harmony-integrated layers instead of per-slide clustering.

Panels (second argument):
    rna_tf  H&E | RNA only | RNA + space | TF only | TF + space | RNA + TF + space | authors' niches |
            SPP1 | NPPB   (as fig_rna_tf_combo.py on 23 Nov)
    ccc     H&E | RNA only | RNA + space | CCC only | CCC + space | RNA + CCC + space | authors' niches |
            SPP1 -> CD44 (our LIANA+ CCC)   (as fig_modality_combo.py ccc_sp025)
Runs (K = 9 over all 27 slides, seed 0; label columns r<spatial weight>_s0):
    RNA      results/kuppe_2022/exp2_joint_harmony        (Harmony RNA, exp1 "authors" covariates)
    TF       results/kuppe_2022/exp2_joint_harmony_tf     (Harmony CollecTRI ULM, tf_harmony.py)
    CCC      results/kuppe_2022/exp2_joint_harmony_ccc_ours (Harmony on our LIANA+ CCC, ccc_harmony.py)
    RNA + TF / RNA + CCC   results/kuppe_2022/exp3_joint_rna_tf, exp3_joint_rna_ccc_ours (second layer 0.5)
Spatial weight: SPACE_W (default 0.75, as yesterday) for RNA and TF; CCC_W (default 0.25) for the
CCC columns, because the CCC scores are spatially smoothed. Override with FIG_SPACE_W / FIG_CCC_W.
Colours: the domains are shared by all 27 slides, so a colour is one domain on every slide. Each
run's domains are matched to the RNA + space domains by overlap over all 88,704 spots (Hungarian),
so equal colours across columns mean best-matching domains. A slide may show only some of the 9.
Also writes, per slide and column, the number of domains present, the edge-weighted coherence and
the ARI to yesterday's per-slide clustering of the same layers (same weights where it exists).

Usage (from Derived_Modalities/):
    python 03_single_modality_clustering/exp2_joint_harmony/fig_joint_three_slides.py data/kuppe_2022 results/kuppe_2022 rna_tf|ccc
Output: results/kuppe_2022/joint_three_slides/fig_joint_<panel>.png/.pdf + _settings.csv
"""

import json
import os
import sys
from pathlib import Path

import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import EllipseCollection
from PIL import Image
from scipy.optimize import linear_sum_assignment
from scipy.spatial import cKDTree
from sklearn.metrics import adjusted_rand_score

DATA, RES, PANEL = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3]
OUT = RES / "joint_three_slides"
OUT.mkdir(parents=True, exist_ok=True)
W = os.environ.get("FIG_SPACE_W", "0.75")
WC = os.environ.get("FIG_CCC_W", "0.25")
MM = 1 / 25.4
mpl.rcParams.update({"font.family": "Arial", "font.size": 6, "axes.titlesize": 6.5, "pdf.fonttype": 42})
SLIDES = {"10X001": "P1 control\n(10X001)", "ACH0012": "P9 IZ\n(ACH0012)", "ACH0013": "P9 IZ rep2\n(ACH0013)"}
SLIDE_COL = {"10X001": "#0072B2", "ACH0013": "#D55E00", "ACH0012": "#E69F00"}
H5AD = {"10X001": "Visium_control_P1.h5ad", "ACH0013": "Visium_GT_IZ_P9_rep2.h5ad", "ACH0012": "Visium_GT_IZ_P9.h5ad"}
NICHE_COL = {f"ctniche_{i}": c for i, c in enumerate(
    ["#D51F26", "#272E6A", "#208A42", "#89288F", "#F47D2B", "#FEE500", "#8A9FD1", "#C06CAB", "#D8A767"], start=1)}
PAL = [c for i, c in enumerate(mpl.colormaps["tab20"].colors) if i not in (14, 15)] + list(mpl.colormaps["Dark2"].colors)
REF = ("RNA\n+ space", "exp2_joint_harmony", f"r{W}_s0")
# (title, folder, column, yesterday's per-slide (folder, column) for the ARI comparison)
PANELS = {
    "rna_tf": [("RNA only", "exp2_joint_harmony", "r0_s0", ("exp0_k_celltype_niche", "r0_s0")),
               REF + (("exp0_k_celltype_niche", "r0.75_s0"),),
               ("TF only", "exp2_joint_harmony_tf", "r0_s0", ("exp0tf_k_celltype_niche", "r0_s0")),
               ("TF\n+ space", "exp2_joint_harmony_tf", f"r{W}_s0", ("exp0tf_k_celltype_niche", "r0.75_s0")),
               ("RNA + TF\n+ space", "exp3_joint_rna_tf", f"r{W}_s0", ("mm_configs", "rna_tf_s0"))],
    "ccc": [("RNA only", "exp2_joint_harmony", "r0_s0", ("exp0_k_celltype_niche", "r0_s0")),
            REF + (("exp0_k_celltype_niche", "r0.75_s0"),),
            ("CCC only", "exp2_joint_harmony_ccc_ours", "r0_s0", ("exp0ccc_k_celltype_niche", "r0_s0")),
            (f"CCC\n+ space {WC}", "exp2_joint_harmony_ccc_ours", f"r{WC}_s0", ("exp0ccc_k_celltype_niche", "r0.25_s0")),
            (f"RNA + CCC\n+ space {WC}", "exp3_joint_rna_ccc_ours", f"r{WC}_s0", ("mm_configs_sp0.25", "rna_ccc_s0"))],
}
RUNS = PANELS[PANEL]


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


def match_to(ref, lab):
    tab = pd.crosstab(pd.Series(lab), pd.Series(ref))
    r, c = linear_sum_assignment(-tab.to_numpy())
    mapping = {tab.index[i]: tab.columns[j] for i, j in zip(r, c)}
    for n, k in enumerate(k for k in tab.index if k not in mapping):
        mapping[k] = f"x{n}"
    return pd.Series(lab).map(mapping).to_numpy()


# joint labels over all spots; match every run to the reference globally
cache = {}
for name, folder, col, _y in RUNS:
    if folder not in cache:
        cache[folder] = pd.read_csv(RES / folder / "exp2_labels.csv.gz", index_col=0, dtype=str)
ref_all = cache[REF[1]][REF[2]]
labs = {}
for name, folder, col, _y in RUNS:
    lab = cache[folder][col].loc[ref_all.index]
    labs[name] = pd.Series(lab.to_numpy() if name == REF[0] else match_to(ref_all.to_numpy(), lab.to_numpy()),
                           index=ref_all.index)
order = list(ref_all.value_counts().index)  # colour by size of the reference domains
cmap = {d: PAL[k % len(PAL)] for k, d in enumerate(order)}

# features
feat = []
if PANEL == "rna_tf":
    for mk in ["SPP1", "NPPB"]:
        feat.append((f"{mk}\nexpression", mk, "log1p(CP10k)", "magma"))
else:
    feat.append(("SPP1 → CD44\nCCC (ours)", "SPP1^CD44", "local score", "viridis"))
fv = {f[0]: {} for f in feat}
for run in SLIDES:
    h = ad.read_h5ad(DATA / "processed1" / H5AD[run])
    for title, key, unit, cm in feat:
        if PANEL == "rna_tf":
            tot = np.asarray(h.raw.X.sum(1)).ravel()
            g = list(h.raw.var_names).index(key)
            fv[title][run] = pd.Series(np.log1p(np.asarray(h.raw.X[:, g].todense()).ravel() / tot * 1e4), index=h.obs_names)
        else:
            c = pd.read_csv(DATA / "modalities" / "ccc_all27_ours" / f"ccc_{run}.csv", index_col=0)
            fv[title][run] = c[key] if key in c else None
vmax = {t: float(np.quantile(np.concatenate([v.to_numpy() for v in d.values() if v is not None]), 0.99))
        for t, d in fv.items()}

cols = ["H&E"] + [r[0] for r in RUNS] + ["Authors'\nniches"] + [f[0] for f in feat]
J_AUTH = 1 + len(RUNS)
fig, axes = plt.subplots(len(SLIDES), len(cols), figsize=(183 * MM, 32 * MM * len(SLIDES)),
                         gridspec_kw={"wspace": 0.04, "hspace": 0.08})
rows = []
for i, (run, label) in enumerate(SLIDES.items()):
    sp = pd.read_csv(DATA / "visium_tables" / f"{run}_spots.csv").set_index("barcode")
    idx = [x for x in ref_all.index if x.startswith(f"{run}_")]
    bc = [x[len(run) + 1:] for x in idx]
    s = sp.loc[bc]
    xy = s[["pxl_col_fullres", "pxl_row_fullres"]].to_numpy(float)
    pos = {(r, c): n for n, (r, c) in enumerate(zip(s.array_row, s.array_col))}
    ei, ej = zip(*[(n, pos[(r + dr, c + dc)]) for (r, c), n in pos.items()
                   for dr, dc in [(0, 2), (0, -2), (1, 1), (1, -1), (-1, 1), (-1, -1)] if (r + dr, c + dc) in pos])
    ei, ej = np.array(ei), np.array(ej)
    ax = axes[i, 0]
    img = np.asarray(Image.open(DATA / "visium_tables" / f"{run}_lowres.png").convert("RGB"))
    f = json.loads((DATA / "visium_tables" / f"{run}_scalefactors.json").read_text())["tissue_lowres_scalef"]
    ax.imshow(img)
    ax.set_xlim(xy[:, 0].min() * f - 8, xy[:, 0].max() * f + 8)
    ax.set_ylim(xy[:, 1].max() * f + 8, xy[:, 1].min() * f - 8)
    ax.axis("off")
    ax.text(-0.06, 0.5, label, transform=ax.transAxes, rotation=90, ha="right", va="center",
            fontweight="bold", color=SLIDE_COL[run], fontsize=7)
    for j, (name, folder, col, (yf, yc)) in enumerate(RUNS, start=1):
        lab = labs[name].loc[idx].to_numpy()
        draw(axes[i, j], xy, [cmap.get(v, "#777777") for v in lab])
        raw = cache[folder][col].loc[idx].to_numpy()
        yday = pd.read_csv(RES / yf / f"labels_{run}.csv.gz", index_col=0, dtype=str)[yc].loc[bc].to_numpy()
        rows.append({"slide": run, "run": name.replace("\n", " "), "folder": folder, "column": col,
                     "domains_on_slide": len(set(raw)), "coherence": round(float((raw[ei] == raw[ej]).mean()), 3),
                     "ARI_vs_per_slide_23nov": round(adjusted_rand_score(yday, raw), 3),
                     "ARI_vs_celltype_niche": round(adjusted_rand_score(s.celltype_niche, raw), 3)})
    draw(axes[i, J_AUTH], xy, [NICHE_COL[v] for v in s.celltype_niche])
    for jf, (title, key, unit, cm) in enumerate(feat, start=J_AUTH + 1):
        v = fv[title][run]
        if v is None:
            draw(axes[i, jf], xy, [(0.88, 0.88, 0.88, 1)] * len(xy))
            axes[i, jf].text(0.5, 0.5, "pair not\nscored", transform=axes[i, jf].transAxes, ha="center", va="center", fontsize=5.5)
            continue
        norm = mpl.colors.Normalize(0, vmax[title], clip=True)
        draw(axes[i, jf], xy, mpl.colormaps[cm](norm(v.loc[bc].to_numpy())))
for j, c in enumerate(cols):
    axes[0, j].set_title(c, fontweight="bold")
for jf, (title, key, unit, cm) in enumerate(feat, start=J_AUTH + 1):
    cax = axes[-1, jf].inset_axes([0.1, -0.12, 0.8, 0.05])
    cb = fig.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(0, vmax[title]), cmap=cm), cax=cax,
                      orientation="horizontal")
    cb.ax.tick_params(labelsize=5, length=1.5, width=0.4)
    cb.set_label(unit, fontsize=5.5)
    cb.outline.set_linewidth(0.4)
tag = f"fig_joint_{PANEL}" + ("" if (W, WC) == ("0.75", "0.25") else f"_w{W}_c{WC}")
fig.savefig(OUT / f"{tag}.png", dpi=250, bbox_inches="tight")
fig.savefig(OUT / f"{tag}.pdf", dpi=300, bbox_inches="tight")
pd.DataFrame(rows).to_csv(OUT / f"{tag}_settings.csv", index=False)
print(pd.DataFrame(rows).to_string(index=False))
