"""Multipanels comparing SpatialLeiden on RNA with the derived layers (CCC, scCellFie
metabolic tasks) and all layers together. Same layout and style as fig_rna_tf_combo.py.

Panels (choose with the second argument):
    ccc   H&E | RNA only | RNA + space | CCC only | CCC + space | RNA + CCC + space |
          authors' niches | SPP1->CD44 score
    ccc_sp025  H&E | RNA + space | CCC only | CCC + space 0.25 / 0.75 | RNA + CCC + space 0.25 / 0.75 |
          authors' niches | SPP1->CD44 score (spatial weight 0.25 = less double-counting of space;
          RNA + CCC at 0.25 from mm_configs.py run with MM_SPATIAL_W=0.25 -> mm_configs_sp0.25)
    met   H&E | RNA only | RNA + space | metabolic only | metabolic + space |
          RNA + metabolic + space | authors' niches | two scCellFie task maps
    all   H&E | RNA + space | TF + space | CCC + space | metabolic + space |
          RNA + TF | RNA + CCC | RNA + metabolic | RNA + all (each + space) | authors' niches
Runs (all K = authors' cell-type niche count, seed 0):
    X only / X + space   experiment-0 sweep on layer X, layer_ratio 0 / 0.75
                         (results/kuppe_2022/exp0{,tf,ccc,met}_k_celltype_niche, r0_s0 / r0.75_s0)
    RNA + X + space      mm_configs.py: RNA 1, each extra layer 0.5, spatial 0.75 (label <config>_s0)
Colours: domains matched (Hungarian, by overlap) to RNA + space. Caveat for CCC: the team's
scores are already spatially smoothed, so "CCC + space" partly counts space twice.

Usage (from Derived_Modalities/):
    python 04_multimodal_weight_grid/proto_rna_tf/fig_modality_combo.py data/kuppe_2022 results/kuppe_2022 ccc|ccc_sp025|met|all
Output: results/kuppe_2022/mm_configs/figures/fig_<panel>_combo.png/.pdf + _settings.csv
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
from PIL import Image
from scipy.optimize import linear_sum_assignment
from scipy.spatial import cKDTree

DATA, RES, PANEL = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3]
OUT = RES / "mm_configs" / "figures"
OUT.mkdir(parents=True, exist_ok=True)
MM = 1 / 25.4
mpl.rcParams.update({"font.family": "Arial", "font.size": 6, "axes.titlesize": 6.5, "pdf.fonttype": 42})
SLIDES = {"10X001": "P1 control\n(10X001)", "ACH0012": "P9 IZ\n(ACH0012)", "ACH0013": "P9 IZ rep2\n(ACH0013)"}
SLIDE_COL = {"10X001": "#0072B2", "ACH0013": "#D55E00", "ACH0012": "#E69F00"}
NICHE_COL = {f"ctniche_{i}": c for i, c in enumerate(
    ["#D51F26", "#272E6A", "#208A42", "#89288F", "#F47D2B", "#FEE500", "#8A9FD1", "#C06CAB", "#D8A767"], start=1)}
PAL = [c for i, c in enumerate(mpl.colormaps["tab20"].colors) if i not in (14, 15)] + list(mpl.colormaps["Dark2"].colors)
RNA0 = ("RNA only", "exp0_k_celltype_niche", "r0_s0")
RNASP = ("RNA\n+ space", "exp0_k_celltype_niche", "r0.75_s0")
PANELS = {
    "ccc": [RNA0, RNASP, ("CCC only", "exp0ccc_k_celltype_niche", "r0_s0"),
            ("CCC\n+ space", "exp0ccc_k_celltype_niche", "r0.75_s0"),
            ("RNA + CCC\n+ space", "mm_configs", "rna_ccc_s0")],
    "ccc_sp025": [RNASP, ("CCC only", "exp0ccc_k_celltype_niche", "r0_s0"),
                  ("CCC + space\n0.25", "exp0ccc_k_celltype_niche", "r0.25_s0"),
                  ("CCC + space\n0.75", "exp0ccc_k_celltype_niche", "r0.75_s0"),
                  ("RNA + CCC\n+ space 0.25", "mm_configs_sp0.25", "rna_ccc_s0"),
                  ("RNA + CCC\n+ space 0.75", "mm_configs", "rna_ccc_s0")],
    "met": [RNA0, RNASP, ("Metabolic\nonly", "exp0met_k_celltype_niche", "r0_s0"),
            ("Metabolic\n+ space", "exp0met_k_celltype_niche", "r0.75_s0"),
            ("RNA + metab.\n+ space", "mm_configs", "rna_met_s0")],
    "all": [RNASP, ("TF\n+ space", "exp0tf_k_celltype_niche", "r0.75_s0"),
            ("CCC\n+ space", "exp0ccc_k_celltype_niche", "r0.75_s0"),
            ("Metabolic\n+ space", "exp0met_k_celltype_niche", "r0.75_s0"),
            ("RNA + TF\n+ space", "mm_configs", "rna_tf_s0"),
            ("RNA + CCC\n+ space", "mm_configs", "rna_ccc_s0"),
            ("RNA + metab.\n+ space", "mm_configs", "rna_met_s0"),
            ("RNA + all\n+ space", "mm_configs", "rna_all_s0")],
}
RUNS = PANELS[PANEL]
H5AD = {"10X001": "Visium_control_P1.h5ad", "ACH0013": "Visium_GT_IZ_P9_rep2.h5ad", "ACH0012": "Visium_GT_IZ_P9.h5ad"}


def features():
    """Feature columns: {title: {run: Series indexed by barcode, or None if absent}}, unit, cmap."""
    if PANEL.startswith("ccc"):
        tabs = {r: pd.read_csv(DATA / "modalities" / f"ccc_{r}.csv", index_col=0) for r in SLIDES}
        return [("SPP1 → CD44\nCCC score", {r: t["SPP1^CD44"] if "SPP1^CD44" in t else None for r, t in tabs.items()},
                 "local score", "viridis")]
    if PANEL == "met":
        tasks = [("Palmitate degradation", "Fatty-acid\noxidation\n(palmitate degr.)"),
                 ("ATP generation from glucose (hypoxic conditions) - glycolysis", "Anaerobic\nglycolysis\n(hypoxic ATP)")]
        out = []
        for task, title in tasks:
            vals = {}
            for r in SLIDES:
                a = ad.read_h5ad(DATA / "modalities" / f"sccellfie_tasks_{r}.h5ad")
                vals[r] = pd.Series(np.asarray(a[:, task].X.todense() if hasattr(a.X, "todense") else a[:, task].X).ravel(),
                                    index=a.obs_names)
            out.append((title, vals, "task activity", "viridis"))
        return out
    return []


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


FEAT = features()
vmax = {}
for title, vals, _, _ in FEAT:
    vmax[title] = float(np.quantile(np.concatenate([v.to_numpy() for v in vals.values() if v is not None]), 0.99))
cols = ["H&E"] + [r[0] for r in RUNS] + ["Authors'\nniches"] + [f[0] for f in FEAT]
J_AUTH = 1 + len(RUNS)
fig, axes = plt.subplots(len(SLIDES), len(cols), figsize=(183 * MM, (183 / len(cols)) * 1.0 * MM * len(SLIDES) + 10 * MM),
                         gridspec_kw={"wspace": 0.04, "hspace": 0.08})
rows = []
for i, (run, label) in enumerate(SLIDES.items()):
    spots = pd.read_csv(DATA / "visium_tables" / f"{run}_spots.csv").set_index("barcode")
    labs = {}
    for name, folder, key in RUNS + [RNASP]:
        labs[name] = pd.read_csv(RES / folder / f"labels_{run}.csv.gz", index_col=0, dtype=str)[key]
    idx = labs[RNASP[0]].index
    for name in labs:
        assert list(labs[name].index) == list(idx), "spot order differs between runs"
    sp = spots.loc[idx]
    xy = sp[["pxl_col_fullres", "pxl_row_fullres"]].to_numpy(float)
    pos = {(r, c): n for n, (r, c) in enumerate(zip(sp.array_row, sp.array_col))}
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
    ref = labs[RNASP[0]].to_numpy()
    cmap = {c: PAL[k % len(PAL)] for k, c in enumerate(sorted(set(ref), key=int))}
    for j, (name, folder, key) in enumerate(RUNS, start=1):
        lab = labs[name].to_numpy()
        m = ref if name == RNASP[0] else match_to(ref, lab)
        draw(axes[i, j], xy, [cmap.get(v, "#777777") for v in m])
        rows.append({"slide": run, "run": name.replace("\n", " "), "source": folder, "key": key,
                     "n_domains": len(set(lab)), "coherence_edge_weighted": round(float((lab[ei] == lab[ej]).mean()), 3)})
    draw(axes[i, J_AUTH], xy, [NICHE_COL[v] for v in sp.celltype_niche])
    for jf, (title, vals, unit, cm) in enumerate(FEAT, start=J_AUTH + 1):
        ax = axes[i, jf]
        v = vals[run]
        if v is None:
            draw(ax, xy, [(0.88, 0.88, 0.88, 1)] * len(xy))
            ax.text(0.5, 0.5, "not in this\nslide's results", transform=ax.transAxes, ha="center", va="center",
                    fontsize=5.5, color="#444444")
            continue
        v = v.loc[idx].to_numpy()
        draw(ax, xy, mpl.colormaps[cm](mpl.colors.Normalize(0, vmax[title], clip=True)(v)))
        rows.append({"slide": run, "run": title.replace("\n", " "), "source": unit, "key": "",
                     "median": round(float(np.median(v)), 3), "share_spots_gt0": round(float((v > 0).mean()), 3)})
for j, c in enumerate(cols):
    axes[0, j].set_title(c, fontweight="bold", fontsize=6)
for jf, (title, _, unit, cm) in enumerate(FEAT, start=J_AUTH + 1):
    cax = axes[-1, jf].inset_axes([0.1, -0.12, 0.8, 0.05])
    cb = fig.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(0, vmax[title]), cmap=cm), cax=cax,
                      orientation="horizontal")
    cb.ax.tick_params(labelsize=5, length=1.5, width=0.4)
    cb.set_label(unit, fontsize=5.5)
    cb.outline.set_linewidth(0.4)
fig.savefig(OUT / f"fig_{PANEL}_combo.png", dpi=250, bbox_inches="tight")
fig.savefig(OUT / f"fig_{PANEL}_combo.pdf", dpi=300, bbox_inches="tight")
pd.DataFrame(rows).to_csv(OUT / f"fig_{PANEL}_combo_settings.csv", index=False)
print(pd.DataFrame(rows).to_string(index=False))
