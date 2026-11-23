"""Slide atlas: everything that ships with each Kuppe et al. 2022 Visium run, one page per run.

Per page: H&E, the authors' cell-type and molecular niches, the dominant cell type per
spot (cardiomyocyte / fibroblast / myeloid / other, as in the paper's Fig. 2e), depth and
QC maps (UMIs, genes, % mitochondrial, estimated cells per spot) and the 11 shipped
cell2location proportions. Descriptive only: nothing is recomputed from expression.

Inputs (from Derived_Modalities/): data/kuppe_2022/visium_tables (made by
01_data_preprocess/kuppe_visium_exploration/run_all.sh) and the cohort table
(results/kuppe_2022/cohort/cohort_table.csv, made by 01_data_preprocess/kuppe_cohort).

Usage (from Derived_Modalities/):
    python 01_data_preprocess/kuppe_slide_atlas/slide_atlas.py data/kuppe_2022 results/kuppe_2022/slide_atlas
"""

import json
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.collections import EllipseCollection
from matplotlib.lines import Line2D
from PIL import Image
from scipy.spatial import cKDTree

DATA, OUT = Path(sys.argv[1]), Path(sys.argv[2])
OUT.mkdir(parents=True, exist_ok=True)
TABLES = DATA / "visium_tables"
MM = 1 / 25.4

mpl.rcParams.update(
    {
        "font.family": "Arial",
        "font.size": 6,
        "axes.titlesize": 6.5,
        "legend.fontsize": 5.5,
        "axes.linewidth": 0.4,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)

ZONES = ["CTRL", "RZ", "BZ", "IZ", "FZ"]
ZONE_LABEL = {
    "CTRL": "control",
    "RZ": "remote",
    "BZ": "border",
    "IZ": "ischaemic",
    "FZ": "fibrotic",
}
ZONE_COL = {
    "CTRL": "#0072B2",
    "RZ": "#009E73",
    "BZ": "#E69F00",
    "IZ": "#D55E00",
    "FZ": "#CC79A7",
}
GROUP_LABEL = {"group_1": "myogenic", "group_2": "ischaemic", "group_3": "fibrotic"}
# authors' niche colours (saezlab/visium_heart, plot_niche_example.R)
CT_COL = {
    1: "#D51F26",
    2: "#272E6A",
    3: "#208A42",
    4: "#89288F",
    5: "#F47D2B",
    6: "#FEE500",
    7: "#8A9FD1",
    8: "#C06CAB",
    9: "#D8A767",
}
# dominant cell type, as in the paper's Fig. 2e
DOM_COL = {
    "Cardiomyocyte": "#C0392B",
    "Fibroblast": "#9FD3EA",
    "Myeloid": "#E6A23C",
    "Other": "#B3B3B3",
}
CT_NAME = {
    "Adipo": "Adipocytes",
    "CM": "Cardiomyocytes",
    "Endo": "Endothelial",
    "Fib": "Fibroblasts",
    "Lymphoid": "Lymphoid",
    "Mast": "Mast",
    "Myeloid": "Myeloid",
    "Neuronal": "Neuronal",
    "PC": "Pericytes",
    "prolif": "Proliferating",
    "vSMCs": "vSMCs",
}

# ------------------------------------------------------------------ data
spots = pd.concat(
    [pd.read_csv(f) for f in sorted(TABLES.glob("*_spots.csv"))], ignore_index=True
)
prop_cols = [c for c in spots.columns if c.startswith("prop_")]
spots["ct_niche"] = spots.celltype_niche.str.removeprefix("ctniche_").astype(int)
spots["mol_niche"] = spots.molecular_niche.str.removeprefix("molniche_").astype(int)
dom = spots[prop_cols].idxmax(axis=1).str.removeprefix("prop_")
spots["dominant"] = dom.map(
    {"CM": "Cardiomyocyte", "Fib": "Fibroblast", "Myeloid": "Myeloid"}
).fillna("Other")

mol_levels = sorted(spots.mol_niche.unique())
mol_pal = [mpl.colormaps["tab20"](i) for i in range(20)]
MOL_COL = {
    k: mol_pal[(2 * i) % 20 if i < 10 else (2 * (i - 10) + 1) % 20]
    for i, k in enumerate(mol_levels)
}

# one colour scale per quantity across all slides (99th percentile), so pages are comparable
VMAX = {c: float(spots[c].quantile(0.99)) for c in prop_cols}
QC = [
    ("nCount_Spatial", "UMIs per spot (log10)", True),
    ("nFeature_Spatial", "Genes per spot (log10)", True),
    ("percent.mt", "% mitochondrial UMIs", False),
    ("c2l_total_abundance", "Estimated cells per spot", False),
]
QC_RANGE = {}
for col, _, log in QC:
    v = np.log10(spots[col].clip(lower=1)) if log else spots[col]
    QC_RANGE[col] = (float(v.quantile(0.01)), float(v.quantile(0.99)))

cohort_file = (
    DATA.parent.parent / "results" / "kuppe_2022" / "cohort" / "cohort_table.csv"
)
days = {}
if cohort_file.exists():
    ct = pd.read_csv(cohort_file)
    days = dict(zip(ct.patient_region_id, ct.days_after_onset))

runs = (
    spots.groupby("hca_sample_id")
    .agg(
        section=("patient_region_id", "first"),
        patient=("patient", "first"),
        zone=("major_labl", "first"),
        group=("patient_group", "first"),
        spots=("barcode", "size"),
        umi=("nCount_Spatial", "median"),
        genes=("nFeature_Spatial", "median"),
    )
    .reset_index()
)
runs["z"] = runs.zone.map(ZONES.index)
runs = runs.sort_values(["z", "section", "hca_sample_id"]).reset_index(drop=True)


# ------------------------------------------------------------------ drawing helpers
def spot_layer(ax, df, colours, sid):
    x, y = df.pxl_col_fullres.to_numpy(float), df.pxl_row_fullres.to_numpy(float)
    pitch = np.median(cKDTree(np.c_[x, y]).query(np.c_[x, y], k=2)[0][:, 1])
    d = np.full(len(x), 1.05 * pitch)
    coll = EllipseCollection(
        d,
        d,
        np.zeros(len(x)),
        units="xy",
        offsets=np.c_[x, y],
        offset_transform=ax.transData,
        facecolors=colours,
        linewidths=0,
    )
    coll.set_rasterized(True)
    ax.add_collection(coll)
    pad = 1.5 * pitch
    ax.set_xlim(x.min() - pad, x.max() + pad)
    ax.set_ylim(y.max() + pad, y.min() - pad)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


def continuous(ax, df, values, vmin, vmax, cmap, title, sid):
    norm = mpl.colors.Normalize(vmin, vmax, clip=True)
    spot_layer(ax, df, mpl.colormaps[cmap](norm(values)), sid)
    ax.set_title(title)
    sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    cb = plt.colorbar(sm, ax=ax, fraction=0.04, pad=0.02, shrink=0.8)
    cb.ax.tick_params(labelsize=5, width=0.3, length=1.5)
    cb.outline.set_linewidth(0.3)


def legend_below(ax, colours, labels, ncol):
    h = [
        Line2D([], [], ls="", marker="o", ms=3.5, mfc=c, mec="none", label=lab)
        for c, lab in zip(colours, labels)
    ]
    ax.legend(
        handles=h,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.02),
        ncol=ncol,
        frameon=False,
        handletextpad=0.1,
        columnspacing=0.6,
        labelspacing=0.2,
    )


# ------------------------------------------------------------------ pages
summary = []
with PdfPages(OUT / "kuppe_visium_slide_atlas.pdf") as pdf:
    # cover page
    fig = plt.figure(figsize=(210 * MM, 297 * MM))
    fig.text(
        0.08,
        0.93,
        "Kuppe et al. 2022: Visium slide atlas",
        fontsize=12,
        fontweight="bold",
    )
    fig.text(
        0.08,
        0.905,
        "Everything that ships with each Visium run, one page per run. Descriptive only.",
        fontsize=8,
    )
    lines = [
        "Per page: H&E (Space Ranger low resolution) · cell-type niches (authors' Louvain clustering of",
        "cell2location proportions; authors' colours) · molecular niches (authors' clustering of spot expression)",
        "· dominant cell type per spot (highest cell2location proportion; as in the paper's Fig. 2e) · UMIs,",
        "genes, % mitochondrial UMIs and estimated cells per spot · the 11 cell2location proportions.",
        "Neither niche set used spot coordinates. Colour scales are shared across all pages (1st-99th",
        "percentile over all 91,517 spots), so slides can be compared directly.",
    ]
    for i, t in enumerate(lines):
        fig.text(0.08, 0.87 - i * 0.017, t, fontsize=7)
    fig.text(0.08, 0.74, "Runs (ordered by zone)", fontsize=8, fontweight="bold")
    for i, r in runs.iterrows():
        d = days.get(r.section)
        dtxt = (
            "control"
            if r.zone == "CTRL"
            else (f"{int(d)} d" if d == d and d is not None else "n.a.")
        )
        fig.text(
            0.08,
            0.72 - i * 0.021,
            f"p. {i + 2:>2}   {r.section:<16} {r.hca_sample_id:<8} {r.patient:<4} "
            f"{ZONE_LABEL[r.zone]:<10} {GROUP_LABEL.get(r.group, r.group):<10} {dtxt:<8} "
            f"{r.spots:>5,} spots   median {r.umi:>6,.0f} UMIs",
            fontsize=7,
            family="monospace",
            color=ZONE_COL[r.zone],
        )
    pdf.savefig(fig)
    plt.close(fig)

    for i, r in runs.iterrows():
        sid = r.hca_sample_id
        df = spots[spots.hca_sample_id == sid]
        fig, axes = plt.subplots(
            5,
            4,
            figsize=(210 * MM, 297 * MM),
            gridspec_kw={
                "hspace": 0.45,
                "wspace": 0.12,
                "left": 0.03,
                "right": 0.97,
                "top": 0.92,
                "bottom": 0.03,
            },
        )
        d = days.get(r.section)
        dtxt = (
            "control donor"
            if r.zone == "CTRL"
            else (
                f"{int(d)} days after symptom onset" if d == d and d is not None else ""
            )
        )
        fig.text(
            0.03,
            0.965,
            f"{r.section}  ·  run {sid}  ·  patient {r.patient}",
            fontsize=10,
            fontweight="bold",
            color=ZONE_COL[r.zone],
        )
        fig.text(
            0.03,
            0.945,
            f"zone: {ZONE_LABEL[r.zone]}  ·  pathologist group: {GROUP_LABEL.get(r.group, r.group)}"
            f"  ·  {dtxt}  ·  {r.spots:,} spots  ·  median {r.umi:,.0f} UMIs, {r.genes:,.0f} genes per spot",
            fontsize=7,
        )

        # row 1: H&E, niches, dominant type
        ax = axes[0, 0]
        img = np.asarray(Image.open(TABLES / f"{sid}_lowres.png").convert("RGB"))
        f = json.loads((TABLES / f"{sid}_scalefactors.json").read_text())[
            "tissue_lowres_scalef"
        ]
        x, y = df.pxl_col_fullres * f, df.pxl_row_fullres * f
        ax.imshow(img)
        ax.set_xlim(x.min() - 10, x.max() + 10)
        ax.set_ylim(y.max() + 10, y.min() - 10)
        ax.set_xticks([])
        ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)
        ax.set_title("H&E")

        spot_layer(axes[0, 1], df, df.ct_niche.map(CT_COL).to_list(), sid)
        axes[0, 1].set_title("Cell-type niche (authors)")
        present = sorted(df.ct_niche.unique())
        legend_below(
            axes[0, 1], [CT_COL[k] for k in present], [str(k) for k in present], ncol=5
        )

        spot_layer(axes[0, 2], df, df.mol_niche.map(MOL_COL).to_list(), sid)
        axes[0, 2].set_title("Molecular niche (authors)")
        present = sorted(df.mol_niche.unique())
        legend_below(
            axes[0, 2], [MOL_COL[k] for k in present], [str(k) for k in present], ncol=6
        )

        spot_layer(axes[0, 3], df, df.dominant.map(DOM_COL).to_list(), sid)
        axes[0, 3].set_title("Dominant cell type")
        legend_below(
            axes[0, 3],
            list(DOM_COL.values()),
            ["CM", "Fib", "Myeloid", "Other"],
            ncol=4,
        )

        # row 2: depth and QC
        for ax, (col, title, log) in zip(axes[1], QC):
            v = np.log10(df[col].clip(lower=1)) if log else df[col]
            continuous(ax, df, v.to_numpy(), *QC_RANGE[col], "viridis", title, sid)

        # rows 3-5: cell2location proportions
        for ax, col in zip(axes[2:].flat, prop_cols):
            continuous(
                ax,
                df,
                df[col].to_numpy(),
                0,
                VMAX[col],
                "magma",
                f"{CT_NAME.get(col[5:], col[5:])} (proportion)",
                sid,
            )
        axes[4, 3].axis("off")
        axes[4, 3].text(
            0.3,
            0.6,
            "Colour scales are shared across\nall pages: 0 to the 99th percentile\n"
            "of each quantity over all spots.\nQC scales: 1st to 99th percentile.",
            fontsize=6,
            va="top",
            transform=axes[4, 3].transAxes,
        )

        pdf.savefig(fig, dpi=150)  # rasterised spot layers; keeps the PDF shareable
        plt.close(fig)
        summary.append(
            {
                "page": i + 2,
                "run": sid,
                "section": r.section,
                "zone": r.zone,
                "spots": int(r.spots),
            }
        )
        print(f"page {i + 2}: {r.section} ({sid})", flush=True)

pd.DataFrame(summary).to_csv(OUT / "slide_atlas_pages.csv", index=False)
print("written", OUT / "kuppe_visium_slide_atlas.pdf")
