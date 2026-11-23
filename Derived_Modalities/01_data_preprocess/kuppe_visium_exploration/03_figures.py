"""Descriptive figures for the Kuppe et al. 2022 Visium sections.

Uses only what ships with the data: the authors' Seurat objects (extracted by
01_extract_seurat_tables.R), Space Ranger files (02_extract_spaceranger_files.py),
the Zenodo h5ad files and the metadata sheets. No new analysis of expression.

Usage: python 03_figures.py <data/kuppe_2022 dir> <output dir>
"""

import json
import sys
from pathlib import Path

import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch
from PIL import Image

Image.MAX_IMAGE_PIXELS = None
DATA, OUT = Path(sys.argv[1]), Path(sys.argv[2])
FIG = OUT / "figures"
FIG.mkdir(parents=True, exist_ok=True)
TABLES = DATA / "visium_tables"

mpl.rcParams.update(
    {
        "font.size": 7,
        "axes.titlesize": 8,
        "axes.labelsize": 7,
        "xtick.labelsize": 6,
        "ytick.labelsize": 6,
        "legend.fontsize": 6,
        "pdf.fonttype": 42,
        "axes.spines.top": False,
        "axes.spines.right": False,
    }
)

# Zones: fixed order and colourblind-safe colours (Okabe-Ito) used in every figure.
ZONES = ["CTRL", "RZ", "BZ", "IZ", "FZ"]
ZONE_NAME = {
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
GROUP_NAME = {"group_1": "myogenic", "group_2": "ischaemic", "group_3": "fibrotic"}
A4W = 7.0  # usable width in inches

# ---------------------------------------------------------------- metadata
vis_meta = pd.read_csv(DATA / "visium_raw" / "metadata-Visium.csv")
he_meta = pd.read_excel(DATA / "he_images" / "metadata.xlsx")


def sheet(path):
    x = pd.read_excel(path, header=1)
    return x[x["sample_id"].notna()]


rna_meta = sheet(DATA / "metadata" / "metadata-RNA-seq.xlsx")
atac_meta = sheet(DATA / "metadata" / "metadata-ATAC-seq.xlsx")

spots = []
for f in sorted(TABLES.glob("*_spots.csv")):
    spots.append(pd.read_csv(f))
spots = pd.concat(spots, ignore_index=True)
spots["zone"] = pd.Categorical(spots["major_labl"], ZONES, ordered=True)
prop_cols = [c for c in spots.columns if c.startswith("prop_")]

scalef = {
    f.name.removesuffix("_scalefactors.json"): json.loads(f.read_text())
    for f in TABLES.glob("*_scalefactors.json")
}
metrics = pd.concat(
    [
        pd.read_csv(f).assign(hca_sample_id=f.name.removesuffix("_metrics_summary.csv"))
        for f in TABLES.glob("*_metrics_summary.csv")
    ],
    ignore_index=True,
)

# one row per Visium section
sec = (
    spots.groupby("hca_sample_id")
    .agg(
        patient_region_id=("patient_region_id", "first"),
        patient=("patient", "first"),
        patient_group=("patient_group", "first"),
        zone=("major_labl", "first"),
        batch=("batch", "first"),
        n_spots=("barcode", "size"),
        median_umi=("nCount_Spatial", "median"),
        median_genes=("nFeature_Spatial", "median"),
        median_mt=("percent.mt", "median"),
        median_cells=("c2l_total_abundance", "median"),
        n_ctniche=("celltype_niche", "nunique"),
        n_molniche=("molecular_niche", "nunique"),
    )
    .reset_index()
)
sec["zone"] = pd.Categorical(sec["zone"], ZONES, ordered=True)
sec["label"] = sec["patient_region_id"] + " (" + sec["hca_sample_id"] + ")"
sec = sec.sort_values(["zone", "patient_region_id", "hca_sample_id"]).reset_index(
    drop=True
)
ORDER = list(sec["hca_sample_id"])
LABEL = dict(zip(sec["hca_sample_id"], sec["label"]))

# Zenodo processed h5ad per section (27 files; one of the two RZ_P6 runs)
h5ad_info = {}
for f in sorted((DATA / "processed1").glob("Visium_*.h5ad")):
    a = ad.read_h5ad(f, backed="r")
    h5ad_info[f.stem.removeprefix("Visium_")] = (a.n_obs, a.n_vars, list(a.obs_names))
    a.file.close()


def region_key(r):  # "GT/IZ_P9_rep2" -> "GT_IZ_P9_rep2" as in the h5ad names
    return r.replace("/", "_")


stats = {"n_sections_seurat": len(sec), "n_spots_total": int(len(spots))}

# ------------------------------------------------- 1. modality availability
regions = sorted(
    set(vis_meta.patient_region_id)
    | set(rna_meta.patient_region_id)
    | set(atac_meta.patient_region_id)
    | set(he_meta.patient_region_id)
)
zone_of = {}
for df in (vis_meta, he_meta):
    zone_of.update(dict(zip(df.patient_region_id, df.major_labl)))
for df in (rna_meta, atac_meta):
    for r, z in zip(df.patient_region_id, df.region_novel):
        zone_of.setdefault(r, {"control": "CTRL"}.get(z, z.split("/")[0]))
patient_of = {}
for df in (vis_meta, rna_meta, atac_meta, he_meta):
    patient_of.update(dict(zip(df.patient_region_id, df.patient)))

avail_cols = [
    "snRNA-seq\n(nuclei)",
    "snATAC-seq\n(nuclei)",
    "Visium h5ad\n(Zenodo processed)",
    "Visium Space Ranger\n+ Seurat",
    "H&E image",
]
rows = []
for r in regions:
    n_runs = int((vis_meta.patient_region_id == r).sum())
    rows.append(
        {
            "region": r,
            "patient": patient_of.get(r, "?"),
            "zone": zone_of.get(r, "?"),
            avail_cols[0]: int(r in set(rna_meta.patient_region_id)),
            avail_cols[1]: int(r in set(atac_meta.patient_region_id)),
            avail_cols[2]: int(region_key(r) in h5ad_info),
            avail_cols[3]: int(sec.patient_region_id.eq(r).sum() > 0),
            "H&E image": int(r in set(he_meta.patient_region_id)),
            "visium_runs": n_runs,
        }
    )
avail = pd.DataFrame(rows)
avail["zone"] = pd.Categorical(
    avail["zone"].map(lambda z: z if z in ZONES else "IZ"), ZONES, ordered=True
)
avail["pnum"] = avail["patient"].str.removeprefix("P").astype(int)
avail = avail.sort_values(["zone", "pnum", "region"]).reset_index(drop=True)
avail.drop(columns="pnum").to_csv(OUT / "modality_availability.csv", index=False)
triple = avail[
    (avail[avail_cols[0]] == 1)
    & (avail[avail_cols[1]] == 1)
    & (avail[avail_cols[3]] == 1)
]
stats["regions_total"] = len(avail)
stats["regions_triple"] = list(triple.region)
stats["regions_visium"] = int((avail[avail_cols[3]] == 1).sum())
stats["regions_rna_only"] = list(
    avail[(avail[avail_cols[0]] == 1) & (avail[avail_cols[3]] == 0)].region
)

fig, ax = plt.subplots(figsize=(A4W, 0.17 * len(avail) + 1.2))
mat = avail[avail_cols].to_numpy()
ax.imshow(
    mat,
    cmap=mpl.colors.ListedColormap(["#F2F2F2", "#404040"]),
    aspect="auto",
    vmin=0,
    vmax=1,
)
ax.set_xticks(range(len(avail_cols)), avail_cols)
ax.xaxis.tick_top()
ax.set_yticks(
    range(len(avail)),
    [
        f"{r}  ({p})" + ("  [2 Visium runs]" if n > 1 else "")
        for r, p, n in zip(avail.region, avail.patient, avail.visium_runs)
    ],
)
for t, z in zip(ax.get_yticklabels(), avail.zone):
    t.set_color(ZONE_COL[z])
n_rna = rna_meta.groupby("patient_region_id")["n_cells"].sum()
n_atac = atac_meta.groupby("patient_region_id")["n_cells"].sum()
n_spots = sec.groupby("patient_region_id")["n_spots"].sum()
h5_spots = {
    r: h5ad_info[region_key(r)][0] for r in avail.region if region_key(r) in h5ad_info
}
for i, r in enumerate(avail.region):
    for j, val in enumerate(
        [n_rna.get(r), n_atac.get(r), h5_spots.get(r), n_spots.get(r), None]
    ):
        if mat[i, j] and val is not None and val == val:
            ax.text(
                j,
                i,
                f"{int(val):,}",
                ha="center",
                va="center",
                fontsize=5.5,
                color="white",
            )
avail["n_nuclei_snRNA"] = avail.region.map(n_rna)
avail["n_nuclei_snATAC"] = avail.region.map(n_atac)
avail["n_spots_visium"] = avail.region.map(n_spots)
avail.drop(columns="pnum").to_csv(OUT / "modality_availability.csv", index=False)
stats["nuclei_rna_by_zone_median"] = (
    avail.groupby("zone", observed=True)["n_nuclei_snRNA"].median().to_dict()
)
stats["nuclei_atac_by_zone_median"] = (
    avail.groupby("zone", observed=True)["n_nuclei_snATAC"].median().to_dict()
)
ax.set_xlim(-0.5, len(avail_cols) - 0.5)
ax.tick_params(length=0)
for s in ax.spines.values():
    s.set_visible(False)
ax.legend(
    handles=[Patch(color=ZONE_COL[z], label=ZONE_NAME[z]) for z in ZONES],
    loc="upper left",
    bbox_to_anchor=(1.12, 1),
    frameon=False,
    title="Zone (label colour)",
)
fig.savefig(FIG / "fig1_availability.pdf", bbox_inches="tight")
plt.close(fig)

# ------------------------------------------------------- 2. cohort design
fig, ax = plt.subplots(figsize=(A4W, 2.6))
pats = sorted(avail.patient.unique(), key=lambda p: int(p[1:]))
for _, r in avail.iterrows():
    x, y = pats.index(r.patient), ZONES.index(r.zone)
    vis = r[avail_cols[3]] == 1
    ax.scatter(
        x,
        y + np.random.default_rng(len(r.region)).uniform(-0.18, 0.18),
        s=38 if vis else 22,
        color=ZONE_COL[r.zone] if vis else "white",
        edgecolor=ZONE_COL[r.zone],
        lw=1,
        zorder=3,
    )
multi = avail[avail[avail_cols[3]] == 1].groupby("patient").size()
multi = multi[multi > 1]
for p in multi.index:
    ax.axvspan(pats.index(p) - 0.45, pats.index(p) + 0.45, color="#EEEEEE", zorder=0)
ax.set_xticks(range(len(pats)), pats, rotation=90)
ax.set_yticks(range(len(ZONES)), [ZONE_NAME[z] for z in ZONES])
ax.set_xlabel("Patient")
ax.legend(
    handles=[
        plt.scatter([], [], s=38, color="grey", label="Visium section"),
        plt.scatter(
            [], [], s=22, color="white", edgecolor="grey", label="single-nucleus only"
        ),
        Patch(color="#EEEEEE", label="patient with >1 Visium section"),
    ],
    loc="upper left",
    bbox_to_anchor=(1.0, 1),
    frameon=False,
)
fig.savefig(FIG / "fig2_cohort.pdf", bbox_inches="tight")
plt.close(fig)
stats["patients_multi_visium"] = {p: int(n) for p, n in multi.items()}
stats["n_patients_visium"] = int(sec.patient.nunique())

# ----------------------------------------------------- 3. QC per section
fig, axes = plt.subplots(
    1,
    4,
    figsize=(A4W, 0.16 * len(sec) + 1.0),
    sharey=True,
    gridspec_kw={"width_ratios": [1, 1.4, 1.4, 1.4]},
)
ypos = np.arange(len(sec))
cols = [ZONE_COL[z] for z in sec.zone]
axes[0].barh(ypos, sec.n_spots, color=cols)
axes[0].set_xlabel("Spots on tissue")
for ax, col, lab, log in [
    (axes[1], "nCount_Spatial", "UMIs per spot", True),
    (axes[2], "nFeature_Spatial", "Genes per spot", True),
    (axes[3], "percent.mt", "% mitochondrial UMIs", False),
]:
    data = [spots.loc[spots.hca_sample_id == s, col].to_numpy() for s in ORDER]
    bp = ax.boxplot(
        data,
        positions=ypos,
        vert=False,
        widths=0.65,
        showfliers=False,
        patch_artist=True,
        medianprops={"color": "black", "lw": 0.8},
        whiskerprops={"lw": 0.6},
        capprops={"lw": 0.6},
        boxprops={"lw": 0.6},
    )
    for b, c in zip(bp["boxes"], cols):
        b.set_facecolor(c)
    ax.set_xlabel(lab)
    if log:
        ax.set_xscale("log")
axes[0].set_yticks(ypos, [LABEL[s] for s in ORDER])
axes[0].invert_yaxis()
axes[3].legend(
    handles=[Patch(color=ZONE_COL[z], label=ZONE_NAME[z]) for z in ZONES],
    loc="upper left",
    bbox_to_anchor=(1.02, 1),
    frameon=False,
)
fig.savefig(FIG / "fig3_qc_per_section.pdf", bbox_inches="tight")
plt.close(fig)

# ------------------------------------------- 4. spot depth and cellularity by zone
fig, axes = plt.subplots(2, 2, figsize=(A4W, 5.2))
for ax, col, lab, log in [
    (axes[0, 0], "nCount_Spatial", "UMIs per spot", True),
    (axes[0, 1], "nFeature_Spatial", "Genes per spot", True),
    (
        axes[1, 0],
        "c2l_total_abundance",
        "Estimated cells per spot\n(cell2location, shipped)",
        False,
    ),
]:
    data = [spots.loc[spots.zone == z, col].to_numpy() for z in ZONES]
    vp = ax.violinplot(
        [np.log10(d) if log else d for d in data], showmedians=True, showextrema=False
    )
    for b, z in zip(vp["bodies"], ZONES):
        b.set_facecolor(ZONE_COL[z])
        b.set_alpha(0.8)
    vp["cmedians"].set_color("black")
    ax.set_xticks(range(1, 6), [ZONE_NAME[z] for z in ZONES])
    ax.set_ylabel(("log10 " if log else "") + lab)
    if not log:
        ax.set_ylim(
            0, np.quantile(spots[col], 0.995)
        )  # long tail clipped at the 99.5th percentile
ax = axes[1, 1]
spots["umi_per_cell"] = spots.nCount_Spatial / spots.c2l_total_abundance
for _, s in sec.iterrows():
    ax.scatter(s.median_cells, s.median_umi, color=ZONE_COL[s.zone], s=18, zorder=3)
ax.set_xlabel("Median estimated cells per spot (per section)")
ax.set_ylabel("Median UMIs per spot (per section)")
ax.legend(
    handles=[Patch(color=ZONE_COL[z], label=ZONE_NAME[z]) for z in ZONES], frameon=False
)
fig.tight_layout()
fig.savefig(FIG / "fig4_depth_by_zone.pdf", bbox_inches="tight")
plt.close(fig)
by_zone = spots.groupby("zone", observed=True).agg(
    n_sections=("hca_sample_id", "nunique"),
    n_spots=("barcode", "size"),
    median_umi=("nCount_Spatial", "median"),
    median_genes=("nFeature_Spatial", "median"),
    median_cells=("c2l_total_abundance", "median"),
    median_umi_per_cell=("umi_per_cell", "median"),
    median_mt=("percent.mt", "median"),
)
by_zone.to_csv(OUT / "depth_by_zone.csv")
stats["by_zone"] = (
    by_zone.round(2).reset_index().astype({"zone": str}).to_dict(orient="records")
)
stats["rho_cells_umi_section"] = float(
    sec[["median_cells", "median_umi"]].corr(method="spearman").iloc[0, 1]
)
stats["rho_cells_umi_spot"] = float(
    spots[["c2l_total_abundance", "nCount_Spatial"]].corr(method="spearman").iloc[0, 1]
)

# ------------------------------------------------ 5. published niches, spatially
ct_levels = sorted(
    spots.celltype_niche.dropna().unique(), key=lambda s: int(s.split("_")[1])
)
mol_levels = sorted(
    spots.molecular_niche.dropna().unique(), key=lambda s: int(s.split("_")[1])
)
ct_col = dict(
    zip(ct_levels, mpl.colormaps["tab10"].colors + mpl.colormaps["Dark2"].colors)
)
mol_col = dict(zip(mol_levels, mpl.colormaps["tab20"].colors))


def lowres(sid):
    return np.asarray(Image.open(TABLES / f"{sid}_lowres.png").convert("RGB"))


def spot_xy(df, sid):
    f = scalef[sid]["tissue_lowres_scalef"]
    return df.pxl_col_fullres.to_numpy() * f, df.pxl_row_fullres.to_numpy() * f


def spot_size(sid, ax_px):
    return (
        0.8 * scalef[sid]["spot_diameter_fullres"] * scalef[sid]["tissue_lowres_scalef"]
    )


# representative section per zone: the one with most spots
rep = {
    z: sec[sec.zone == z].sort_values("n_spots").iloc[-1].hca_sample_id
    for z in ZONES
    if (sec.zone == z).any()
}
stats["representative"] = {z: LABEL[s] for z, s in rep.items()}
fig, axes = plt.subplots(len(rep), 3, figsize=(A4W, 2.2 * len(rep)))
for i, (z, sid) in enumerate(rep.items()):
    df = spots[spots.hca_sample_id == sid]
    img = lowres(sid)
    x, y = spot_xy(df, sid)
    pad = 15
    lim = (x.min() - pad, x.max() + pad, y.max() + pad, y.min() - pad)
    for j, (title, colour) in enumerate(
        [
            ("H&E", None),
            ("Cell-type niche", df.celltype_niche.map(ct_col)),
            ("Molecular niche", df.molecular_niche.map(mol_col)),
        ]
    ):
        ax = axes[i, j]
        ax.imshow(img if colour is None else np.full_like(img, 255))
        if colour is not None:
            ax.scatter(x, y, c=list(colour), s=2.2, linewidths=0)
        ax.set_xlim(lim[0], lim[1])
        ax.set_ylim(lim[2], lim[3])
        ax.set_xticks([])
        ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)
        if i == 0:
            ax.set_title(title)
    axes[i, 0].set_ylabel(f"{ZONE_NAME[z]}\n{LABEL[sid]}", fontsize=7)
fig.legend(
    handles=[Patch(color=ct_col[k], label=k) for k in ct_levels],
    loc="upper center",
    bbox_to_anchor=(0.5, 0.06),
    ncol=9,
    frameon=False,
    title="Cell-type niches",
)
fig.legend(
    handles=[Patch(color=mol_col[k], label=k) for k in mol_levels],
    loc="upper center",
    bbox_to_anchor=(0.5, 0.025),
    ncol=12,
    frameon=False,
    title="Molecular niches",
)
fig.savefig(FIG / "fig5_niches_representative.pdf", bbox_inches="tight", dpi=300)
plt.close(fig)

# 5b. spatial coherence of each niche label: share of a spot's grid neighbours with the
# same label, against the same share after shuffling labels within the section.
HEX = [(0, 2), (0, -2), (1, 1), (1, -1), (-1, 1), (-1, -1)]
rng = np.random.default_rng(0)


def same_label_share(df, col):
    key = {(r, c): lab for r, c, lab in zip(df.array_row, df.array_col, df[col])}
    same = tot = 0
    for (r, c), lab in key.items():
        for dr, dc in HEX:
            n = key.get((r + dr, c + dc))
            if n is not None:
                tot += 1
                same += n == lab
    return same / tot if tot else np.nan


coh = []
for sid in ORDER:
    df = spots[spots.hca_sample_id == sid]
    for col in ["celltype_niche", "molecular_niche"]:
        obs = same_label_share(df, col)
        shuf = np.mean(
            [
                same_label_share(
                    df.assign(**{col: rng.permutation(df[col].to_numpy())}), col
                )
                for _ in range(20)
            ]
        )
        coh.append(
            {"hca_sample_id": sid, "annotation": col, "observed": obs, "shuffled": shuf}
        )
coh = pd.DataFrame(coh).merge(sec[["hca_sample_id", "zone", "label"]])
coh.to_csv(OUT / "niche_coherence.csv", index=False)
fig, axes = plt.subplots(1, 2, figsize=(A4W, 0.16 * len(sec) + 1.0), sharey=True)
for ax, col, title in [
    (axes[0], "celltype_niche", "Cell-type niche"),
    (axes[1], "molecular_niche", "Molecular niche"),
]:
    d = coh[coh.annotation == col].set_index("hca_sample_id").loc[ORDER]
    ax.hlines(ypos, d.shuffled, d.observed, color="#BBBBBB", lw=1)
    ax.scatter(
        d.shuffled,
        ypos,
        color="white",
        edgecolor="#555555",
        s=14,
        zorder=3,
        label="labels shuffled",
    )
    ax.scatter(
        d.observed,
        ypos,
        color=[ZONE_COL[z] for z in d.zone],
        s=18,
        zorder=3,
        label="observed",
    )
    ax.set_xlim(0, 1)
    ax.set_xlabel("Share of neighbouring spots with the same label")
    ax.set_title(title)
axes[0].set_yticks(ypos, [LABEL[s] for s in ORDER])
axes[0].invert_yaxis()
axes[1].legend(loc="lower right", frameon=False)
fig.savefig(FIG / "fig5b_niche_coherence.pdf", bbox_inches="tight")
plt.close(fig)
stats["coherence_median"] = (
    coh.groupby("annotation")[["observed", "shuffled"]].median().round(2).to_dict()
)

# ---------------------------------------------------------- 6. composition
comp = spots.groupby("hca_sample_id")[prop_cols].mean().loc[ORDER]
comp.to_csv(OUT / "composition_per_section.csv")
ct_names = [c.removeprefix("prop_") for c in prop_cols]
pal = dict(
    zip(
        ct_names,
        mpl.colormaps["tab20"].colors[::2] + mpl.colormaps["tab20"].colors[1::2],
    )
)
fig, ax = plt.subplots(figsize=(A4W, 0.16 * len(sec) + 1.0))
left = np.zeros(len(comp))
for c, n in zip(prop_cols, ct_names):
    ax.barh(ypos, comp[c], left=left, color=pal[n], label=n, height=0.8)
    left += comp[c].to_numpy()
ax.set_yticks(ypos, [LABEL[s] for s in ORDER])
for t, z in zip(ax.get_yticklabels(), sec.zone):
    t.set_color(ZONE_COL[z])
ax.invert_yaxis()
ax.set_xlim(0, 1)
ax.set_xlabel("Mean cell2location proportion across spots (shipped)")
ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1), frameon=False)
fig.savefig(FIG / "fig6_composition.pdf", bbox_inches="tight")
plt.close(fig)
stats["cm_fraction_by_zone"] = (
    comp.join(sec.set_index("hca_sample_id").zone)
    .groupby("zone", observed=True)["prop_CM"]
    .median()
    .round(2)
    .to_dict()
)
stats["fib_fraction_by_zone"] = (
    comp.join(sec.set_index("hca_sample_id").zone)
    .groupby("zone", observed=True)["prop_Fib"]
    .median()
    .round(2)
    .to_dict()
)

# --------------------------------------------------------- 7. H&E inventory
ncol = 6
nrow = int(np.ceil(len(ORDER) / ncol))
fig, axes = plt.subplots(nrow, ncol, figsize=(A4W, 1.35 * nrow))
for ax in axes.flat:
    ax.axis("off")
for ax, sid in zip(axes.flat, ORDER):
    ax.imshow(lowres(sid))
    z = sec.set_index("hca_sample_id").zone[sid]
    ax.set_title(LABEL[sid].replace(" (", "\n("), fontsize=5.5, color=ZONE_COL[z])
fig.savefig(FIG / "fig7_he_thumbnails.pdf", bbox_inches="tight", dpi=200)
plt.close(fig)

he_files = sorted((DATA / "he_images").glob("*.tif"))
he_tab = []
for f in he_files:
    with Image.open(f) as im:
        he_tab.append({"file": f.name, "width": im.size[0], "height": im.size[1]})
he_tab = pd.DataFrame(he_tab)
he_meta["file"] = he_meta.sample_id.map(
    lambda s: (
        f"{s.split('_')[1]}_2019{s.split('_')[0]}.tif"
        if s.startswith("AKK")
        else f"{s}.tif"
    )
)
he_tab = he_meta.merge(he_tab, on="file", how="outer")
he_tab.to_csv(OUT / "he_inventory.csv", index=False)
stats["he_images"] = len(he_files)
stats["he_unmatched"] = int(
    he_tab.width.isna().sum() + he_tab.patient_region_id.isna().sum()
)
stats["scalefactors"] = len(scalef)

# ---------------------------------------------------------- section table
sec["h5ad"] = sec.patient_region_id.map(
    lambda r: "yes" if region_key(r) in h5ad_info else "no"
)
sec["h5ad_spots"] = sec.patient_region_id.map(
    lambda r: h5ad_info.get(region_key(r), (np.nan,))[0]
)
sec["he_image"] = sec.patient_region_id.isin(set(he_meta.patient_region_id)).map(
    {True: "yes", False: "no"}
)
sec = sec.merge(
    metrics[["hca_sample_id", "Mean Reads per Spot", "Sequencing Saturation"]],
    how="left",
)
sec.to_csv(OUT / "visium_sections.csv", index=False)
stats["rz_p6_in_h5ad_spots"] = h5ad_info.get("RZ_P6", (None,))[0]
(OUT / "stats.json").write_text(json.dumps(stats, indent=1, default=str))
print(json.dumps(stats, indent=1, default=str))
