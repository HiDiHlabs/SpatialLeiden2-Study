"""Cohort overview figure for Kuppe et al. 2022 (Nature 608, 766-777).

One row per patient (23). Left: which zones each patient contributed and which
data types exist for each sample. Right: days after symptom onset. Descriptive
only; values come from the paper's Supplementary Table 1 (clinical data), the
authors' snRNA/snATAC sample sheets (Zenodo 6578617) and the Visium sample sheet
(Zenodo 6580069).

Supplementary Table 1 is not on Zenodo. Download it once into
data/kuppe_2022/metadata/ as Kuppe2022_SupplementaryTable1_clinical.xlsx from
https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-022-05060-x/MediaObjects/41586_2022_5060_MOESM6_ESM.xlsx

Usage (from Derived_Modalities/):
    python 01_data_preprocess/kuppe_cohort/cohort_figure.py data/kuppe_2022 results/kuppe_2022/cohort
Writes kuppe_cohort.pdf, kuppe_cohort.png and cohort_table.csv.
"""

import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D

DATA, OUT = Path(sys.argv[1]), Path(sys.argv[2])
OUT.mkdir(parents=True, exist_ok=True)
MM = 1 / 25.4

mpl.rcParams.update(
    {
        "font.family": "Arial",
        "font.size": 6,
        "axes.labelsize": 6,
        "axes.titlesize": 6,
        "xtick.labelsize": 5.5,
        "ytick.labelsize": 6,
        "legend.fontsize": 5.5,
        "axes.linewidth": 0.5,
        "xtick.major.width": 0.5,
        "ytick.major.width": 0.5,
        "xtick.minor.width": 0.4,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)

# Zones: fixed order and colourblind-safe colours (Wong / Okabe-Ito), as in the exploration report.
ZONES = ["CTRL", "RZ", "BZ", "IZ", "FZ"]
ZONE_LABEL = {
    "CTRL": "Control",
    "RZ": "Remote",
    "BZ": "Border",
    "IZ": "Ischaemic",
    "FZ": "Fibrotic",
}
ZONE_COL = {
    "CTRL": "#0072B2",
    "RZ": "#009E73",
    "BZ": "#E69F00",
    "IZ": "#D55E00",
    "FZ": "#CC79A7",
}

# ------------------------------------------------------------------ data
st1 = pd.read_excel(
    DATA / "metadata" / "Kuppe2022_SupplementaryTable1_clinical.xlsx", header=1
)
st1.columns = st1.columns.str.strip()
st1 = st1[st1["patient"].astype(str).str.startswith("P")].copy()
st1["region_zone_description"] = st1["region_zone_description"].str.strip()
st1["zone"] = st1["region_zone_description"]
st1["days"] = pd.to_numeric(st1["days after infarction"], errors="coerce")


def sheet(name):
    x = pd.read_excel(DATA / "metadata" / name, header=1)
    return set(x.loc[x["sample_id"].notna(), "patient_region_id"])


rna, atac = sheet("metadata-RNA-seq.xlsx"), sheet("metadata-ATAC-seq.xlsx")
vis = pd.read_csv(DATA / "visium_raw" / "metadata-Visium.csv")
vis_runs = vis.groupby("patient_region_id")["hca_sample_id"].apply(
    lambda s: ";".join(sorted(s))
)


def group(r):
    if r.disease == "control":
        return "Control donor"
    if r.patient_region_id.endswith(("_ext1", "_ext2", "_ext3")):
        return "Validation (snRNA-seq only)"
    if r.disease == "AMI":
        return "Acute MI"
    return "Transplant recipient"


st1["group"] = st1.apply(group, axis=1)
st1["snRNA"] = st1.patient_region_id.isin(rna) | st1.group.eq(
    "Validation (snRNA-seq only)"
)
st1["snATAC"] = st1.patient_region_id.isin(atac)
st1["visium"] = st1.patient_region_id.isin(vis_runs.index)
st1["visium_runs"] = st1.patient_region_id.map(vis_runs)

# per-sample Visium depth from the exploration step, if it has been run
sec_file = (
    DATA.parent.parent
    / "results"
    / "kuppe_2022"
    / "visium_exploration"
    / "visium_sections.csv"
)
if sec_file.exists():
    sec = (
        pd.read_csv(sec_file)
        .groupby("patient_region_id")
        .agg(
            n_spots=("n_spots", "sum"),
            median_counts_per_spot=("median_umi", "median"),
            median_genes_per_spot=("median_genes", "median"),
        )
    )
    st1 = st1.join(sec, on="patient_region_id")

table = st1.rename(
    columns={
        "patient": "donor",
        "region_novel": "region",
        "age (years)": "age",
        "gender": "sex",
        "days": "days_after_onset",
    }
)
table["modalities"] = table.apply(
    lambda r: "+".join(
        m
        for m, k in [
            ("snRNA-seq", "snRNA"),
            ("snATAC-seq", "snATAC"),
            ("Visium", "visium"),
        ]
        if r[k]
    ),
    axis=1,
)
keep = [
    "donor",
    "patient_region_id",
    "region",
    "zone",
    "group",
    "disease",
    "age",
    "sex",
    "days_after_onset",
    "modalities",
    "visium_runs",
    "n_spots",
    "median_counts_per_spot",
    "median_genes_per_spot",
]
table[[c for c in keep if c in table.columns]].to_csv(
    OUT / "cohort_table.csv", index=False
)

# ------------------------------------------------------------- row order
GROUPS = [
    "Control donor",
    "Acute MI",
    "Transplant recipient",
    "Validation (snRNA-seq only)",
]
donors = (
    st1.groupby("patient")
    .agg(group=("group", "first"), days=("days", "max"))
    .reset_index()
)
donors["g"] = donors.group.map(GROUPS.index)
donors["pnum"] = donors.patient.str[1:].astype(int)
donors = donors.sort_values(["g", "days", "pnum"]).reset_index(drop=True)
row = {p: i for i, p in enumerate(donors.patient)}
n = len(donors)

# ------------------------------------------------------------------ figure
fig = plt.figure(figsize=(183 * MM, 96 * MM))
gs = fig.add_gridspec(
    1,
    2,
    width_ratios=[1.0, 1.05],
    wspace=0.08,
    left=0.20,
    right=0.80,
    top=0.90,
    bottom=0.17,
)
ax = fig.add_subplot(gs[0])
axt = fig.add_subplot(gs[1], sharey=ax)

# group bands, shared by both panels
for gi, gname in enumerate(GROUPS):
    rows = donors.index[donors.group == gname]
    if len(rows) == 0:
        continue
    lo, hi = rows.min() - 0.5, rows.max() + 0.5
    for a in (ax, axt):
        if gi % 2 == 0:
            a.axhspan(lo, hi, color="#F2F2F2", zorder=0, lw=0)
    ax.annotate(
        gname.replace(" (", "\n(")
        .replace("recipient", "\nrecipients")
        .replace("donor", "donors"),
        xy=(-0.19, (lo + hi) / 2),
        xycoords=("axes fraction", "data"),
        ha="right",
        va="center",
        fontsize=6,
        rotation=0,
        linespacing=1.1,
    )

# a: zones x data types
for (p, z), grp in st1.groupby(["patient", "zone"]):
    y, x = row[p], ZONES.index(z)
    full = grp.snRNA.all() and grp.snATAC.all() and grp.visium.all()
    ax.scatter(
        x,
        y,
        s=34,
        marker="o",
        zorder=3,
        linewidths=0.8,
        facecolor=ZONE_COL[z] if full else "white",
        edgecolor=ZONE_COL[z],
    )
    if len(grp) > 1:
        ax.text(x + 0.22, y, f"×{len(grp)}", va="center", ha="left", fontsize=5.5)
ax.set_xticks(range(len(ZONES)), [ZONE_LABEL[z] for z in ZONES])
ax.xaxis.tick_top()
ax.set_xlim(-0.6, len(ZONES) - 0.4)
ax.set_ylim(n - 0.5, -0.5)
ax.set_yticks(range(n), donors.patient)
ax.tick_params(length=0)
for s in ax.spines.values():
    s.set_visible(False)
# sample counts per zone, below the grid
counts = st1.groupby("zone").size().reindex(ZONES)
vcounts = st1[st1.visium].groupby("zone").size().reindex(ZONES).fillna(0).astype(int)
for i, z in enumerate(ZONES):
    ax.text(i, n + 0.4, f"{counts[z]}", ha="center", va="top", fontsize=6)
    ax.text(i, n + 1.5, f"{vcounts[z]}", ha="center", va="top", fontsize=6)
ax.text(-0.75, n + 0.4, "Samples", ha="right", va="top", fontsize=6)
ax.text(-0.75, n + 1.5, "with Visium", ha="right", va="top", fontsize=6)

# b: days after symptom onset (log scale)
d = donors.dropna(subset=["days"])
axt.scatter(d.days, d.index, s=14, color="#333333", zorder=3, linewidths=0)
for i, r in d.iterrows():
    axt.text(r.days * 1.12, i, f"{int(r.days)} d", va="center", fontsize=5.5)
for i, r in donors[donors.days.isna()].iterrows():
    axt.text(1.2, i, "n.a.", va="center", fontsize=5.5, color="#777777")
axt.set_xscale("log")
axt.set_xlim(1, 400)
axt.set_xticks(
    [1, 2, 5, 10, 20, 50, 100, 200], ["1", "2", "5", "10", "20", "50", "100", "200"]
)
axt.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
axt.axvline(15, color="#999999", lw=0.5, ls=(0, (2, 2)), zorder=1)
axt.text(14, -0.9, "acute", ha="right", va="bottom", fontsize=5.5, color="#555555")
axt.text(16, -0.9, "later stage", ha="left", va="bottom", fontsize=5.5, color="#555555")
axt.set_xlabel("Days after symptom onset")
axt.tick_params(axis="y", left=False, labelleft=False)
for s in ("top", "right", "left"):
    axt.spines[s].set_visible(False)

# legend, outside the panels
handles = [
    Line2D(
        [],
        [],
        ls="",
        marker="o",
        ms=5,
        mfc="#555555",
        mec="#555555",
        mew=0.8,
        label="snRNA-seq + snATAC-seq\n+ Visium",
    ),
    Line2D(
        [],
        [],
        ls="",
        marker="o",
        ms=5,
        mfc="white",
        mec="#555555",
        mew=0.8,
        label="snRNA-seq only",
    ),
    Line2D([], [], ls="", marker="", label="×2  two samples\n      of one zone"),
]
fig.legend(
    handles=handles,
    loc="upper left",
    bbox_to_anchor=(0.815, 0.90),
    frameon=False,
    title="Data per sample",
    title_fontsize=6,
    alignment="left",
    handletextpad=0.4,
    labelspacing=0.9,
)

fig.text(0.008, 0.965, "a", fontsize=8, fontweight="bold", va="top")
fig.text(0.505, 0.965, "b", fontsize=8, fontweight="bold", va="top")

fig.savefig(OUT / "kuppe_cohort.pdf")
fig.savefig(OUT / "kuppe_cohort.png", dpi=450)
print(
    f"{len(st1)} samples, {n} donors;", counts.to_dict(), "| Visium:", vcounts.to_dict()
)
print(donors.group.value_counts().to_dict())
