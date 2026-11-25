"""Which joint SpatialLeiden domains are enriched in control vs patient tissue?

Domains: one joint clustering of all 27 slides (default RNA + space, weight 0.75, seed 0;
results/kuppe_2022/exp2_joint_harmony, column r0.75_s0). The unit of analysis is the slide
(n = 27), not the spot: per slide, the share of its spots in each domain.
Groups from the authors' spot tables: zone (major_labl: CTRL, RZ, BZ, IZ, FZ) and patient group
(patient_group: myogenic, ischaemic, fibrotic). Control = CTRL slides.
Tests (per domain, on the slide shares): Kruskal-Wallis across zones; one-sided Mann-Whitney
"higher in X than in CTRL" for each patient zone and "higher in CTRL than in all patient
slides"; BH correction over all tests of a family. With 4 control slides, p-values are coarse:
read the effect sizes (mean share per zone) first.
Domain identity: median cell2location proportions of each domain's spots (z-scored across domains).

Usage (from Derived_Modalities/):
    python 03_single_modality_clustering/exp2_joint_harmony/domain_enrichment_by_zone.py data/kuppe_2022 results/kuppe_2022 \
        [labels_folder=exp2_joint_harmony] [column=r0.75_s0]
Output: results/kuppe_2022/joint_domain_enrichment/<folder>_<column>_{shares.csv,tests.csv,fig.png,fig.pdf}
"""

import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import kruskal, mannwhitneyu
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "common"))
from plot_patient_region_umap import PATIENT_COLOURS, REGION_COLOURS  # noqa: E402  authors' palettes (ED Fig. 3a-b)

DATA, RES = Path(sys.argv[1]), Path(sys.argv[2])
FOLDER = sys.argv[3] if len(sys.argv) > 3 else "exp2_joint_harmony"
COL = sys.argv[4] if len(sys.argv) > 4 else "r0.75_s0"
OUT = RES / "joint_domain_enrichment"
OUT.mkdir(parents=True, exist_ok=True)
TAG = f"{FOLDER}_{COL}"
MM = 1 / 25.4
mpl.rcParams.update({"font.family": "Arial", "font.size": 6, "axes.titlesize": 6.5, "pdf.fonttype": 42,
                     "axes.spines.top": False, "axes.spines.right": False})
ZONES = ["CTRL", "RZ", "BZ", "IZ", "FZ"]
ZONE_COL = {"CTRL": REGION_COLOURS["control"], **{z: REGION_COLOURS[z] for z in ["RZ", "BZ", "IZ", "FZ"]}}  # authors' palette

lab = pd.read_csv(RES / FOLDER / "exp2_labels.csv.gz", index_col=0, dtype=str)[COL]
run = lab.index.str.split("_", n=1).str[0]
bc = lab.index.str.split("_", n=1).str[1]
spots = []
for r in pd.unique(run):
    t = pd.read_csv(DATA / "visium_tables" / f"{r}_spots.csv").set_index("barcode")
    t.index = r + "_" + t.index
    t["run"] = r
    spots.append(t)
spots = pd.concat(spots).loc[lab.index]
spots["domain"] = lab.to_numpy()
# domains numbered D1.. by total size
order = spots.domain.value_counts().index
name = {d: f"D{i + 1}" for i, d in enumerate(order)}
spots["domain"] = spots.domain.map(name)
doms = [name[d] for d in order]

meta = spots.groupby("run")[["major_labl", "patient_group", "patient", "patient_region_id"]].first()
shares = pd.crosstab(spots.run, spots.domain, normalize="index")[doms]
shares = meta.join(shares)
shares.to_csv(OUT / f"{TAG}_shares.csv")

tests = []
ctrl = shares.major_labl == "CTRL"
for d in doms:
    groups = [shares.loc[shares.major_labl == z, d].to_numpy() for z in ZONES if (shares.major_labl == z).any()]
    kw = kruskal(*groups).pvalue if all(len(g) > 1 for g in groups) else np.nan
    row = {"domain": d, "n_spots": int((spots.domain == d).sum()), "kruskal_zone_p": kw}
    for z in ZONES:
        row[f"mean_share_{z}"] = shares.loc[shares.major_labl == z, d].mean()
    for z in ZONES[1:]:
        row[f"p_{z}_gt_CTRL"] = mannwhitneyu(shares.loc[shares.major_labl == z, d], shares.loc[ctrl, d],
                                             alternative="greater").pvalue
    row["p_CTRL_gt_patients"] = mannwhitneyu(shares.loc[ctrl, d], shares.loc[~ctrl, d], alternative="greater").pvalue
    row["log2_patient_over_ctrl"] = np.log2((shares.loc[~ctrl, d].mean() + 1e-3) / (shares.loc[ctrl, d].mean() + 1e-3))
    for g in ["myogenic", "ischaemic", "fibrotic"]:
        row[f"mean_share_{g}"] = shares.loc[shares.patient_group == g, d].mean()
    tests.append(row)
tests = pd.DataFrame(tests)
pcols = [c for c in tests.columns if c.startswith("p_")]
adj = multipletests(tests[pcols].to_numpy().ravel(), method="fdr_bh")[1].reshape(len(tests), len(pcols))
for i, c in enumerate(pcols):
    tests[c.replace("p_", "padj_")] = adj[:, i]
tests["padj_kruskal_zone"] = multipletests(tests.kruskal_zone_p.fillna(1), method="fdr_bh")[1]
tests.to_csv(OUT / f"{TAG}_tests.csv", index=False)

# composition identity of each domain
props = [c for c in spots.columns if c.startswith("prop_")]
comp = spots.groupby("domain")[props].median().loc[doms]
compz = (comp - comp.mean()) / comp.std(ddof=0).replace(0, 1)
top = {d: ", ".join(c.replace("prop_", "") for c in comp.loc[d].sort_values(ascending=False).index[:2]) for d in doms}

# figure: a mean share per zone; b per-slide shares; c composition
fig = plt.figure(figsize=(183 * MM, 95 * MM))
gs = fig.add_gridspec(1, 3, width_ratios=[0.9, 2.2, 1.1], wspace=0.55)
ax = fig.add_subplot(gs[0])
M = tests.set_index("domain")[[f"mean_share_{z}" for z in ZONES]].to_numpy()
im = ax.imshow(M, cmap="Blues", aspect="auto", vmin=0)
for (i, j), v in np.ndenumerate(M):
    ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=4.8, color="white" if v > 0.6 * M.max() else "black")
ax.set_xticks(range(len(ZONES)), ZONES)
ax.set_yticks(range(len(doms)), [f"{d}  ({top[d]})" for d in doms])
ax.set_title("a  Mean share of the slide", fontweight="bold", loc="left")
cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.08, location="bottom")
cb.ax.tick_params(labelsize=5)
ax2 = fig.add_subplot(gs[1])
rng = np.random.default_rng(0)
for i, d in enumerate(doms):
    for k, z in enumerate(ZONES):
        v = shares.loc[shares.major_labl == z, d].to_numpy()
        x = i + (k - 2) * 0.16
        ax2.scatter(x + rng.uniform(-0.03, 0.03, len(v)), v, s=5, color=ZONE_COL[z], linewidths=0)
        ax2.plot([x - 0.06, x + 0.06], [v.mean()] * 2, color=ZONE_COL[z], lw=1)
    t = tests.set_index("domain").loc[d]
    mark = "▲ patients" if t.padj_CTRL_gt_patients > 0.1 and min(t[f"padj_{z}_gt_CTRL"] for z in ZONES[1:]) < 0.1 else (
        "▲ control" if t.padj_CTRL_gt_patients < 0.1 else "")
    ax2.text(i, 0.97, mark.replace(" ", "\n"), ha="center", va="top", fontsize=5, transform=ax2.get_xaxis_transform())
ax2.set_xticks(range(len(doms)), doms)
ax2.set_ylabel("share of the slide's spots")
ax2.legend(handles=[mpl.lines.Line2D([], [], marker="o", ls="", color=ZONE_COL[z], label=z, markersize=3) for z in ZONES],
           frameon=False, fontsize=5, ncol=5, loc="upper center", bbox_to_anchor=(0.5, -0.08))
ax2.set_title("b  Per slide (dots) and zone mean (bar); ▲ = BH-adjusted p < 0.1", fontweight="bold", loc="left")
ax3 = fig.add_subplot(gs[2])
ax3.imshow(compz.to_numpy().clip(-2.5, 2.5), cmap="RdBu_r", vmin=-2.5, vmax=2.5, aspect="auto")
ax3.set_xticks(range(len(props)), [p.replace("prop_", "") for p in props], rotation=90)
ax3.set_yticks(range(len(doms)), doms)
ax3.set_title("c  Cell-type composition\n(median, scaled across domains)", fontweight="bold", loc="left")
for a in (ax, ax2, ax3):
    a.tick_params(labelsize=5.5, length=1.5, width=0.4)
fig.suptitle(f"Joint SpatialLeiden domains ({FOLDER}, {COL}) across 27 slides", fontsize=7, y=1.02)
fig.savefig(OUT / f"{TAG}_fig.png", dpi=250, bbox_inches="tight")
fig.savefig(OUT / f"{TAG}_fig.pdf", dpi=300, bbox_inches="tight")
pd.set_option("display.width", 220)
show = ["domain", "n_spots"] + [f"mean_share_{z}" for z in ZONES] + ["log2_patient_over_ctrl", "padj_CTRL_gt_patients"] + \
       [f"padj_{z}_gt_CTRL" for z in ZONES[1:]] + ["padj_kruskal_zone"]
t = tests[show].copy()
t.insert(1, "top_celltypes", t.domain.map(top))
print(shares.major_labl.value_counts().to_dict())
print(t.round(3).to_string(index=False))
