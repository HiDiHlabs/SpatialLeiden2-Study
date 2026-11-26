"""Zone enrichment beyond the per-domain permutation test: scCODA (compositional) and Milo (clustering-free).

scCODA (pertpy.tl.Sccoda, Bayesian compositional model, slide = sample, 27 slides):
  counts of spots per domain per slide for the joint clusterings RNA + space (exp2_joint_harmony r0.75_s0),
  RNA + TF + space (exp3_joint_rna_tf r0.75_s0) and RNA + CCC + space (exp3_joint_rna_ccc_ours r0.25_s0);
  formula C(zone, Treatment('CTRL')), reference domain chosen automatically, NUTS (seed 0), FDR 0.2.
  Output: credible effects (log fold change vs control) per domain and zone.
Milo (pertpy.tl.Milo, neighbourhood differential abundance):
  graph = 30-NN on the RNA Harmony embedding (exp1 "authors"); neighbourhoods from 10 % of spots (seed 0);
  counts per slide; one test per patient zone vs control (subset of slides, design ~zone), solver pydeseq2.
  Neighbourhoods are then annotated with the RNA + TF + space domains, with membership of the
  stressed-CM domain and of the grid-search region families (majority of the neighbourhood's spots), so the
  logFC / SpatialFDR distribution per domain or region shows whether the enrichment holds without clusters.
Caveats: 4 control and 3 BZ slides; patient is not modelled (several patients contribute several slides).

Usage (from Derived_Modalities/, env pertpy):
    python 08_grid_search/sccoda_milo.py data/kuppe_2022 results/kuppe_2022
Output: results/kuppe_2022/sccoda_milo/ (sccoda_effects.csv, milo_nhoods.csv, milo_summary.csv, figure)
"""

import sys
import warnings
from pathlib import Path

import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import mudata as md
import numpy as np
import pandas as pd
import pertpy as pt
import scanpy as sc

warnings.filterwarnings("ignore")
DATA, RES = Path(sys.argv[1]), Path(sys.argv[2])
OUT = RES / "sccoda_milo"
OUT.mkdir(parents=True, exist_ok=True)
ZONES = ["CTRL", "RZ", "BZ", "IZ", "FZ"]
MM = 1 / 25.4
mpl.rcParams.update({"font.family": "Arial", "font.size": 5, "axes.titlesize": 6, "pdf.fonttype": 42})

e = np.load(RES / "exp1_harmony_covariates" / "embeddings.npz", allow_pickle=True)
names = pd.Index(e["obs_names"].astype(str))
run = names.str.split("_", n=1).str[0]
sp = []
for r in pd.unique(run):
    t = pd.read_csv(DATA / "visium_tables" / f"{r}_spots.csv", usecols=["barcode", "major_labl", "patient"])
    t.index = r + "_" + t.barcode
    sp.append(t)
sp = pd.concat(sp).loc[names]


def load(folder, col):
    lab = pd.read_csv(RES / folder / "exp2_labels.csv.gz", index_col=0, dtype=str)[col].loc[names]
    order = lab.value_counts().index
    return lab.map({d: f"D{i + 1}" for i, d in enumerate(order)})


CLUST = {"RNA + space": load("exp2_joint_harmony", "r0.75_s0"), "RNA + TF + space": load("exp3_joint_rna_tf", "r0.75_s0"),
         "RNA + CCC + space": load("exp3_joint_rna_ccc_ours", "r0.25_s0")}
STRESSED = {"RNA + TF + space": "D5", "RNA + CCC + space": "D6"}

# ---------------- scCODA ----------------
rows = []
for name, lab in CLUST.items():
    a = ad.AnnData(obs=pd.DataFrame({"sample": run, "zone": pd.Categorical(sp.major_labl.to_numpy(), categories=ZONES),
                                     "domain": lab.to_numpy()}, index=names))
    keep = a.obs.domain.map(a.obs.domain.value_counts()) >= 200  # drop near-empty domains
    a = a[keep.to_numpy()].copy()
    sccoda = pt.tl.Sccoda()
    mdata = sccoda.load(a, type="cell_level", generate_sample_level=True, cell_type_identifier="domain",
                        sample_identifier="sample", covariate_obs=["zone"])
    mdata = sccoda.prepare(mdata, formula="C(zone, Treatment('CTRL'))", reference_cell_type="automatic")
    sccoda.run_nuts(mdata, rng_key=0, num_samples=10000, num_warmup=1000)
    sccoda.set_fdr(mdata, est_fdr=0.2)
    eff = sccoda.get_effect_df(mdata).reset_index()
    eff["clustering"] = name
    eff["reference"] = mdata["coda"].uns["scCODA_params"]["reference_cell_type"]
    rows.append(eff)
    print(name, "reference", eff.reference.iloc[0], flush=True)
eff = pd.concat(rows)
eff.to_csv(OUT / "sccoda_effects.csv", index=False)
cred = eff[eff["Final Parameter"] != 0]
print(cred[["clustering", "Covariate", "Cell Type", "Final Parameter", "log2-fold change"]].to_string(index=False), flush=True)

# ---------------- Milo ----------------
a = ad.AnnData(obs=pd.DataFrame({"sample": run, "zone": pd.Categorical(sp.major_labl.to_numpy(), categories=ZONES),
                                 "patient": sp.patient.to_numpy()}, index=names))
a.obsm["X_harmony"] = e["authors"].astype(np.float32)
for name, lab in CLUST.items():
    a.obs[name] = lab.to_numpy()
cons = pd.read_csv(RES / "grid_search" / "family_consensus_spots.csv.gz", index_col=0).loc[names]
fams = pd.read_csv(RES / "grid_search" / "candidate_families.csv")
top_fams = fams.sort_values(["n_configs", "n_instances"], ascending=False).family.head(12).tolist()
a.obs["region"] = "other"
for f in top_fams[::-1]:  # larger / more robust families written last win
    a.obs.loc[cons[f"family_{f}"].to_numpy(), "region"] = f"family_{f}"
a.obs.loc[(CLUST["RNA + TF + space"] == "D5").to_numpy(), "region"] = "stressed_CM"
sc.pp.neighbors(a, n_neighbors=30, use_rep="X_harmony", random_state=0)
milo = pt.tl.Milo()
mdata = milo.load(a)
milo.make_nhoods(mdata["rna"], prop=0.1, seed=0)
mdata = milo.count_nhoods(mdata, sample_col="sample")
print("neighbourhoods:", mdata["milo"].n_vars, flush=True)
res = []
for z in ZONES[1:]:
    samples = sorted(set(run[sp.major_labl.isin(["CTRL", z]).to_numpy()]))
    mdata["rna"].obs["zone_test"] = pd.Categorical(np.where(mdata["rna"].obs.zone == z, z, "CTRL"), categories=["CTRL", z])
    milo.da_nhoods(mdata, design="~zone_test", subset_samples=samples, solver="pydeseq2")
    v = mdata["milo"].var[["logFC", "PValue", "SpatialFDR"]].copy()
    v["contrast"] = f"{z} vs CTRL"
    for col in ["region", "RNA + TF + space"]:
        milo.annotate_nhoods(mdata, anno_col=col)
        v[col] = mdata["milo"].var["nhood_annotation"].astype(str).to_numpy()
        v[col + "_frac"] = mdata["milo"].var["nhood_annotation_frac"].to_numpy()
    res.append(v.reset_index(names="nhood"))
    print(z, "done", flush=True)
nh = pd.concat(res)
nh.to_csv(OUT / "milo_nhoods.csv", index=False)
summ = nh.groupby(["contrast", "region"]).agg(n_nhoods=("logFC", "size"), median_logFC=("logFC", "median"),
                                              frac_FDR_lt_0_1_up=("SpatialFDR", lambda s: 0.0)).reset_index()
up = nh.assign(up=(nh.SpatialFDR < 0.1) & (nh.logFC > 0), down=(nh.SpatialFDR < 0.1) & (nh.logFC < 0))
summ = up.groupby(["contrast", "region"]).agg(n_nhoods=("logFC", "size"), median_logFC=("logFC", "median"),
                                              share_up=("up", "mean"), share_down=("down", "mean")).reset_index()
summ.to_csv(OUT / "milo_summary.csv", index=False)
print(summ.round(2).to_string(index=False))

# ---------------- figure: Milo beeswarm per region and contrast + scCODA heatmap ----------------
regions = ["stressed_CM"] + [f"family_{f}" for f in top_fams] + ["other"]
fig = plt.figure(figsize=(183 * MM, 120 * MM))
gs = fig.add_gridspec(1, 5, width_ratios=[1, 1, 1, 1, 1.4], wspace=0.15)
rng = np.random.default_rng(0)
for j, z in enumerate(ZONES[1:]):
    ax = fig.add_subplot(gs[j])
    d = nh[nh.contrast == f"{z} vs CTRL"]
    for i, rg in enumerate(regions):
        dd = d[d.region == rg]
        if dd.empty:
            continue
        y = i + rng.uniform(-0.3, 0.3, len(dd))
        col = np.where(dd.SpatialFDR < 0.1, np.where(dd.logFC > 0, "#c0392b", "#2c7fb8"), "#cfcfcf")
        ax.scatter(dd.logFC, y, s=1.2, c=col, linewidths=0, rasterized=True)
    ax.axvline(0, color="black", lw=0.4)
    ax.set_yticks(range(len(regions)), regions if j == 0 else [""] * len(regions))
    ax.set_ylim(len(regions) - 0.5, -0.5)
    ax.set_xlabel("log fold change")
    ax.set_title(f"Milo: {z} vs control", fontweight="bold")
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)
    ax.tick_params(length=1.5, width=0.4)
ax = fig.add_subplot(gs[4])
piv = eff.assign(key=eff.clustering.str.replace(" + space", "") + " " + eff["Cell Type"].astype(str)) \
         .pivot_table(index="key", columns="Covariate", values="Final Parameter")
piv.columns = [c.split("T.")[-1].rstrip("]") for c in piv.columns]
piv = piv.reindex(columns=ZONES[1:])
lim = np.nanmax(np.abs(piv.to_numpy())) or 1
ax.imshow(piv.to_numpy(), cmap="RdBu_r", vmin=-lim, vmax=lim, aspect="auto")
ax.set_xticks(range(piv.shape[1]), piv.columns)
ax.set_yticks(range(piv.shape[0]), piv.index, fontsize=4)
ax.set_title("scCODA credible effects\n(vs control; 0 = not credible)", fontweight="bold")
ax.tick_params(length=1.5, width=0.4)
fig.savefig(OUT / "fig_sccoda_milo.png", dpi=250, bbox_inches="tight")
fig.savefig(OUT / "fig_sccoda_milo.pdf", dpi=300, bbox_inches="tight")
