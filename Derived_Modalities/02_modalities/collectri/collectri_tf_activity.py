"""TF activity per Visium spot with decoupler (ULM) and the CollecTRI network.

Slides: P1 control (10X001), P9 IZ rep2 (ACH0013) and P9 IZ (ACH0012), the slides
of experiment 0. Expression: raw counts from the Zenodo h5ad, joined by barcode to the
authors' per-spot table, then normalize_total (10,000) + log1p on all genes (no
scaling, no gene filtering beyond what the file holds).

Method: decoupler 2.x, dc.mt.ulm (univariate linear model) with tmin = 5 targets,
CollecTRI (human) from dc.op.collectri. ULM rather than the paper's wmean: the paper
never computed spot-level TF activity from CollecTRI, so there is nothing to match, and
ULM is the method decoupler recommends for TF activity.

Depth check (the two slides differ ~2x in genes per spot): per slide we record how many
TFs pass tmin, how many targets each has in the gene universe, and, per spot, how many
of a TF's targets are detected (count > 0). If P1 detects far fewer targets per spot,
its TF layer is noisier, and "no structure" cannot be told apart from "no signal".

Outputs
  data/kuppe_2022/modalities/collectri_ulm_<run>.h5ad   spots x TFs (X = ULM score),
      layers["padj"], obs = spot metadata in h5ad barcode order, obsm["spatial"]
  results/kuppe_2022/collectri/collectri_qc.csv          per-slide summary
  results/kuppe_2022/collectri/collectri_tf_table_<run>.csv  per-TF targets, detection, score spread
  results/kuppe_2022/collectri/figures/                   maps of selected TFs

Usage (from Derived_Modalities/; env with decoupler >= 2.1, e.g. masld_atlas):
    python 02_modalities/collectri/collectri_tf_activity.py data/kuppe_2022 results/kuppe_2022/collectri
"""

import sys
from pathlib import Path

import anndata as ad
import decoupler as dc
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib.collections import EllipseCollection
from scipy import sparse
from scipy.spatial import cKDTree

DATA, RES = Path(sys.argv[1]), Path(sys.argv[2])
MOD = DATA / "modalities"
FIG = RES / "figures"
for d in (MOD, RES, FIG):
    d.mkdir(parents=True, exist_ok=True)

SLIDES = {
    "10X001": ("Visium_control_P1.h5ad", "P1 control (10X001)"),
    "ACH0013": ("Visium_GT_IZ_P9_rep2.h5ad", "P9 IZ rep2 (ACH0013)"),
    "ACH0012": ("Visium_GT_IZ_P9.h5ad", "P9 IZ (ACH0012)"),
}
TMIN = 5
# TFs shown in the maps: inflammation, fibrosis, hypoxia, myeloid, cardiomyocyte identity
SHOW = ["NFKB1", "STAT3", "SMAD3", "HIF1A", "SPI1", "MEF2C", "GATA4", "TEAD1"]

print("decoupler", dc.__version__)
net = dc.op.collectri(organism="human")
net.to_csv(RES / "collectri_network_used.csv.gz", index=False)
print(f"CollecTRI: {net.source.nunique()} TFs, {len(net)} edges")

qc_rows, score_store = [], {}
for run, (h5_name, label) in SLIDES.items():
    h5 = ad.read_h5ad(DATA / "processed1" / h5_name)
    spots = pd.read_csv(DATA / "visium_tables" / f"{run}_spots.csv").set_index(
        "barcode"
    )
    assert set(spots.index) == set(h5.obs_names), f"{run}: barcodes differ"
    spots = spots.loc[h5.obs_names]
    a = ad.AnnData(
        X=sparse.csr_matrix(h5.raw.X),
        obs=spots.copy(),
        var=pd.DataFrame(index=h5.raw.var_names),
    )
    a.obsm["spatial"] = spots[["pxl_col_fullres", "pxl_row_fullres"]].to_numpy(float)
    counts = a.X.copy()
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)

    # --- depth / target-coverage diagnostics
    genes = set(a.var_names)
    net_in = net[net.target.isin(genes)]
    n_targets = net_in.groupby("source").size()
    kept = n_targets[n_targets >= TMIN].index
    detected = (counts > 0).astype(np.float32).tocsc()
    gidx = {g: i for i, g in enumerate(a.var_names)}
    per_tf = []
    for tf in kept:
        cols = [gidx[g] for g in net_in.loc[net_in.source == tf, "target"]]
        det = np.asarray(detected[:, cols].sum(1)).ravel()  # targets detected per spot
        per_tf.append(
            {
                "tf": tf,
                "n_targets_in_slide": len(cols),
                "median_targets_detected_per_spot": float(np.median(det)),
                "frac_spots_ge5_targets_detected": float((det >= TMIN).mean()),
            }
        )
    per_tf = pd.DataFrame(per_tf).set_index("tf")

    # --- ULM
    dc.mt.ulm(data=a, net=net, tmin=TMIN)
    scores = dc.pp.get_obsm(a, key="score_ulm")
    padj = a.obsm["padj_ulm"]
    out = ad.AnnData(
        X=scores.X.astype(np.float32),
        obs=a.obs.copy(),
        var=pd.DataFrame(index=scores.var_names),
    )
    out.layers["padj"] = np.asarray(
        padj.loc[out.obs_names, out.var_names], dtype=np.float32
    )
    out.obsm["spatial"] = a.obsm["spatial"]
    out.uns["collectri_ulm"] = {
        "decoupler_version": dc.__version__,
        "method": "ulm",
        "tmin": TMIN,
        "network": "CollecTRI human (dc.op.collectri)",
        "normalization": "normalize_total 1e4 + log1p, all genes in the file",
        "n_genes_in_slide": int(a.n_vars),
    }
    out.write_h5ad(MOD / f"collectri_ulm_{run}.h5ad", compression="gzip")
    score_store[run] = (out, label)

    sd = pd.Series(out.X.std(0), index=out.var_names)
    per_tf["score_sd_across_spots"] = sd.reindex(per_tf.index)
    per_tf["frac_spots_padj_lt_0.05"] = pd.Series(
        (out.layers["padj"] < 0.05).mean(0), index=out.var_names
    ).reindex(per_tf.index)
    per_tf.to_csv(RES / f"collectri_tf_table_{run}.csv")
    qc_rows.append(
        {
            "slide": run,
            "label": label,
            "spots": a.n_obs,
            "genes_in_file": a.n_vars,
            "median_genes_per_spot": float(
                np.median(np.asarray((counts > 0).sum(1)).ravel())
            ),
            "collectri_tfs_total": net.source.nunique(),
            "tfs_with_ge1_target_in_slide": int((n_targets >= 1).sum()),
            "tfs_passing_tmin": int(len(kept)),
            "tfs_scored": int(out.n_vars),
            "median_targets_per_tf": float(per_tf.n_targets_in_slide.median()),
            "median_targets_detected_per_spot_per_tf": float(
                per_tf.median_targets_detected_per_spot.median()
            ),
            "median_frac_spots_ge5_targets_detected": float(
                per_tf.frac_spots_ge5_targets_detected.median()
            ),
            "frac_scores_padj_lt_0.05": float((out.layers["padj"] < 0.05).mean()),
        }
    )
    print(qc_rows[-1])

qc = pd.DataFrame(qc_rows)
common = set.intersection(*(set(o.var_names) for o, _ in score_store.values()))
qc["tfs_scored_in_all_three"] = len(common)
qc.to_csv(RES / "collectri_qc.csv", index=False)

# --- maps of selected TFs, one shared symmetric colour scale per TF across slides
mpl.rcParams.update(
    {"font.family": "Arial", "font.size": 6, "axes.titlesize": 6.5, "pdf.fonttype": 42}
)
show = [t for t in SHOW if t in common]
fig, axes = plt.subplots(
    len(score_store),
    len(show),
    figsize=(183 / 25.4, 1.7 * len(score_store)),
    gridspec_kw={"wspace": 0.05, "hspace": 0.15},
)
for j, tf in enumerate(show):
    vmax = max(np.quantile(np.abs(o[:, tf].X), 0.99) for o, _ in score_store.values())
    norm = mpl.colors.TwoSlopeNorm(0, -vmax, vmax)
    for i, (run, (o, label)) in enumerate(score_store.items()):
        ax = axes[i, j]
        x, y = o.obsm["spatial"].T
        pitch = np.median(cKDTree(np.c_[x, y]).query(np.c_[x, y], k=2)[0][:, 1])
        d = np.full(len(x), 1.05 * pitch)
        coll = EllipseCollection(
            d,
            d,
            np.zeros(len(x)),
            units="xy",
            offsets=np.c_[x, y],
            offset_transform=ax.transData,
            facecolors=mpl.colormaps["RdBu_r"](norm(o[:, tf].X.ravel())),
            linewidths=0,
        )
        coll.set_rasterized(True)
        ax.add_collection(coll)
        ax.set_xlim(x.min() - pitch, x.max() + pitch)
        ax.set_ylim(y.max() + pitch, y.min() - pitch)
        ax.set_aspect("equal")
        ax.axis("off")
        if i == 0:
            ax.set_title(tf, fontweight="bold")
        if j == 0:
            ax.text(
                -0.08,
                0.5,
                label,
                transform=ax.transAxes,
                rotation=90,
                ha="right",
                va="center",
            )
fig.suptitle(
    "TF activity (decoupler ULM, CollecTRI). Red = active, blue = inactive; "
    "one colour scale per TF, shared across slides",
    fontsize=6.5,
    y=0.995,
)
fig.savefig(FIG / "collectri_selected_tfs.png", dpi=250, bbox_inches="tight")
fig.savefig(FIG / "collectri_selected_tfs.pdf", dpi=300, bbox_inches="tight")
print("done;", len(common), "TFs scored on all three slides")
