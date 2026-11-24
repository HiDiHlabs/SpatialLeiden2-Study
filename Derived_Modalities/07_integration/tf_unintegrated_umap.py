"""Unintegrated UMAPs of the CollecTRI TF layer across all 27 sections: is there slide / patient
structure that would need Harmony?

TF matrix: ULM scores of the 622 TFs scored in all 27 sections (collectri_all27). Two scalings,
both without any batch correction:
  pooled     each TF z-scored over all 88,704 spots (fully uncorrected)
  per-slide  each TF z-scored within its section (= our TF layer recipe, preprocess_tf); this
             removes per-section mean and scale, a mild correction in itself
Then PCA 30, 15-NN, UMAP (seed 0). Reference row: the uncorrected RNA embedding of experiment 1
(variant `none`).

Mixing metric: LISI (harmonypy.compute_lisi, perplexity 30) on slide, patient and Visium batch =
effective number of slides / patients / batches in each spot's neighbourhood (1 = all neighbours
from one; higher = mixed). Reported as the median over spots, with the RNA `none` and `authors`
(Harmony) embeddings of experiment 1 for scale.

Usage (env spatialleiden, from Derived_Modalities/):
    python 07_integration/tf_unintegrated_umap.py data/kuppe_2022 results/kuppe_2022/tf_unintegrated
"""
import sys
from pathlib import Path

import anndata as ad
import harmonypy
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

SEED = 0
data, out = Path(sys.argv[1]), Path(sys.argv[2])
out.mkdir(parents=True, exist_ok=True)
MOD = data / "modalities" / "collectri_all27"
EXP1 = data.parent.parent / "results" / data.name / "exp1_harmony_covariates"
REN = {"batch": "visium_batch", "major_labl": "zone", "hca_sample_id": "sample"}

runs = sorted(p.stem.replace("collectri_ulm_", "") for p in MOD.glob("collectri_ulm_*.h5ad"))
assert len(runs) == 27
parts = [ad.read_h5ad(MOD / f"collectri_ulm_{r}.h5ad") for r in runs]
common = sorted(set.intersection(*(set(p.var_names) for p in parts)))
X = []
for r, p in zip(runs, parts):
    p.obs_names = [f"{r}_{b}" for b in p.obs_names]
    X.append(p[:, common].X)
obs = pd.concat([p.obs[["batch", "major_labl", "hca_sample_id", "patient"]] for p in parts]).rename(columns=REN).astype(str)
X = np.vstack(X).astype(np.float64)
sid = obs["sample"].to_numpy()
print(f"{X.shape[0]} spots x {X.shape[1]} TFs", flush=True)


def zs(m):
    sd = m.std(0)
    sd[sd == 0] = 1
    return (m - m.mean(0)) / sd


versions = {"TF, pooled z-score": zs(X)}
Xs = np.empty_like(X)
for s in np.unique(sid):
    Xs[sid == s] = zs(X[sid == s])
versions["TF, per-slide z-score"] = Xs

emb, umaps = {}, {}
for name, m in versions.items():
    a = ad.AnnData(m.astype(np.float32), obs=obs.copy())
    sc.pp.pca(a, n_comps=30, random_state=SEED)
    sc.pp.neighbors(a, n_neighbors=15, use_rep="X_pca", random_state=SEED)
    sc.tl.umap(a, random_state=SEED)
    emb[name], umaps[name] = a.obsm["X_pca"], a.obsm["X_umap"]
    print(f"{name}: UMAP done", flush=True)

# RNA reference (experiment 1), aligned to the same spots
e1 = np.load(EXP1 / "embeddings.npz", allow_pickle=True)
u1 = np.load(EXP1 / "umaps.npz")
pos = pd.Series(np.arange(len(e1["obs_names"])), index=e1["obs_names"]).loc[obs.index].to_numpy()
umaps["RNA, uncorrected (exp. 1 none)"] = u1["none"][pos]
emb["RNA, uncorrected (exp. 1 none)"] = e1["none"][pos]
emb["RNA, Harmony authors (exp. 1)"] = e1["authors"][pos]

rows = []
for name, z in emb.items():
    lisi = harmonypy.compute_lisi(z, obs[["sample", "patient", "visium_batch"]].reset_index(drop=True),
                                  ["sample", "patient", "visium_batch"])
    rows.append({"embedding": name, "LISI_slide_median": np.median(lisi[:, 0]),
                 "LISI_patient_median": np.median(lisi[:, 1]), "LISI_batch_median": np.median(lisi[:, 2]),
                 "frac_spots_LISI_slide_lt_1.5": float((lisi[:, 0] < 1.5).mean())})
    print(rows[-1], flush=True)
mix = pd.DataFrame(rows)
mix.to_csv(out / "tf_unintegrated_mixing.csv", index=False)

order = np.random.default_rng(SEED).permutation(len(obs))
cols = [("sample", "slide (27)"), ("patient", "patient (20)"), ("zone", "zone")]
cm = list(plt.get_cmap("tab20").colors) + list(plt.get_cmap("tab20b").colors)
zcol = {"CTRL": "#4CAF50", "FZ": "#009db2", "IZ": "#b50084", "RZ": "#ffc107", "BZ": "#f57f17"}
fig, ax = plt.subplots(len(umaps), 3, figsize=(17, 5.2 * len(umaps)))
for i, (name, u) in enumerate(umaps.items()):
    lis = mix.set_index("embedding").loc[name]
    for j, (c, t) in enumerate(cols):
        cats = sorted(obs[c].unique())
        col = zcol if c == "zone" else {k: cm[q % len(cm)] for q, k in enumerate(cats)}
        ax[i, j].scatter(*u[order].T, c=[col[k] for k in obs[c].to_numpy()[order]], s=0.2, rasterized=True)
        ax[i, j].set(xticks=[], yticks=[])
        sub = {"sample": f"median LISI slide {lis.LISI_slide_median:.2f}",
               "patient": f"median LISI patient {lis.LISI_patient_median:.2f}",
               "zone": f"median LISI batch {lis.LISI_batch_median:.2f}"}[c]
        ax[i, j].set_title(f"{name}\ncoloured by {t} ({sub})", fontsize=10)
        if c != "sample":
            ax[i, j].legend(handles=[plt.Line2D([], [], ls="", marker="o", color=col[k], label=k) for k in
                                     (sorted(cats, key=lambda p: int(p[1:])) if c == "patient" else cats)],
                            fontsize=6, ncol=2 if c == "patient" else 1, loc="best", markerscale=1.5)
fig.suptitle("Unintegrated embeddings: does slide / patient structure need Harmony?", fontsize=13)
fig.tight_layout()
fig.savefig(out / "tf_unintegrated_umaps.png", dpi=110)
print(mix.round(3).to_string(index=False))
