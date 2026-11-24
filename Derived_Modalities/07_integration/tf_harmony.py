"""Harmony on the CollecTRI TF layer across all 27 sections, compared with the unintegrated layer.

TF matrix: ULM scores of the 622 TFs scored in all 27 sections (collectri_all27), each TF z-scored
over all spots pooled (no per-section centring, which also erases zone; see tf_unintegrated_umap.py),
PCA 30 (seed 0) computed once. Two embeddings of the same PCs:
  unintegrated   the PCs as they are
  harmony        harmonypy on the PCs with the covariates fixed by experiment 1
                 (patient + sample + visium_batch; 07_integration/exp1_harmony_covariates/RESULT.md)
Each: 15-NN, UMAP (seed 0).

Metrics (as in experiment 1), pooled and per zone:
  LISI on slide / patient / Visium batch (median effective number in a spot's neighbourhood;
  1 = all neighbours from one), iLISI batch scaled to [0, 1], cLISI zone scaled (1 = pure),
  silhouette of zone on the embedding (stratified subsample, 3,000 spots per zone, seed 0).

Outputs (results/kuppe_2022/tf_harmony/): tf_embeddings.npz (obs_names, X_pca, X_harmony, UMAPs),
tf_obs.csv.gz, tf_harmony_metrics.csv, tf_harmony_umaps.png/.pdf.

Usage (env spatialleiden, from Derived_Modalities/):
    python 07_integration/tf_harmony.py data/kuppe_2022 results/kuppe_2022/tf_harmony
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
from sklearn.metrics import silhouette_samples

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))
from plot_patient_region_umap import (PATIENT_COLOURS, REGION_COLOURS, legend,  # noqa: E402
                                      patient_order, umap_arrows)

SEED = 0
COVARIATES = ["patient", "sample", "visium_batch"]  # experiment 1 decision
ZONES = ["CTRL", "RZ", "BZ", "IZ", "FZ"]
SIL_PER_ZONE = 3000
REN = {"batch": "visium_batch", "major_labl": "zone", "hca_sample_id": "sample"}

data, out = Path(sys.argv[1]), Path(sys.argv[2])
out.mkdir(parents=True, exist_ok=True)
MOD = data / "modalities" / "collectri_all27"

runs = sorted(p.stem.replace("collectri_ulm_", "") for p in MOD.glob("collectri_ulm_*.h5ad"))
assert len(runs) == 27
parts = [ad.read_h5ad(MOD / f"collectri_ulm_{r}.h5ad") for r in runs]
common = sorted(set.intersection(*(set(p.var_names) for p in parts)))
for r, p in zip(runs, parts):
    p.obs_names = [f"{r}_{b}" for b in p.obs_names]
X = np.vstack([p[:, common].X for p in parts]).astype(np.float64)
obs = pd.concat([p.obs[["batch", "major_labl", "hca_sample_id", "patient", "celltype_niche"]] for p in parts]
                ).rename(columns=REN).astype(str)
sd = X.std(0)
X = (X - X.mean(0)) / np.where(sd > 0, sd, 1)
print(f"{X.shape[0]} spots x {X.shape[1]} TFs (pooled z-score)", flush=True)

a = ad.AnnData(X.astype(np.float32), obs=obs.copy())
sc.pp.pca(a, n_comps=30, random_state=SEED)
pcs = a.obsm["X_pca"]
ho = harmonypy.run_harmony(pcs, obs, COVARIATES, random_state=SEED, verbose=False)
z = np.asarray(ho.Z_corr)
emb = {"unintegrated": pcs, "harmony": z.T if z.shape[0] != len(obs) else z}
print("Harmony done", flush=True)

umaps = {}
for name, e in emb.items():
    b = ad.AnnData(obs=obs.copy())
    b.obsm["X_e"] = e
    sc.pp.neighbors(b, n_neighbors=15, use_rep="X_e", random_state=SEED)
    sc.tl.umap(b, random_state=SEED)
    umaps[name] = b.obsm["X_umap"]
    print(f"UMAP {name} done", flush=True)
np.savez_compressed(out / "tf_embeddings.npz", obs_names=np.asarray(obs.index), X_pca=emb["unintegrated"],
                    X_harmony=emb["harmony"], umap_unintegrated=umaps["unintegrated"], umap_harmony=umaps["harmony"])
obs.to_csv(out / "tf_obs.csv.gz")

rng = np.random.default_rng(SEED)
sil_idx = np.concatenate([rng.choice(np.where(obs.zone == zn)[0], min(SIL_PER_ZONE, int((obs.zone == zn).sum())),
                                     replace=False) for zn in ZONES])
rows = []
for name, e in emb.items():
    lisi = harmonypy.compute_lisi(e, obs[["sample", "patient", "visium_batch", "zone"]].reset_index(drop=True),
                                  ["sample", "patient", "visium_batch", "zone"])
    sil = np.full(len(obs), np.nan)
    sil[sil_idx] = silhouette_samples(e[sil_idx], obs.zone.to_numpy()[sil_idx])
    for zn in ["all"] + ZONES:
        m = np.ones(len(obs), bool) if zn == "all" else (obs.zone == zn).to_numpy()
        rows.append({"embedding": name, "zone": zn, "n_spots": int(m.sum()),
                     "LISI_slide_median": float(np.median(lisi[m, 0])),
                     "LISI_patient_median": float(np.median(lisi[m, 1])),
                     "LISI_batch_median": float(np.median(lisi[m, 2])),
                     "iLISI_batch_scaled": float(np.mean(lisi[m, 2] - 1)),
                     "cLISI_zone_scaled": float(np.mean((len(ZONES) - lisi[m, 3]) / (len(ZONES) - 1))),
                     "silhouette_zone": float(np.nanmean(sil[m])),
                     "frac_spots_LISI_slide_lt_1.5": float((lisi[m, 0] < 1.5).mean())})
    print(f"metrics {name} done", flush=True)
met = pd.DataFrame(rows)
met.to_csv(out / "tf_harmony_metrics.csv", index=False)

# ---- figure: rows = unintegrated / Harmony; columns = slide, patient, zone ------------------
order = rng.permutation(len(obs))
slides = sorted(obs["sample"].unique())
scm = list(plt.get_cmap("tab20").colors) + list(plt.get_cmap("tab20b").colors)
scol = {s: scm[i % len(scm)] for i, s in enumerate(slides)}
pats = patient_order(obs.patient.unique())
zone_lab = obs.zone.replace({"CTRL": "control"}).to_numpy()
regs = [r for r in REGION_COLOURS if r in set(zone_lab)]
fig, ax = plt.subplots(2, 3, figsize=(19, 11.5))
for i, name in enumerate(["unintegrated", "harmony"]):
    u = umaps[name][order]
    m = met[(met.embedding == name) & (met.zone == "all")].iloc[0]
    title = "TF, unintegrated" if name == "unintegrated" else "TF, Harmony (patient + sample + batch)"
    for j, (lab, cols, keys, what, stat) in enumerate([
            (obs["sample"].to_numpy(), scol, slides, "slide", f"LISI slide {m.LISI_slide_median:.2f}"),
            (obs.patient.to_numpy(), PATIENT_COLOURS, pats, "patient", f"LISI patient {m.LISI_patient_median:.2f}"),
            (zone_lab, REGION_COLOURS, regs, "region", f"silhouette zone {m.silhouette_zone:.3f}")]):
        ax[i, j].scatter(u[:, 0], u[:, 1], c=[cols[k] for k in lab[order]], s=0.6, linewidths=0, rasterized=True)
        umap_arrows(ax[i, j], f"{title}\n{what} ({stat})")
        legend(ax[i, j], cols, keys, ncol=2 if len(keys) > 22 else 1)
        if what == "slide":
            ax[i, j].get_legend().remove()
fig.tight_layout()
fig.savefig(out / "tf_harmony_umaps.png", dpi=130, bbox_inches="tight")
fig.savefig(out / "tf_harmony_umaps.pdf", bbox_inches="tight")
print(met[met.zone == "all"].round(3).to_string(index=False))
print(met[met.zone != "all"].pivot(index="zone", columns="embedding",
                                   values=["LISI_slide_median", "silhouette_zone"]).round(3).to_string())
