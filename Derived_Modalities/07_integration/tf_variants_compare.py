"""Three TF-layer variants on both axes (batch mixing vs zone preservation), same decision rule as
experiment 1. Adds the per-slide z-score variant to the two of tf_harmony.py, with identical metrics
(the none / harmony rows are read from tf_harmony_metrics.csv; same definitions and subsample).

  none        each TF z-scored over all spots pooled, PCA 30            (tf_harmony.py, X_pca)
  harmony     the same PCs, Harmony on patient + sample + visium_batch  (tf_harmony.py, X_harmony)
  per_slide   each TF z-scored within its section, PCA 30 (computed here, seed 0)

Metrics as in tf_harmony.py: LISI on slide / patient / batch (median), scaled iLISI batch, scaled
cLISI zone, silhouette of zone on the same stratified subsample (3,000 spots per zone, seed 0);
pooled and per zone. Writes tf_variants_compare.csv (all rows) and tf_variants_summary.csv (pooled).

Usage (env spatialleiden, from Derived_Modalities/, after tf_harmony.py):
    python 07_integration/tf_variants_compare.py data/kuppe_2022 results/kuppe_2022/tf_harmony
"""
import sys
from pathlib import Path

import anndata as ad
import harmonypy
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.metrics import silhouette_samples

SEED = 0
ZONES = ["CTRL", "RZ", "BZ", "IZ", "FZ"]
SIL_PER_ZONE = 3000
data, out = Path(sys.argv[1]), Path(sys.argv[2])
MOD = data / "modalities" / "collectri_all27"

e = np.load(out / "tf_embeddings.npz", allow_pickle=True)
obs = pd.read_csv(out / "tf_obs.csv.gz", index_col=0).astype(str)
assert list(obs.index) == list(e["obs_names"])

# per-slide z-score of the same 622 TFs, same spot order
runs = sorted(p.stem.replace("collectri_ulm_", "") for p in MOD.glob("collectri_ulm_*.h5ad"))
parts = [ad.read_h5ad(MOD / f"collectri_ulm_{r}.h5ad") for r in runs]
common = sorted(set.intersection(*(set(p.var_names) for p in parts)))
blocks = []
for r, p in zip(runs, parts):
    m = p[:, common].X.astype(np.float64)
    sd = m.std(0)
    blocks.append(pd.DataFrame((m - m.mean(0)) / np.where(sd > 0, sd, 1), index=[f"{r}_{b}" for b in p.obs_names]))
Xs = pd.concat(blocks).loc[obs.index].to_numpy(np.float32)
b = ad.AnnData(Xs)
sc.pp.pca(b, n_comps=30, random_state=SEED)

emb = {"per_slide": b.obsm["X_pca"]}  # none / harmony: metrics already in tf_harmony_metrics.csv
rng = np.random.default_rng(SEED)  # same draws as tf_harmony.py
sil_idx = np.concatenate([rng.choice(np.where(obs.zone == zn)[0], min(SIL_PER_ZONE, int((obs.zone == zn).sum())),
                                     replace=False) for zn in ZONES])
rows = []
for name, z in emb.items():
    lisi = harmonypy.compute_lisi(z, obs[["sample", "patient", "visium_batch", "zone"]].reset_index(drop=True),
                                  ["sample", "patient", "visium_batch", "zone"])
    sil = np.full(len(obs), np.nan)
    sil[sil_idx] = silhouette_samples(z[sil_idx], obs.zone.to_numpy()[sil_idx])
    for zn in ["all"] + ZONES:
        m = np.ones(len(obs), bool) if zn == "all" else (obs.zone == zn).to_numpy()
        rows.append({"variant": name, "zone": zn, "n_spots": int(m.sum()),
                     "LISI_slide_median": float(np.median(lisi[m, 0])),
                     "LISI_patient_median": float(np.median(lisi[m, 1])),
                     "LISI_batch_median": float(np.median(lisi[m, 2])),
                     "iLISI_batch_scaled": float(np.mean(lisi[m, 2] - 1)),
                     "cLISI_zone_scaled": float(np.mean((len(ZONES) - lisi[m, 3]) / (len(ZONES) - 1))),
                     "silhouette_zone": float(np.nanmean(sil[m])),
                     "frac_spots_LISI_slide_lt_1.5": float((lisi[m, 0] < 1.5).mean())})
    print(f"{name} done", flush=True)
prev = pd.read_csv(out / "tf_harmony_metrics.csv").rename(columns={"embedding": "variant"})
prev["variant"] = prev["variant"].replace({"unintegrated": "none"})
df = pd.concat([prev, pd.DataFrame(rows)], ignore_index=True)
df.to_csv(out / "tf_variants_compare.csv", index=False)
df[df.zone == "all"].to_csv(out / "tf_variants_summary.csv", index=False)
print(df[df.zone == "all"].drop(columns=["zone"]).round(3).to_string(index=False))
for col in ["silhouette_zone", "cLISI_zone_scaled", "LISI_slide_median"]:
    print(f"\n{col} by zone")
    print(df[df.zone != "all"].pivot(index="variant", columns="zone", values=col)[ZONES].loc[["none", "harmony", "per_slide"]].round(3).to_string())
