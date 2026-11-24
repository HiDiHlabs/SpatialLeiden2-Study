"""CCC layer embedding across all 27 sections for the joint SpatialLeiden runs, built like the
TF layer (07_integration/tf_harmony.py, RESULT_tf_harmony.md).

Input: the colleague's CCC ligand-receptor local scores for all 27 slides, as either
  * a folder of per-slide tables ccc_<run>.csv (index = spot barcode, columns
    "LIGAND^RECEPTOR"), or
  * one .h5ad (spots x pairs) with an obs column naming the Visium run (--sample-col,
    default "hca_sample_id") and obs_names = barcodes or "<run>_<barcode>".
<run> is the Visium run id used everywhere else (10X001, ACH0012, ...; see
data/kuppe_2022/visium_raw/metadata-Visium.csv, column hca_sample_id).

Recipe:
  * pairs present (non-missing) on all 27 slides;
  * each pair z-scored over all spots pooled; PCA 30 (seed 0);
  * mode "harmony" (for the NON-integrated scores, recommended): harmonypy on the PCs with the
    covariates of experiment 1 (patient + sample + visium_batch), seed 0, as for RNA and TF;
  * mode "pca" (for scores the colleague already integrated): PCA only, no second correction.

Outputs (in OUT): ccc_embeddings.npz with obs_names "<run>_<barcode>" and X_harmony (mode
harmony) or X_pca (both modes); ccc_pairs.txt (pairs used).

Usage (env spatialleiden, from Derived_Modalities/):
    python 07_integration/ccc_harmony.py harmony data/kuppe_2022 <ccc_input> results/kuppe_2022/ccc_harmony
    python 07_integration/ccc_harmony.py pca     data/kuppe_2022 <ccc_integrated_input> results/kuppe_2022/ccc_integrated
"""
import argparse
from pathlib import Path

import anndata as ad
import harmonypy
import numpy as np
import pandas as pd
import scanpy as sc

SEED = 0
COVARIATES = ["patient", "sample", "visium_batch"]  # experiment 1 decision
MIN_PAIRS = 50

p = argparse.ArgumentParser()
p.add_argument("mode", choices=["harmony", "pca"])
p.add_argument("data", type=Path)
p.add_argument("ccc", type=Path, help="folder of ccc_<run>.csv, or one .h5ad")
p.add_argument("out", type=Path)
p.add_argument("--sample-col", default="hca_sample_id")
args = p.parse_args()
args.out.mkdir(parents=True, exist_ok=True)

# ---- read: dict run -> DataFrame (barcode x pair) ----
if args.ccc.is_dir():
    parts = {f.stem.removeprefix("ccc_"): pd.read_csv(f, index_col=0)
             for f in sorted(args.ccc.glob("ccc_*.csv"))}
else:
    a = ad.read_h5ad(args.ccc)
    X = a.X.toarray() if hasattr(a.X, "toarray") else np.asarray(a.X)
    df = pd.DataFrame(X, index=a.obs_names.astype(str), columns=a.var_names.astype(str))
    smp = a.obs[args.sample_col].astype(str)
    parts = {}
    for r in smp.unique():
        x = df[smp.to_numpy() == r]
        x.index = [i.removeprefix(f"{r}_") for i in x.index]
        parts[r] = x
runs = sorted(parts)
assert len(runs) == 27, f"{len(runs)} slides: {runs}"

common = sorted(set.intersection(*(set(x.columns[x.notna().all()]) for x in parts.values())))
print(f"pairs per slide {min(x.shape[1] for x in parts.values())}-"
      f"{max(x.shape[1] for x in parts.values())}; on all 27 slides: {len(common)}", flush=True)
assert len(common) >= MIN_PAIRS, "too few pairs shared by all slides"

X, obs = [], []
for r in runs:
    spots = pd.read_csv(args.data / "visium_tables" / f"{r}_spots.csv").set_index("barcode")
    x = parts[r]
    missing = set(spots.index) - set(x.index)
    assert not missing, f"{r}: {len(missing)} spots of the spot table have no CCC scores"
    X.append(x.loc[spots.index, common].to_numpy(np.float64))
    o = spots[["batch", "major_labl", "hca_sample_id", "patient"]].astype(str)
    o.index = [f"{r}_{b}" for b in spots.index]
    obs.append(o)
X = np.vstack(X)
obs = pd.concat(obs).rename(columns={"batch": "visium_batch", "major_labl": "zone",
                                     "hca_sample_id": "sample"})
sd = X.std(0)
X = (X - X.mean(0)) / np.where(sd > 0, sd, 1)
print(f"{X.shape[0]} spots x {X.shape[1]} pairs (pooled z-score)", flush=True)

a = ad.AnnData(X.astype(np.float32), obs=obs.copy())
sc.pp.pca(a, n_comps=30, random_state=SEED)
emb = {"X_pca": a.obsm["X_pca"]}
if args.mode == "harmony":
    ho = harmonypy.run_harmony(emb["X_pca"], obs, COVARIATES, random_state=SEED, verbose=False)
    z = np.asarray(ho.Z_corr)
    emb["X_harmony"] = z.T if z.shape[0] != len(obs) else z
np.savez_compressed(args.out / "ccc_embeddings.npz", obs_names=np.asarray(obs.index), **emb)
(args.out / "ccc_pairs.txt").write_text("\n".join(common) + "\n")
print(f"done: {list(emb)}", flush=True)
