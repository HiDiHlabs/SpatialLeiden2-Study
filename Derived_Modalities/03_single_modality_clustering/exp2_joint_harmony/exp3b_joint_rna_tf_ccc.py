"""Joint multimodal SpatialLeiden over all 27 slides with THREE expression-type layers: RNA + TF + CCC
(+ the spatial graph). Same recipe as exp3_joint_multimodal_sweep.py (which takes one extra layer):
  layers  rna  15-NN on the RNA Harmony embedding (exp1_harmony_covariates/embeddings.npz, "authors")
          tf   15-NN on the TF Harmony embedding (tf_harmony/tf_embeddings.npz, "X_harmony")
          ccc  15-NN on our LIANA+ CCC Harmony embedding (ccc_harmony_ours/ccc_embeddings.npz, "X_harmony")
          spatial  squidpy grid, 6 neighbours, per slide
  weights rna 1, tf 0.5, ccc 0.5 (as rna_all in mm_configs.py), spatial = the swept weight
  K = 9; seed(s) from EXP_SEEDS (default 0); spatial weights from EXP_RATIOS (default 0,0.75).
  Holding K: the shared expression resolution is searched at spatial weight 0 with the multimodal
  call itself; at weight > 0 the spatial resolution is searched (start 0.4); if that fails, the
  expression resolution is searched with the spatial resolution fixed (fallback, recorded).
Outputs (OUT): exp2_labels.csv.gz (columns r<weight>_s<seed>, same format as the other joint runs),
exp3b_sweep.csv, exp3b.log.

Usage (from Derived_Modalities/, env spatialleiden_joint):
    python 03_single_modality_clustering/exp2_joint_harmony/exp3b_joint_rna_tf_ccc.py data/kuppe_2022 results/kuppe_2022 \
        results/kuppe_2022/exp3b_joint_rna_tf_ccc [workers]
"""

import os

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys  # noqa: E402
import time  # noqa: E402
from multiprocessing import get_context  # noqa: E402
from pathlib import Path  # noqa: E402

import anndata as ad  # noqa: E402
import mudata as md  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scanpy as sc  # noqa: E402
import spatialleiden as sl  # noqa: E402
import squidpy as sq  # noqa: E402

DATA, RES, OUT = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
N_WORKERS = int(sys.argv[4]) if len(sys.argv) > 4 else 4
OUT.mkdir(parents=True, exist_ok=True)
K = 9
SEEDS = [int(x) for x in os.environ.get("EXP_SEEDS", "0").split(",")]
RATIOS = [float(x) if float(x) else 0 for x in os.environ.get("EXP_RATIOS", "0,0.75").split(",")]
W = {"rna": 1.0, "tf": 0.5, "ccc": 0.5}
LAYERS = {"rna": ("exp1_harmony_covariates/embeddings.npz", "authors"),
          "tf": ("tf_harmony/tf_embeddings.npz", "X_harmony"),
          "ccc": ("ccc_harmony_ours/ccc_embeddings.npz", "X_harmony")}
LOG = open(OUT / "exp3b.log", "a")


def log(msg):
    print(msg, flush=True)
    LOG.write(msg + "\n")
    LOG.flush()


def search(n_at, target, start=1.0, step=0.1, n_iter=30):
    """Step-halving search, identical to experiments 0-3. Returns (resolution, reached)."""
    res, n, i = start, n_at(start), 1
    while n != target:
        if i >= n_iter:
            return res, False
        sign = 1 if n < target else -1
        res += step * sign
        n = n_at(res)
        if (1 if n < target else -1) != sign:
            step /= 2
        i += 1
    return res, True


t0 = time.time()
emb = {}
names = None
for k, (f, key) in LAYERS.items():
    e = np.load(RES / f, allow_pickle=True)
    idx = pd.Index(e["obs_names"].astype(str))
    names = idx if names is None else names
    emb[k] = pd.DataFrame(e[key], index=idx).loc[names].to_numpy(np.float32)
sample = names.str.split("_", n=1).str[0]
frames = []
for s in pd.unique(sample):
    t = pd.read_csv(DATA / "visium_tables" / f"{s}_spots.csv", usecols=["barcode", "pxl_row_fullres", "pxl_col_fullres"])
    t.index = s + "_" + t.barcode
    frames.append(t)
spots = pd.concat(frames).loc[names]
mods = {}
for k, X in emb.items():
    a = ad.AnnData(obs=pd.DataFrame(index=names))
    a.obsm["X_e"] = X
    sc.pp.neighbors(a, n_neighbors=15, use_rep="X_e", random_state=0)
    mods[k] = ad.AnnData(obs=pd.DataFrame(index=names), obsp={"connectivities": a.obsp["connectivities"]})
g = ad.AnnData(obs=pd.DataFrame({"sample": pd.Categorical(sample)}, index=names))
g.obsm["spatial"] = spots[["pxl_col_fullres", "pxl_row_fullres"]].to_numpy(float)
sq.gr.spatial_neighbors(g, library_key="sample", coord_type="grid", n_neighs=6)
M = md.MuData(mods)
M.obsp["spatial_connectivities"] = g.obsp["spatial_connectivities"]
log(f"layers {list(W)} weights {W}; seeds {SEEDS}; spatial {RATIOS}; {len(names)} spots; {time.time() - t0:.0f}s")


def cluster(res_expr, res_sp, w_sp, seed):
    sl.spatialleiden_multimodal(
        M, resolution={**{k: res_expr for k in W}, "spatial": res_sp},
        layer_weights={**W, "spatial": float(w_sp)},
        directed={**{k: False for k in W}, "spatial": False},
        n_iterations=-1, random_state=seed, key_added="tmp")
    return M.obs["tmp"].astype(str).to_numpy()


def job(args):
    seed, w, r_expr = args
    t = time.time()
    if w == 0:
        r_e, ok = search(lambda r: len(set(cluster(r, 1.0, 0, seed))), K)
        return seed, w, r_e, np.nan, not ok, False, cluster(r_e, 1.0, 0, seed), time.time() - t
    r_sp, ok = search(lambda r: len(set(cluster(r_expr, r, w, seed))), K, start=0.4)
    r_e, fb = r_expr, False
    if not ok:
        r_e, ok = search(lambda r: len(set(cluster(r, r_sp, w, seed))), K, start=r_expr)
        fb = True
    return seed, w, r_e, r_sp, not ok, fb, cluster(r_e, r_sp, w, seed), time.time() - t


results = []
base = {}
for seed in SEEDS:
    r = job((seed, 0, None))
    results.append(r)
    base[seed] = r[2]
    log(f"spatial 0 seed {seed}: res_expr {r[2]:.4f} failed={r[4]} K={len(set(r[6]))} {r[7]:.0f}s")
rest = [(s, w, base[s]) for s in SEEDS for w in RATIOS if w != 0]
if rest:
    with get_context("fork").Pool(min(N_WORKERS, len(rest))) as pool:
        for r in pool.imap_unordered(job, rest):
            results.append(r)
            log(f"spatial {r[1]} seed {r[0]}: res=({r[2]:.4f}, {r[3]:.4f}) failed={r[4]} fallback={r[5]} "
                f"K={len(set(r[6]))} {r[7]:.0f}s")
labels = {f"r{w}_s{s}": lab for s, w, *_x, lab, _t in [(r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7]) for r in results]}
cols = [f"r{w}_s{s}" for w in RATIOS for s in SEEDS if f"r{w}_s{s}" in labels]
pd.DataFrame({c: labels[c] for c in cols}, index=names).to_csv(OUT / "exp2_labels.csv.gz")
pd.DataFrame([{"seed": r[0], "spatial_weight": r[1], "res_expression": r[2], "res_spatial": r[3], "search_failed": r[4],
               "expression_fallback": r[5], "n_clusters": len(set(r[6])), "seconds": round(r[7])} for r in results]
             ).to_csv(OUT / "exp3b_sweep.csv", index=False)
log(f"done in {time.time() - t0:.0f}s")
