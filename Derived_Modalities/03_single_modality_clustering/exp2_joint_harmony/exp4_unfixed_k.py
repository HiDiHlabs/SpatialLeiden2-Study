"""Experiment 4: joint SpatialLeiden (all 27 slides) WITHOUT holding K fixed across spatial
weights. Complements experiments 2/3, where K = 9 is searched at every weight.

Per layer: the expression resolution is the one giving K = 9 at spatial weight 0 (taken from
experiment 2's seed-0 search, or searched here with spatialleiden(layer_ratio=0) if not given);
the spatial resolution is fixed at 1.0 (package default). Then ONE run per spatial weight, so K
is free to change with the weight. Seed 0, directed (False, False), n_iterations -1, graphs as
in experiment 2 (15-NN on the embedding; squidpy grid per slide, 6 neighbours).

Embedding sources:
  npz:<path>:<key>      obs_names "<run>_<barcode>" (exp1 RNA, tf_harmony)
  h5ad:<path>:<obsm key>  read with h5py (works for files written by newer anndata). Spots are
                        mapped to Visium runs by identical barcode sets per obs "sample";
                        barcodes come from obs "original_spot" (the colleague's CCC file).

Outputs (OUT): exp2_labels.csv.gz (columns r<weight>_s0; same format as experiment 2, so the
figure scripts read it), exp4_k_per_weight.csv (K overall and per slide), exp4.log.

Usage (from Derived_Modalities/, spatialleiden env):
    python 03_single_modality_clustering/exp2_joint_harmony/exp4_unfixed_k.py data/kuppe_2022 \
        npz:results/kuppe_2022/exp1_harmony_covariates/embeddings.npz:authors \
        results/kuppe_2022/exp4_unfixed_k_rna 0.5
    (last argument: expression resolution; "search" to search K = 9 at weight 0)
"""
import os

os.environ.setdefault("OMP_NUM_THREADS", "2")
import sys  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import anndata as ad  # noqa: E402
import h5py  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scanpy as sc  # noqa: E402
import spatialleiden as sl  # noqa: E402
import squidpy as sq  # noqa: E402

DATA, SRC, OUT = Path(sys.argv[1]), sys.argv[2], Path(sys.argv[3])
RES = sys.argv[4] if len(sys.argv) > 4 else "search"
WEIGHTS = [float(w) if float(w) else 0 for w in os.environ.get("EXP_RATIOS", "0,0.5,1.0,1.5").split(",")]
SEED, K, RES_SP = 0, 9, 1.0
OUT.mkdir(parents=True, exist_ok=True)
LOG = open(OUT / "exp4.log", "a")


def log(m):
    print(m, flush=True)
    LOG.write(m + "\n")
    LOG.flush()


def dec(a):
    return np.array([x.decode() if isinstance(x, bytes) else x for x in a], dtype=object)


def cat(f, col):
    return dec(f[f"obs/{col}/categories"][:])[f[f"obs/{col}/codes"][:]]


def load(src):
    kind, path, key = src.split(":", 2)
    if kind == "npz":
        e = np.load(path, allow_pickle=True)
        return pd.DataFrame(e[key], index=e["obs_names"].astype(str))
    f = h5py.File(path, "r")
    X = f[f"obsm/{key}"][:]
    smp, bc = cat(f, "sample"), cat(f, "original_spot")
    tables = {p.name.removesuffix("_spots.csv"): set(pd.read_csv(p, usecols=["barcode"]).barcode)
              for p in (DATA / "visium_tables").glob("*_spots.csv")}
    run_of = {}
    for s in np.unique(smp):
        b = set(bc[smp == s])
        hits = [r for r, t in tables.items() if t == b]
        assert len(hits) == 1, f"{s}: spot tables with the same barcodes: {hits}"
        run_of[s] = hits[0]
    log("sample -> run: " + ", ".join(f"{s}={r}" for s, r in sorted(run_of.items())))
    names = [f"{run_of[s]}_{b}" for s, b in zip(smp, bc)]
    return pd.DataFrame(X, index=names)


t0 = time.time()
E = load(SRC)
names = E.index
sample = names.str.split("_", n=1).str[0]
xy = []
for s in pd.unique(sample):
    d = pd.read_csv(DATA / "visium_tables" / f"{s}_spots.csv",
                    usecols=["barcode", "pxl_row_fullres", "pxl_col_fullres"])
    d.index = s + "_" + d.barcode
    xy.append(d)
xy = pd.concat(xy)
assert set(xy.index) == set(names), "embedding spots differ from the spot tables"
xy = xy.loc[names]
A = ad.AnnData(obs=pd.DataFrame({"sample": pd.Categorical(sample)}, index=names))
A.obsm["X_e"] = E.to_numpy(np.float32)
A.obsm["spatial"] = xy[["pxl_col_fullres", "pxl_row_fullres"]].to_numpy(float)
sc.pp.neighbors(A, n_neighbors=15, use_rep="X_e", random_state=0)
sq.gr.spatial_neighbors(A, library_key="sample", coord_type="grid", n_neighs=6)
log(f"{SRC}: {A.n_obs} spots, {A.obs['sample'].nunique()} slides, {time.time() - t0:.0f}s")


def run(res_lat, w):
    sl.spatialleiden(A, resolution=(res_lat, RES_SP), layer_ratio=w, directed=(False, False),
                     n_iterations=-1, random_state=SEED, key_added="t")
    return A.obs["t"].astype(str).to_numpy()


if RES == "search":  # step-halving, as in experiments 0 and 2
    res, step, n, i = 1.0, 0.1, len(set(run(1.0, 0))), 1
    while n != K and i < 30:
        sign = 1 if n < K else -1
        res += step * sign
        n = len(set(run(res, 0)))
        if (1 if n < K else -1) != sign:
            step /= 2
        i += 1
    log(f"searched expression resolution {res!r} -> K {n} at weight 0 ({i} runs)")
else:
    res = float(RES)

labels, rows = {}, []
for w in WEIGHTS:
    t = time.time()
    lab = run(res, w)
    labels[f"r{w}_s{SEED}"] = lab
    per = {s: len(set(lab[sample == s])) for s in pd.unique(sample)}
    rows.append({"weight": w, "res_expression": res, "res_spatial": RES_SP, "K": len(set(lab)),
                 **{f"K_{s}": k for s, k in per.items()}})
    log(f"weight {w}: K = {len(set(lab))} ({time.time() - t:.0f}s)")
pd.DataFrame(labels, index=names).to_csv(OUT / "exp2_labels.csv.gz")
pd.DataFrame(rows).to_csv(OUT / "exp4_k_per_weight.csv", index=False)
log(f"done in {time.time() - t0:.0f}s")
