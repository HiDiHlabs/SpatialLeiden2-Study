"""Experiment 5: joint multimodal SpatialLeiden (all 27 slides) WITHOUT holding K fixed across
spatial weights; the multimodal counterpart of experiment 4 (exp4_unfixed_k.py).

Layers: RNA (weight 1) plus any of CCC / TF (weight 0.5 each, as in
04_multimodal_weight_grid/proto_rna_tf/mm_configs.py) plus the spatial grid (per slide).
The shared expression resolution is searched once so that spatial weight 0 gives K = 9
(step-halving, as in experiments 0-3); the spatial resolution is fixed at 1.0; then ONE run
per spatial weight (0, 0.5, 1, 1.5), seed 0, so K is free to change with the weight.

Layer sources are given as name=src, with src as in exp4_unfixed_k.py
(npz:<path>:<key> or h5ad:<path>:<obsm key>).

Outputs (OUT): exp2_labels.csv.gz (columns r<weight>_s0), exp5_k_per_weight.csv, exp5.log.

Usage (from Derived_Modalities/, spatialleiden env with mudata):
    python 03_single_modality_clustering/exp2_joint_harmony/exp5_unfixed_k_multimodal.py \
        data/kuppe_2022 results/kuppe_2022/exp5_unfixed_k_rna_ccc \
        rna=npz:results/kuppe_2022/exp1_harmony_covariates/embeddings.npz:authors \
        ccc=h5ad:/path/CCC_all27_TFmatched_Harmony.h5ad:X_pca_harmony_TFmatched
"""
import os

os.environ.setdefault("OMP_NUM_THREADS", "2")
import importlib.util  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import anndata as ad  # noqa: E402
import mudata as md  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scanpy as sc  # noqa: E402
import spatialleiden as sl  # noqa: E402
import squidpy as sq  # noqa: E402

DATA, OUT = Path(sys.argv[1]), Path(sys.argv[2])
LAYERS = dict(a.split("=", 1) for a in sys.argv[3:])
assert "rna" in LAYERS, "an rna=... layer is required"
W_EXTRA = float(os.environ.get("W_EXTRA", "0.5"))
WEIGHTS = [float(w) if float(w) else 0 for w in os.environ.get("EXP_RATIOS", "0,0.5,1.0,1.5").split(",")]
SEED, K, RES_SP = 0, 9, 1.0
OUT.mkdir(parents=True, exist_ok=True)
LOG = open(OUT / "exp5.log", "a")


def log(m):
    print(m, flush=True)
    LOG.write(m + "\n")
    LOG.flush()


# reuse exp4's loader (npz / h5ad via h5py, run mapping by barcode sets)
_spec = importlib.util.spec_from_file_location("exp4_loader", Path(__file__).with_name("exp4_loader.py"))
_l = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_l)

t0 = time.time()
emb = {n: _l.load(src, DATA, log) for n, src in LAYERS.items()}
names = emb["rna"].index
for n, e in emb.items():
    assert set(e.index) == set(names), f"{n}: spots differ from rna"
sample = names.str.split("_", n=1).str[0]
xy = []
for s in pd.unique(sample):
    d = pd.read_csv(DATA / "visium_tables" / f"{s}_spots.csv",
                    usecols=["barcode", "pxl_row_fullres", "pxl_col_fullres"])
    d.index = s + "_" + d.barcode
    xy.append(d)
xy = pd.concat(xy).loc[names]
mods = {}
for n, e in emb.items():
    a = ad.AnnData(obs=pd.DataFrame(index=names))
    a.obsm["X_e"] = e.loc[names].to_numpy(np.float32)
    sc.pp.neighbors(a, n_neighbors=15, use_rep="X_e", random_state=0)
    mods[n] = ad.AnnData(obs=pd.DataFrame(index=names), obsp={"connectivities": a.obsp["connectivities"]})
g = ad.AnnData(obs=pd.DataFrame({"sample": pd.Categorical(sample)}, index=names))
g.obsm["spatial"] = xy[["pxl_col_fullres", "pxl_row_fullres"]].to_numpy(float)
sq.gr.spatial_neighbors(g, library_key="sample", coord_type="grid", n_neighs=6)
M = md.MuData(mods)
M.obsp["spatial_connectivities"] = g.obsp["spatial_connectivities"]
EXPR = list(emb)
log(f"layers {LAYERS}; extra weight {W_EXTRA}; {len(names)} spots; {time.time() - t0:.0f}s")


def run(res, w):
    sl.spatialleiden_multimodal(
        M, resolution={**{n: res for n in EXPR}, "spatial": RES_SP},
        layer_weights={**{n: (1.0 if n == "rna" else W_EXTRA) for n in EXPR}, "spatial": float(w)},
        directed={**{n: False for n in EXPR}, "spatial": False},
        n_iterations=-1, random_state=SEED, key_added="t")
    return M.obs["t"].astype(str).to_numpy()


if os.environ.get("RES_EXPR"):  # expression resolution already known (e.g. from a parallel scan)
    res, step, n, i = float(os.environ["RES_EXPR"]), 0.1, K, 30
else:
    res, step, n, i = 1.0, 0.1, len(set(run(1.0, 0))), 1
while n != K and i < 30:
    sign = 1 if n < K else -1
    res += step * sign
    n = len(set(run(res, 0)))
    if (1 if n < K else -1) != sign:
        step /= 2
    i += 1
log(f"searched expression resolution {res!r} -> K {n} at weight 0 ({i} runs)")

def one(w):
    t = time.time()
    return w, run(res, w), time.time() - t


labels, rows = {}, []
from multiprocessing import get_context  # noqa: E402
with get_context("fork").Pool(len(WEIGHTS)) as pool:  # weights in parallel
    out = sorted(pool.map(one, WEIGHTS), key=lambda r: r[0])
for w, lab, secs in out:
    labels[f"r{w}_s{SEED}"] = lab
    per = {s: len(set(lab[sample == s])) for s in pd.unique(sample)}
    rows.append({"weight": w, "res_expression": res, "res_spatial": RES_SP, "K": len(set(lab)),
                 **{f"K_{s}": k for s, k in per.items()}})
    log(f"weight {w}: K = {len(set(lab))} ({secs:.0f}s)")
pd.DataFrame(labels, index=names).to_csv(OUT / "exp2_labels.csv.gz")
pd.DataFrame(rows).to_csv(OUT / "exp5_k_per_weight.csv", index=False)
log(f"done in {time.time() - t0:.0f}s")
