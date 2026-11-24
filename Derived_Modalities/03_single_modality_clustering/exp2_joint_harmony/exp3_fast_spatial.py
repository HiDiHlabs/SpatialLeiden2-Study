"""Fast finish of experiment 3 (multimodal RNA + second layer + space, K = 9, seed 0).

Same model and search rule as exp3_joint_multimodal_sweep.py, but for a deadline:
the shared expression resolution at spatial weight 0 is taken from that script's log
(argument RES_EXPR), and each spatial-resolution search starts from a given value
(default: the seed-0 values RNA alone needed for K = 9 in experiment 2) with step 0.05,
instead of 0.4 with step 0.1. One process per weight, in parallel.

Writes exp2_labels.csv.gz (columns r<weight>_s0) and exp3_fast.csv into OUT.

Usage (from Derived_Modalities/):
    python 03_single_modality_clustering/exp2_joint_harmony/exp3_fast_spatial.py data/kuppe_2022 \
        results/kuppe_2022/exp1_harmony_covariates/embeddings.npz:authors \
        ccc=h5ad:/home/gage11/CCC_all27_TFmatched_Harmony.h5ad:X_pca_harmony_TFmatched \
        results/kuppe_2022/exp3_joint_multimodal_rna_ccc <RES_EXPR> [0.5]
"""
import os

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import importlib.util  # noqa: E402
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

DATA = Path(sys.argv[1])
RNA_SRC = "npz:" + sys.argv[2]
L2, L2_SRC = sys.argv[3].split("=", 1)
L2_SRC = L2_SRC if L2_SRC.startswith(("npz:", "h5ad:")) else "npz:" + L2_SRC
OUT = Path(sys.argv[4])
RES_EXPR = float(sys.argv[5])
W2 = float(sys.argv[6]) if len(sys.argv) > 6 else 0.5
K, SEED = 9, 0
START = {0.5: 1.5, 1.0: 0.98125, 1.5: 0.65}  # RNA-only seed-0 spatial resolutions (exp2)
OUT.mkdir(parents=True, exist_ok=True)

spec = importlib.util.spec_from_file_location("exp4_loader", Path(__file__).with_name("exp4_loader.py"))
ldr = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ldr)

t0 = time.time()
rna = ldr.load(RNA_SRC, DATA)
names = rna.index
l2 = ldr.load(L2_SRC, DATA).loc[names]
sample = names.str.split("_", n=1).str[0]
xy = []
for s in pd.unique(sample):
    d = pd.read_csv(DATA / "visium_tables" / f"{s}_spots.csv",
                    usecols=["barcode", "pxl_row_fullres", "pxl_col_fullres"])
    d.index = s + "_" + d.barcode
    xy.append(d)
xy = pd.concat(xy).loc[names]
mods = {}
for n, X in [("rna", rna), (L2, l2)]:
    a = ad.AnnData(obs=pd.DataFrame(index=names))
    a.obsm["X_e"] = X.to_numpy(np.float32)
    sc.pp.neighbors(a, n_neighbors=15, use_rep="X_e", random_state=0)
    mods[n] = ad.AnnData(obs=pd.DataFrame(index=names), obsp={"connectivities": a.obsp["connectivities"]})
g = ad.AnnData(obs=pd.DataFrame({"sample": pd.Categorical(sample)}, index=names))
g.obsm["spatial"] = xy[["pxl_col_fullres", "pxl_row_fullres"]].to_numpy(float)
sq.gr.spatial_neighbors(g, library_key="sample", coord_type="grid", n_neighs=6)
M = md.MuData(mods)
M.obsp["spatial_connectivities"] = g.obsp["spatial_connectivities"]
print(f"graphs built {time.time() - t0:.0f}s", flush=True)


def cluster(res_sp, w):
    sl.spatialleiden_multimodal(
        M, resolution={"rna": RES_EXPR, L2: RES_EXPR, "spatial": res_sp},
        layer_weights={"rna": 1.0, L2: W2, "spatial": float(w)},
        directed={"rna": False, L2: False, "spatial": False},
        n_iterations=-1, random_state=SEED, key_added="t")
    return M.obs["t"].astype(str).to_numpy()


def job(w):
    t = time.time()
    if w == 0:
        lab = cluster(1.0, 0)
        return w, 1.0, len(set(lab)) == K, lab, time.time() - t, 1
    res, step = START.get(w, 0.4), 0.05
    lab = cluster(res, w)
    n, i = len(set(lab)), 1
    while n != K and i < 30:
        sign = 1 if n < K else -1
        res += step * sign
        lab = cluster(res, w)
        n = len(set(lab))
        if (1 if n < K else -1) != sign:
            step /= 2
        i += 1
    return w, res, n == K, lab, time.time() - t, i


def scan(res):
    """K at spatial weight 0 for one expression resolution (FAST_SCAN mode)."""
    global RES_EXPR
    RES_EXPR = res
    return res, len(set(cluster(1.0, 0)))


if __name__ == "__main__" and os.environ.get("FAST_SCAN"):
    # parallel scan of the shared expression resolution at weight 0 (instead of a serial search)
    cand = [float(x) for x in os.environ["FAST_SCAN"].split(",")]
    with get_context("fork").Pool(len(cand)) as pool:
        for r, k in pool.imap_unordered(scan, cand):
            print(f"scan res_expr {r}: K = {k}", flush=True)
    sys.exit(0)

if __name__ == "__main__":
    with get_context("fork").Pool(4) as pool:
        res = []
        for r in pool.imap_unordered(job, [0, 0.5, 1.0, 1.5]):
            print(f"weight {r[0]}: res_spatial {r[1]!r} reached K={r[2]} ({r[5]} runs, {r[4]:.0f}s)", flush=True)
            res.append(r)
    res.sort(key=lambda r: r[0])
    pd.DataFrame({f"r{w}_s0": lab for w, _, _, lab, _, _ in res}, index=names).to_csv(OUT / "exp2_labels.csv.gz")
    pd.DataFrame([{"weight": w, "res_expression": RES_EXPR, "res_spatial": r, "reached_k": ok,
                   "n_runs": n, "seconds": round(s)} for w, r, ok, _, s, n in res]).to_csv(OUT / "exp3_fast.csv", index=False)
    print(f"done {time.time() - t0:.0f}s", flush=True)
