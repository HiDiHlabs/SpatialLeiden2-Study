"""Systematic grid of joint multimodal SpatialLeiden runs over all 27 Kuppe slides (strategy section 4.1 of
handoff/HANDOFF_2026-11-25_joint_multimodal_findings_and_strategy.md). Resumable: every run is one file,
existing files are skipped, so the script can be stopped and restarted at any time.

Layers (as in the joint runs of 24 Nov): 15-NN graphs on the Harmony embeddings of
  rna  results/kuppe_2022/exp1_harmony_covariates/embeddings.npz ("authors")
  tf   results/kuppe_2022/tf_harmony/tf_embeddings.npz ("X_harmony")
  ccc  results/kuppe_2022/ccc_harmony_ours/ccc_embeddings.npz ("X_harmony"; our LIANA+ layer)
plus the Visium grid (6 neighbours, per slide). All undirected, n_iterations -1.
Each configuration builds a MuData with ONLY its layers (a layer at weight 0 still changes the optimiser).

Grid
  configs   rna; rna+tf, rna+ccc, rna+tf+ccc with extra-layer weight 0.25 / 0.5 / 1 (RNA = 1)  -> 10
  target K  9, 15, 25, 40 (sets the expression resolution; K itself is then free)
  spatial   0, 0.25, 0.75, 1.2, 2
  seeds     0-4 (GRID_SEEDS to restrict)
Phase "calibrate": per (config, K), seed 0, search the shared expression resolution that gives exactly K
  domains at spatial weight 0 (step-halving search as in experiments 0-3) -> calibration.csv.
Phase "grid": per (config, K, spatial weight, seed) ONE run with that expression resolution and the
  spatial resolution fixed at 1.0 (package default; the preprint's discovery mode) -> runs/<id>.npz
  (int16 labels in the fixed spot order of spots.csv.gz) and a row in runs.csv.

Usage (from Derived_Modalities/, env spatialleiden_joint):
    python 08_grid_search/grid_search.py data/kuppe_2022 results/kuppe_2022 results/kuppe_2022/grid_search calibrate [workers]
    python 08_grid_search/grid_search.py data/kuppe_2022 results/kuppe_2022 results/kuppe_2022/grid_search grid [workers]
Environment: GRID_SEEDS (e.g. "0"), GRID_K (e.g. "9,15"), GRID_CONFIGS (e.g. "rna,rna+tf_0.5"), GRID_WEIGHTS.
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

DATA, RES, OUT, PHASE = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]), sys.argv[4]
N_WORKERS = int(sys.argv[5]) if len(sys.argv) > 5 else 28
(OUT / "runs").mkdir(parents=True, exist_ok=True)
LAYERS = {"rna": ("exp1_harmony_covariates/embeddings.npz", "authors"),
          "tf": ("tf_harmony/tf_embeddings.npz", "X_harmony"),
          "ccc": ("ccc_harmony_ours/ccc_embeddings.npz", "X_harmony")}
CONFIGS = {"rna": {"rna": 1.0}}
for w in (0.25, 0.5, 1.0):
    CONFIGS[f"rna+tf_{w}"] = {"rna": 1.0, "tf": w}
    CONFIGS[f"rna+ccc_{w}"] = {"rna": 1.0, "ccc": w}
    CONFIGS[f"rna+tf+ccc_{w}"] = {"rna": 1.0, "tf": w, "ccc": w}


def env_list(name, default, cast):
    v = os.environ.get(name)
    return [cast(x) for x in v.split(",")] if v else default


K_TARGETS = env_list("GRID_K", [9, 15, 25, 40], int)
WEIGHTS = env_list("GRID_WEIGHTS", [0.0, 0.25, 0.75, 1.2, 2.0], float)
SEEDS = env_list("GRID_SEEDS", [0, 1, 2, 3, 4], int)
CFG_NAMES = env_list("GRID_CONFIGS", list(CONFIGS), str)
RES_SP = 1.0
LOG = open(OUT / f"{PHASE}.log", "a")


def log(msg):
    print(msg, flush=True)
    LOG.write(time.strftime("%H:%M:%S ") + msg + "\n")
    LOG.flush()


# ---------- graphs, built once in the parent; workers inherit them via fork ----------
t0 = time.time()
emb, names = {}, None
for k, (f, key) in LAYERS.items():
    e = np.load(RES / f, allow_pickle=True)
    idx = pd.Index(e["obs_names"].astype(str))
    names = idx if names is None else names
    emb[k] = pd.DataFrame(e[key], index=idx).loc[names].to_numpy(np.float32)
sample = names.str.split("_", n=1).str[0]
frames = []
for s in pd.unique(sample):
    t = pd.read_csv(DATA / "visium_tables" / f"{s}_spots.csv",
                    usecols=["barcode", "pxl_row_fullres", "pxl_col_fullres", "major_labl", "patient",
                             "celltype_niche", "molecular_niche"])
    t.index = s + "_" + t.barcode
    frames.append(t)
spots = pd.concat(frames).loc[names]
if not (OUT / "spots.csv.gz").exists():
    spots.assign(sample=sample).drop(columns=["barcode"]).to_csv(OUT / "spots.csv.gz")
GRAPHS = {}
for k, X in emb.items():
    a = ad.AnnData(obs=pd.DataFrame(index=names))
    a.obsm["X_e"] = X
    sc.pp.neighbors(a, n_neighbors=15, use_rep="X_e", random_state=0)
    GRAPHS[k] = a.obsp["connectivities"]
g = ad.AnnData(obs=pd.DataFrame({"sample": pd.Categorical(sample)}, index=names))
g.obsm["spatial"] = spots[["pxl_col_fullres", "pxl_row_fullres"]].to_numpy(float)
sq.gr.spatial_neighbors(g, library_key="sample", coord_type="grid", n_neighs=6)
SPATIAL = g.obsp["spatial_connectivities"]
log(f"phase {PHASE}: graphs built for {len(names)} spots in {time.time() - t0:.0f}s")


def mudata(cfg):
    m = md.MuData({k: ad.AnnData(obs=pd.DataFrame(index=names), obsp={"connectivities": GRAPHS[k]}) for k in CONFIGS[cfg]})
    m.obsp["spatial_connectivities"] = SPATIAL
    return m


def cluster(m, cfg, res_expr, w_sp, seed):
    lw = CONFIGS[cfg]
    sl.spatialleiden_multimodal(
        m, resolution={**{k: res_expr for k in lw}, "spatial": RES_SP},
        layer_weights={**lw, "spatial": float(w_sp)}, directed={**{k: False for k in lw}, "spatial": False},
        n_iterations=-1, random_state=seed, key_added="t")
    return m.obs["t"].astype(int).to_numpy()


def search(n_at, target, start=1.0, step=0.1, n_iter=30):
    res, n, i = start, n_at(start), 1
    while n != target:
        if i >= n_iter:
            return res, n, False
        sign = 1 if n < target else -1
        res += step * sign
        if res <= 0:
            res, step = step / 2, step / 2
        n = n_at(res)
        if (1 if n < target else -1) != sign:
            step /= 2
        i += 1
    return res, n, True


def calib_job(args):
    cfg, k = args
    t = time.time()
    m = mudata(cfg)
    start = {9: 0.8, 15: 1.2, 25: 1.8, 40: 2.6}.get(k, 1.0)
    res, n, ok = search(lambda r: len(set(cluster(m, cfg, r, 0.0, 0))), k, start=start, step=0.2)
    return {"config": cfg, "target_k": k, "res_expression": res, "k_reached": n, "reached": ok,
            "seconds": round(time.time() - t)}


def grid_job(args):
    cfg, k, res, w, seed, rid = args
    f = OUT / "runs" / f"{rid}.npz"
    if f.exists():
        return None
    t = time.time()
    lab = cluster(mudata(cfg), cfg, res, w, seed)
    np.savez_compressed(f, labels=lab.astype(np.int16))
    return {"run_id": rid, "config": cfg, "target_k": k, "res_expression": res, "res_spatial": RES_SP,
            "spatial_weight": w, "seed": seed, "n_domains": int(len(set(lab))), "seconds": round(time.time() - t)}


if PHASE == "calibrate":
    cal_f = OUT / "calibration.csv"
    done = pd.read_csv(cal_f) if cal_f.exists() else pd.DataFrame(columns=["config", "target_k"])
    todo = [(c, k) for c in CFG_NAMES for k in K_TARGETS
            if not ((done.config == c) & (done.target_k == k)).any()]
    log(f"calibration jobs: {len(todo)}")
    with get_context("fork").Pool(min(N_WORKERS, max(len(todo), 1))) as pool:
        for r in pool.imap_unordered(calib_job, todo):
            log(f"calibrated {r['config']} K={r['target_k']}: res {r['res_expression']:.4f} -> {r['k_reached']} "
                f"({'ok' if r['reached'] else 'NOT reached'}, {r['seconds']}s)")
            pd.DataFrame([r]).to_csv(cal_f, mode="a", header=not cal_f.exists(), index=False)
elif PHASE == "grid":
    cal = pd.read_csv(OUT / "calibration.csv").drop_duplicates(["config", "target_k"], keep="last")
    jobs = []
    for _, c in cal.iterrows():
        if c.config not in CFG_NAMES or c.target_k not in K_TARGETS:
            continue
        for w in WEIGHTS:
            for s in SEEDS:
                rid = f"{c.config}__K{c.target_k}__w{w}__s{s}"
                if not (OUT / "runs" / f"{rid}.npz").exists():
                    jobs.append((c.config, int(c.target_k), float(c.res_expression), w, s, rid))
    # largest graphs (3 layers) and highest K first, so the slowest jobs do not end up last
    jobs.sort(key=lambda j: (-len(CONFIGS[j[0]]), -j[1]))
    log(f"grid jobs to run: {len(jobs)}")
    runs_f = OUT / "runs.csv"
    with get_context("fork").Pool(N_WORKERS) as pool:
        for i, r in enumerate(pool.imap_unordered(grid_job, jobs), 1):
            if r is None:
                continue
            pd.DataFrame([r]).to_csv(runs_f, mode="a", header=not runs_f.exists(), index=False)
            if i % 10 == 0 or i == len(jobs):
                log(f"{i}/{len(jobs)} done; last {r['run_id']} -> {r['n_domains']} domains ({r['seconds']}s)")
log("phase finished")
