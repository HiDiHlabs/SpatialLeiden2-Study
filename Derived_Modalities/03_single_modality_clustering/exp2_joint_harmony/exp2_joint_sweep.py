"""Experiment 2: joint SpatialLeiden over all 27 Kuppe 2022 Visium slides on Harmony-integrated
RNA, sweeping the spatial layer's weight (`layer_ratio`) at fixed K = 9.

Same question as experiment 0 (03_single_modality_clustering/exp0_layer_ratio_sweep/), now on
the integrated expression space, so the domains are shared across slides:
  * Expression graph: 15-NN (scanpy, seed 0) on the 30 Harmony PCs with the authors'
    covariates (patient + sample + batch), from experiment 1
    (results/kuppe_2022/exp1_harmony_covariates/embeddings.npz, key "authors").
  * Spatial graph: squidpy grid graph, 6 neighbours, built per slide (library_key = sample),
    so no spatial edge crosses slides. As in the preprint's
    clustering/spatialleiden_multisample.py (--use_batch, grid).
  * directed (False, False), n_iterations -1, as in experiment 0.
  * K = 9 (the authors' cell-type niches, pooled over 27 slides; also exp1's K).
  * Holding K fixed, as in experiment 0: at ratio 0 the expression resolution is searched with
    spatialleiden(layer_ratio=0) itself (per seed); at ratio > 0 it stays at that value and
    the spatial resolution is searched (search_resolution_spatial), up to 30 iterations.
  * layer_ratio 0, 0.25, 0.5, 0.75, 1, 1.5, 2, 2.5; seeds 0-4.

Outputs (results/kuppe_2022/exp2_joint_harmony/):
  exp2_labels.csv.gz   spots x runs (columns r<ratio>_s<seed>), index "<sample>_<barcode>"
  exp2_sweep.csv       one row per (ratio, seed): K, resolutions, search outcome, pooled
                       coherence gap, ARI to ratio 0 and to the authors' niches, sizes
  exp2_per_slide.csv   one row per (ratio, seed, slide): coherence gap, ARI to niches,
                       number of domains present on the slide
  exp2_seed_stability.csv  pairwise ARI among seeds per ratio
  exp2.log

Usage (from Derived_Modalities/, spatialleiden env; args 4-5 optional: embedding key, workers):
    python 03_single_modality_clustering/exp2_joint_harmony/exp2_joint_sweep.py \
        data/kuppe_2022 results/kuppe_2022/exp1_harmony_covariates/embeddings.npz \
        results/kuppe_2022/exp2_joint_harmony
"""

import os

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_v, "1")  # one thread per job; parallelism from the job pool
import sys  # noqa: E402
import time  # noqa: E402
import warnings  # noqa: E402
from itertools import combinations  # noqa: E402
from multiprocessing import get_context  # noqa: E402
from pathlib import Path  # noqa: E402

import anndata as ad  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scanpy as sc  # noqa: E402
import spatialleiden as sl  # noqa: E402
import squidpy as sq  # noqa: E402
from sklearn.metrics import adjusted_rand_score  # noqa: E402

DATA, EMB, OUT = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
# embedding key in EMB: "authors" for RNA (exp1), "X_harmony" for TF activity
# (results/kuppe_2022/tf_harmony/tf_embeddings.npz, 07_integration/RESULT_tf_harmony.md)
EMB_KEY = sys.argv[4] if len(sys.argv) > 4 else "authors"
N_WORKERS = int(sys.argv[5]) if len(sys.argv) > 5 else 30
OUT.mkdir(parents=True, exist_ok=True)
K = 9
RATIOS = [0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5]
SEEDS = [0, 1, 2, 3, 4]
# quick prototype: EXP_SEEDS=0 EXP_RATIOS=0,0.5,1.0,1.5
if os.environ.get("EXP_SEEDS"):
    SEEDS = [int(x) for x in os.environ["EXP_SEEDS"].split(",")]
if os.environ.get("EXP_RATIOS"):
    RATIOS = [float(x) if float(x) else 0 for x in os.environ["EXP_RATIOS"].split(",")]
N_SEARCH = 30
N_SHUFFLE = 20
DIRECTED = (False, False)
LOG = open(OUT / "exp2.log", "a")


def log(msg):
    print(msg, flush=True)
    LOG.write(msg + "\n")
    LOG.flush()


def search_resolution_steps(n_clusters_at, target, start=1.0, step=0.1, n_iter=N_SEARCH):
    """Step-halving search, identical to experiment 0 (spatialleiden 0.4's private
    _search_resolution). Returns (resolution, reached_target)."""
    res = start
    n = n_clusters_at(res)
    i = 1
    while n != target:
        if i >= n_iter:
            return res, False
        sign = 1 if n < target else -1
        res += step * sign
        n = n_clusters_at(res)
        if (1 if n < target else -1) != sign:
            step /= 2
        i += 1
    return res, True


# ---------- build the joint object once; worker processes inherit it via fork ----------
t0 = time.time()
e = np.load(EMB, allow_pickle=True)
names = pd.Index(e["obs_names"].astype(str))
sample = names.str.split("_", n=1).str[0]
barcode = names.str.split("_", n=1).str[1]
frames = []
for s in pd.unique(sample):
    t = pd.read_csv(DATA / "visium_tables" / f"{s}_spots.csv",
                    usecols=["barcode", "pxl_row_fullres", "pxl_col_fullres",
                             "celltype_niche", "molecular_niche", "major_labl"])
    t.index = s + "_" + t.barcode
    frames.append(t)
spots = pd.concat(frames)
assert names.isin(spots.index).all()
spots = spots.loc[names]
A = ad.AnnData(obs=pd.DataFrame({"sample": pd.Categorical(sample), "barcode": barcode,
                                 "zone": spots.major_labl.to_numpy(),
                                 "celltype_niche": spots.celltype_niche.to_numpy(),
                                 "molecular_niche": spots.molecular_niche.to_numpy()},
                                index=names))
A.obsm["X_harmony"] = e[EMB_KEY]
A.obsm["spatial"] = spots[["pxl_col_fullres", "pxl_row_fullres"]].to_numpy(float)
sc.pp.neighbors(A, n_neighbors=15, use_rep="X_harmony", random_state=0)
sq.gr.spatial_neighbors(A, library_key="sample", coord_type="grid", n_neighs=6)
EI, EJ = A.obsp["spatial_connectivities"].nonzero()
assert (A.obs["sample"].to_numpy()[EI] == A.obs["sample"].to_numpy()[EJ]).all()
log(f"embedding {EMB}:{EMB_KEY}; joint object: {A.n_obs} spots, {A.obs['sample'].nunique()} slides, "
    f"{A.obsp['connectivities'].nnz} kNN entries, {len(EI)} spatial entries, "
    f"{time.time() - t0:.0f}s")


def run(resolution, ratio, seed):
    # each worker process has its own (forked, copy-on-write) A, so "tmp" never collides
    sl.spatialleiden(A, resolution=resolution, layer_ratio=ratio, directed=DIRECTED,
                     n_iterations=-1, random_state=seed, key_added="tmp")
    return A.obs["tmp"].astype(str).to_numpy()


def job_ratio0(seed):
    t = time.time()
    r_lat, ok = search_resolution_steps(lambda r: len(set(run((r, 1.0), 0, seed))), K)
    lab = run((r_lat, 1.0), 0, seed)
    return seed, 0, r_lat, np.nan, not ok, lab, time.time() - t


def job_ratio(args):
    seed, ratio, r_lat = args
    t = time.time()
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        r_sp = sl.search_resolution_spatial(
            A, K, n_iter=N_SEARCH, resolution=(r_lat, 1.0), layer_ratio=ratio,
            directed=DIRECTED, n_iterations=-1, random_state=seed, key_added="tmp")
    failed = any("Correct resolution not found" in str(x.message) for x in w)
    lab = run((r_lat, r_sp), ratio, seed)  # stored labels are exactly this run
    return seed, ratio, r_lat, r_sp, failed, lab, time.time() - t


def coherence_gap(lab, i, j):
    obs = float((lab[i] == lab[j]).mean())
    rng = np.random.default_rng(0)
    # shuffle within slide: permute labels separately per sample
    smp = A.obs["sample"].to_numpy()
    shuf = []
    for _ in range(N_SHUFFLE):
        p = lab.copy()
        for s in np.unique(smp[i]):
            m = smp == s
            p[m] = rng.permutation(lab[m])
        shuf.append(float((p[i] == p[j]).mean()))
    return obs, obs - float(np.mean(shuf))


if __name__ == "__main__":
    ctx = get_context("fork")
    labels, results = {}, []
    with ctx.Pool(len(SEEDS)) as pool:
        for res in pool.imap_unordered(job_ratio0, SEEDS):
            results.append(res)
            log(f"ratio 0 seed {res[0]}: res_latent {res[2]:.4f} failed={res[4]} "
                f"K={len(set(res[5]))} {res[6]:.0f}s")
    r_lat_by_seed = {r[0]: r[2] for r in results}
    jobs = [(s, r, r_lat_by_seed[s]) for s in SEEDS for r in RATIOS if r > 0]
    with ctx.Pool(min(len(jobs), N_WORKERS)) as pool:
        for res in pool.imap_unordered(job_ratio, jobs):
            results.append(res)
            log(f"ratio {res[1]} seed {res[0]}: res=({res[2]:.4f}, {res[3]:.4f}) "
                f"failed={res[4]} K={len(set(res[5]))} {res[6]:.0f}s")

    for seed, ratio, *_rest in results:
        labels[f"r{ratio}_s{seed}"] = _rest[3]
    cols = [f"r{r}_s{s}" for r in RATIOS for s in SEEDS]
    lab_df = pd.DataFrame({c: labels[c] for c in cols}, index=A.obs_names)
    lab_df.to_csv(OUT / "exp2_labels.csv.gz")

    ct, mol = A.obs.celltype_niche.to_numpy(), A.obs.molecular_niche.to_numpy()
    smp = A.obs["sample"].to_numpy()
    rows, per_slide = [], []
    for seed, ratio, r_lat, r_sp, failed, lab, secs in sorted(results, key=lambda x: (x[1], x[0])):
        sizes = np.sort(np.bincount(pd.factorize(lab)[0]))[::-1]
        coh, gap = coherence_gap(lab, EI, EJ)
        rows.append({
            "layer_ratio": ratio, "seed": seed, "n_clusters": len(sizes), "target_k": K,
            "search_failed": failed, "res_latent": r_lat, "res_spatial": r_sp,
            "coherence": coh, "coherence_gap": gap,
            "ARI_to_ratio0": adjusted_rand_score(labels[f"r0_s{seed}"], lab),
            "ARI_celltype_niche": adjusted_rand_score(ct, lab),
            "ARI_molecular_niche": adjusted_rand_score(mol, lab),
            "ARI_sample": adjusted_rand_score(smp, lab),
            "largest_cluster_frac": sizes[0] / len(lab),
            "smallest_cluster_frac": sizes[-1] / len(lab),
            "cluster_sizes": ";".join(map(str, sizes)), "runtime_s": secs,
        })
        for s in np.unique(smp):
            m = smp == s
            em = m[EI]
            same = lab[EI[em]] == lab[EJ[em]]
            rng = np.random.default_rng(0)
            sub = lab[m]
            # within-slide shuffled coherence
            idx = np.flatnonzero(m)
            pos = np.full(len(lab), -1)
            pos[idx] = np.arange(len(idx))
            ii, jj = pos[EI[em]], pos[EJ[em]]
            sh = np.mean([(lambda p: (p[ii] == p[jj]).mean())(rng.permutation(sub))
                          for _ in range(N_SHUFFLE)])
            per_slide.append({
                "layer_ratio": ratio, "seed": seed, "sample": s,
                "n_domains_present": len(set(sub)),
                "coherence": same.mean(), "coherence_gap": same.mean() - sh,
                "ARI_celltype_niche": adjusted_rand_score(ct[m], sub),
                "ARI_molecular_niche": adjusted_rand_score(mol[m], sub),
                "ARI_to_ratio0": adjusted_rand_score(labels[f"r0_s{seed}"][m], sub),
            })
    pd.DataFrame(rows).to_csv(OUT / "exp2_sweep.csv", index=False)
    pd.DataFrame(per_slide).to_csv(OUT / "exp2_per_slide.csv", index=False)
    stab = []
    for r in RATIOS:
        aris = [adjusted_rand_score(labels[f"r{r}_s{a}"], labels[f"r{r}_s{b}"])
                for a, b in combinations(SEEDS, 2)]
        if not aris:  # one seed (quick prototype)
            continue
        stab.append({"layer_ratio": r, "seed_ARI_mean": np.mean(aris),
                     "seed_ARI_min": np.min(aris), "seed_ARI_max": np.max(aris)})
    pd.DataFrame(stab).to_csv(OUT / "exp2_seed_stability.csv", index=False)
    log(pd.DataFrame(rows).groupby("layer_ratio")[
        ["n_clusters", "coherence_gap", "ARI_to_ratio0", "ARI_celltype_niche",
         "ARI_molecular_niche", "ARI_sample", "largest_cluster_frac"]].mean().round(3).to_string())
    if stab:
        log(pd.DataFrame(stab).round(3).to_string())
    log(f"done in {time.time() - t0:.0f}s")
