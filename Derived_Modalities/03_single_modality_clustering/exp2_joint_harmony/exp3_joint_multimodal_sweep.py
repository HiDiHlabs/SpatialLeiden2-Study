"""Experiment 3: joint multimodal SpatialLeiden over all 27 Kuppe 2022 Visium slides on
Harmony-integrated RNA + a second layer (CCC or TF activity) + the spatial graph, sweeping the spatial layer's weight at K = 9.

Same design as experiment 2 (exp2_joint_sweep.py in this folder), with a second expression-type
layer:
  * layers (spatialleiden.spatialleiden_multimodal):
      rna      15-NN (scanpy, seed 0) on the RNA Harmony embedding (experiment 1, key "authors")
      ccc      15-NN on the CCC embedding (07_integration/ccc_harmony.py, key "X_harmony")
      spatial  squidpy grid graph, 6 neighbours, per slide (library_key = sample)
  * layer weights: rna 1, ccc W_CCC (default 0.5, as the rna_ccc configuration of
    04_multimodal_weight_grid/proto_rna_tf/mm_configs.py), spatial = the swept ratio
    0, 0.25, 0.5, 0.75, 1, 1.5, 2, 2.5 (the same values as experiment 2's layer_ratio).
  * directed False for all layers, n_iterations -1, seeds 0-4, K = 9.
  * holding K, as in mm_configs.py: per seed, the shared expression resolution (rna and ccc)
    is searched at spatial weight 0 with the multimodal call itself; at weight > 0 it stays
    fixed and the spatial resolution is searched (start 0.4). If that cannot reach K, the
    expression resolution is searched with the spatial resolution fixed (fallback, recorded
    in exp3_sweep.csv, column expression_fallback).

Outputs (OUT): exp2_labels.csv.gz (same name and format as experiment 2, so
01_data_preprocess/kuppe_fig1/fig1hij.py reads it unchanged), exp3_sweep.csv,
exp3_per_slide.csv, exp3_seed_stability.csv, exp3.log.

Usage (from Derived_Modalities/, spatialleiden env with mudata):
    python 03_single_modality_clustering/exp2_joint_harmony/exp3_joint_multimodal_sweep.py \
        data/kuppe_2022 \
        results/kuppe_2022/exp1_harmony_covariates/embeddings.npz:authors \
        ccc=results/kuppe_2022/ccc_harmony/ccc_embeddings.npz:X_harmony \
        results/kuppe_2022/exp3_joint_multimodal_rna_ccc [w_second=0.5] [workers=30]
RNA + TF: second argument tf=results/kuppe_2022/tf_harmony/tf_embeddings.npz:X_harmony.
Quick prototype (seed 0, four spatial weights): prefix EXP_SEEDS=0 EXP_RATIOS=0,0.5,1.0,1.5.
With fewer seeds, exp3_seed_stability.csv is empty.
"""

import os

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys  # noqa: E402
import time  # noqa: E402
from itertools import combinations  # noqa: E402
from multiprocessing import get_context  # noqa: E402
from pathlib import Path  # noqa: E402

import anndata as ad  # noqa: E402
import mudata as md  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scanpy as sc  # noqa: E402
import spatialleiden as sl  # noqa: E402
import squidpy as sq  # noqa: E402
from sklearn.metrics import adjusted_rand_score  # noqa: E402

DATA = Path(sys.argv[1])
RNA_EMB, RNA_KEY = sys.argv[2].rsplit(":", 1)
# second layer as name=path:key, e.g. ccc=results/.../ccc_embeddings.npz:X_harmony
#                               or tf=results/kuppe_2022/tf_harmony/tf_embeddings.npz:X_harmony
L2, _l2 = sys.argv[3].split("=", 1) if "=" in sys.argv[3] else ("ccc", sys.argv[3])
CCC_EMB, CCC_KEY = _l2.rsplit(":", 1)
OUT = Path(sys.argv[4])
W_CCC = float(sys.argv[5]) if len(sys.argv) > 5 else 0.5
N_WORKERS = int(sys.argv[6]) if len(sys.argv) > 6 else 30
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
MODS = ["rna", L2]
LOG = open(OUT / "exp3.log", "a")


def log(msg):
    print(msg, flush=True)
    LOG.write(msg + "\n")
    LOG.flush()


def search_resolution_steps(n_clusters_at, target, start=1.0, step=0.1, n_iter=N_SEARCH):
    """Step-halving search, identical to experiments 0 and 2. Returns (resolution, reached)."""
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


def load_emb(path, key, names=None):
    if path.startswith("h5ad:"):  # colleague's CCC h5ad (newer anndata): read via exp4_loader
        import importlib.util
        spec = importlib.util.spec_from_file_location(
            "exp4_loader", Path(__file__).with_name("exp4_loader.py"))
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        X = mod.load(f"{path}:{key}", DATA, print)
        return X if names is None else X.loc[names]
    e = np.load(path, allow_pickle=True)
    idx = pd.Index(e["obs_names"].astype(str))
    X = pd.DataFrame(e[key], index=idx)
    return X if names is None else X.loc[names]


# ---------- joint MuData, built once; workers inherit it via fork ----------
t0 = time.time()
rna = load_emb(RNA_EMB, RNA_KEY)
names = rna.index
ccc = load_emb(CCC_EMB, CCC_KEY)
assert set(ccc.index) == set(names), f"RNA and {L2} embeddings cover different spots"
ccc = ccc.loc[names]
sample = names.str.split("_", n=1).str[0]
frames = []
for s in pd.unique(sample):
    t = pd.read_csv(DATA / "visium_tables" / f"{s}_spots.csv",
                    usecols=["barcode", "pxl_row_fullres", "pxl_col_fullres",
                             "celltype_niche", "molecular_niche"])
    t.index = s + "_" + t.barcode
    frames.append(t)
spots = pd.concat(frames).loc[names]
obs = pd.DataFrame({"sample": pd.Categorical(sample),
                    "celltype_niche": spots.celltype_niche.to_numpy(),
                    "molecular_niche": spots.molecular_niche.to_numpy()}, index=names)
mods = {}
for n, X in [("rna", rna), (L2, ccc)]:
    a = ad.AnnData(obs=obs.copy())
    a.obsm["X_emb"] = X.to_numpy(np.float32)
    sc.pp.neighbors(a, n_neighbors=15, use_rep="X_emb", random_state=0)
    mods[n] = ad.AnnData(obs=pd.DataFrame(index=names),
                         obsp={"connectivities": a.obsp["connectivities"]})
g = ad.AnnData(obs=obs.copy())
g.obsm["spatial"] = spots[["pxl_col_fullres", "pxl_row_fullres"]].to_numpy(float)
sq.gr.spatial_neighbors(g, library_key="sample", coord_type="grid", n_neighs=6)
M = md.MuData(mods)
M.obsp["spatial_connectivities"] = g.obsp["spatial_connectivities"]
EI, EJ = g.obsp["spatial_connectivities"].nonzero()
SMP = obs["sample"].to_numpy()
assert (SMP[EI] == SMP[EJ]).all()
log(f"RNA {RNA_EMB}:{RNA_KEY}; {L2} {CCC_EMB}:{CCC_KEY}; w_{L2} {W_CCC}; seeds {SEEDS}; spatial {RATIOS}; {len(names)} spots, "
    f"{obs['sample'].nunique()} slides; {time.time() - t0:.0f}s")


def cluster(res_expr, res_sp, w_sp, seed):
    sl.spatialleiden_multimodal(
        M, resolution={"rna": res_expr, L2: res_expr, "spatial": res_sp},
        layer_weights={"rna": 1.0, L2: W_CCC, "spatial": float(w_sp)},
        directed={"rna": False, L2: False, "spatial": False},
        n_iterations=-1, random_state=seed, key_added="tmp")
    return M.obs["tmp"].astype(str).to_numpy()


def job_ratio0(seed):
    t = time.time()
    r_expr, ok = search_resolution_steps(lambda r: len(set(cluster(r, 1.0, 0, seed))), K)
    return seed, 0, r_expr, np.nan, not ok, False, cluster(r_expr, 1.0, 0, seed), time.time() - t


def job_ratio(args):
    seed, w, r_expr = args
    t = time.time()
    r_sp, ok = search_resolution_steps(lambda r: len(set(cluster(r_expr, r, w, seed))), K,
                                       start=0.4)
    fallback = False
    if not ok:  # expression layers dominate: adjust the shared expression resolution instead
        r_expr, ok = search_resolution_steps(lambda r: len(set(cluster(r, r_sp, w, seed))), K,
                                             start=r_expr)
        fallback = True
    return seed, w, r_expr, r_sp, not ok, fallback, cluster(r_expr, r_sp, w, seed), time.time() - t


def coherence_gap(lab, mask=None):
    em = np.ones(len(EI), bool) if mask is None else mask[EI]
    i, j = EI[em], EJ[em]
    obs_ = float((lab[i] == lab[j]).mean())
    rng = np.random.default_rng(0)
    shuf = []
    for _ in range(N_SHUFFLE):
        p = lab.copy()
        for s in np.unique(SMP[i]):
            m = SMP == s
            p[m] = rng.permutation(lab[m])
        shuf.append(float((p[i] == p[j]).mean()))
    return obs_, obs_ - float(np.mean(shuf))


if __name__ == "__main__":
    ctx = get_context("fork")
    results = []
    with ctx.Pool(len(SEEDS)) as pool:
        for res in pool.imap_unordered(job_ratio0, SEEDS):
            results.append(res)
            log(f"spatial 0 seed {res[0]}: res_expr {res[2]:.4f} failed={res[4]} "
                f"K={len(set(res[6]))} {res[7]:.0f}s")
    r0 = {r[0]: r[2] for r in results}
    jobs = [(s, w, r0[s]) for s in SEEDS for w in RATIOS if w > 0]
    with ctx.Pool(min(len(jobs), N_WORKERS)) as pool:
        for res in pool.imap_unordered(job_ratio, jobs):
            results.append(res)
            log(f"spatial {res[1]} seed {res[0]}: res=({res[2]:.4f}, {res[3]:.4f}) "
                f"failed={res[4]} fallback={res[5]} K={len(set(res[6]))} {res[7]:.0f}s")

    labels = {f"r{w}_s{s}": lab for s, w, _a, _b, _c, _d, lab, _t in results}
    cols = [f"r{w}_s{s}" for w in RATIOS for s in SEEDS]
    pd.DataFrame({c: labels[c] for c in cols}, index=names).to_csv(OUT / "exp2_labels.csv.gz")

    ct, mol = obs.celltype_niche.to_numpy(), obs.molecular_niche.to_numpy()
    rows, per_slide = [], []
    for seed, w, r_e, r_s, failed, fb, lab, secs in sorted(results, key=lambda x: (x[1], x[0])):
        sizes = np.sort(np.bincount(pd.factorize(lab)[0]))[::-1]
        coh, gap = coherence_gap(lab)
        base = labels[f"r0_s{seed}"]
        rows.append({"w_spatial": w, "second_layer": L2, "w_second": W_CCC, "seed": seed, "n_clusters": len(sizes),
                     "target_k": K, "search_failed": failed, "expression_fallback": fb,
                     "res_expression": r_e, "res_spatial": r_s, "coherence": coh,
                     "coherence_gap": gap, "ARI_to_spatial0": adjusted_rand_score(base, lab),
                     "ARI_celltype_niche": adjusted_rand_score(ct, lab),
                     "ARI_molecular_niche": adjusted_rand_score(mol, lab),
                     "ARI_sample": adjusted_rand_score(SMP, lab),
                     "largest_cluster_frac": sizes[0] / len(lab),
                     "cluster_sizes": ";".join(map(str, sizes)), "runtime_s": secs})
        for s in np.unique(SMP):
            m = SMP == s
            c_s, g_s = coherence_gap(lab, m)
            per_slide.append({"w_spatial": w, "seed": seed, "sample": s,
                              "n_domains_present": len(set(lab[m])),
                              "coherence": c_s, "coherence_gap": g_s,
                              "ARI_celltype_niche": adjusted_rand_score(ct[m], lab[m]),
                              "ARI_molecular_niche": adjusted_rand_score(mol[m], lab[m]),
                              "ARI_to_spatial0": adjusted_rand_score(base[m], lab[m])})
    pd.DataFrame(rows).to_csv(OUT / "exp3_sweep.csv", index=False)
    pd.DataFrame(per_slide).to_csv(OUT / "exp3_per_slide.csv", index=False)
    stab = []
    for w in RATIOS:
        aris = [adjusted_rand_score(labels[f"r{w}_s{a}"], labels[f"r{w}_s{b}"])
                for a, b in combinations(SEEDS, 2)]
        if not aris:  # one seed (quick prototype): no stability to report
            continue
        stab.append({"w_spatial": w, "seed_ARI_mean": np.mean(aris),
                     "seed_ARI_min": np.min(aris), "seed_ARI_max": np.max(aris)})
    pd.DataFrame(stab).to_csv(OUT / "exp3_seed_stability.csv", index=False)
    log(pd.DataFrame(rows).groupby("w_spatial")[
        ["n_clusters", "coherence_gap", "ARI_to_spatial0", "ARI_celltype_niche",
         "ARI_molecular_niche", "ARI_sample", "largest_cluster_frac"]].mean().round(3).to_string())
    if stab:
        log(pd.DataFrame(stab).round(3).to_string())
    log(f"done in {time.time() - t0:.0f}s")
