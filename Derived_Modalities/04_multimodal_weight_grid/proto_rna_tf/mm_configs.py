"""Named multimodal SpatialLeiden configurations: RNA plus any of TF activity, CCC and
scCellFie metabolic tasks, plus space. Companion to mm_proto.py (RNA + TF weight grid).

All graphs come from experiment 0's loaders (exp0_sweep.py), so every layer is built the
same way as in the single-layer sweeps: RNA as in experiment 0; TF, CCC and metabolic
layers z-scored per feature, PCA (<= 30), 15-NN; Visium grid for space; all undirected.

Configurations (weights; RNA = 1 is the reference, extra layers 0.5, space 0.75 as in the
RNA + TF prototype figure):
    rna_tf        RNA 1, TF 0.5, space 0.75
    rna_ccc       RNA 1, CCC 0.5, space 0.75
    rna_met       RNA 1, metabolic 0.5, space 0.75
    rna_all       RNA 1, TF 0.5, CCC 0.5, metabolic 0.5, space 0.75
K = authors' cell-type niche count per slide; 5 seeds.
Holding K: the shared expression resolution (RNA and the extra layers) is searched once per
seed with all extra layers and space at weight 0 (= experiment-0 ratio-0 baseline); then the
spatial resolution is searched. If that cannot reach K (a dominant extra layer), the
expression resolution is searched with the spatial resolution fixed (fallback, recorded).

Options via environment: MM_SPATIAL_W (default 0.75), MM_CONFIGS (comma list, default all).
Usage (env spatialleiden, from Derived_Modalities/):
    python 04_multimodal_weight_grid/proto_rna_tf/mm_configs.py data/kuppe_2022 results/kuppe_2022/mm_configs [jobs]
"""

import importlib.util
import os
import sys
import time
from itertools import combinations
from multiprocessing import get_context
from pathlib import Path

HERE = Path(__file__).resolve()
EXP0 = HERE.parents[2] / "03_single_modality_clustering" / "exp0_layer_ratio_sweep" / "exp0_sweep.py"
spec = importlib.util.spec_from_file_location("exp0_sweep", EXP0)
e0 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(e0)  # sets one-thread BLAS limits before numpy loads

import anndata as ad  # noqa: E402
import mudata as md  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import spatialleiden as sl  # noqa: E402
from sklearn.metrics import adjusted_rand_score  # noqa: E402

SLIDES = ["ACH0013", "ACH0012", "10X001"]
SEEDS = [0, 1, 2, 3, 4]
SPATIAL_W = float(os.environ.get("MM_SPATIAL_W", 0.75))  # env vars reach the spawned workers
CONFIGS = {"rna_tf": {"tf": 0.5}, "rna_ccc": {"ccc": 0.5}, "rna_met": {"met": 0.5},
           "rna_all": {"tf": 0.5, "ccc": 0.5, "met": 0.5}}
if os.environ.get("MM_CONFIGS"):
    CONFIGS = {c: CONFIGS[c] for c in os.environ["MM_CONFIGS"].split(",")}
EXTRA = ["tf", "ccc", "met"]


def layer(data, run, name):
    if name == "rna":
        a = e0.load_slide(data, run)
        e0.preprocess(a)
    elif name == "tf":
        a = e0.load_tf_layer(data, run)
        e0.preprocess_tf(a)
    elif name == "ccc":
        a = e0.load_ccc_layer(data, run)
        e0.preprocess_ccc(a)
    else:
        a = e0.load_sccellfie_layer(data, run)
        e0.preprocess_ccc(a)
    return a


def build(data, run):
    layers = {n: layer(data, run, n) for n in ["rna"] + EXTRA}
    names = list(layers["rna"].obs_names)
    for n, a in layers.items():
        assert list(a.obs_names) == names, f"{run}: spot order differs in layer {n}"
    mods = {n: ad.AnnData(obs=pd.DataFrame(index=names), obsp={"connectivities": a.obsp["connectivities"]})
            for n, a in layers.items()}
    m = md.MuData(mods)
    m.obsp["spatial_connectivities"] = layers["rna"].obsp["spatial_connectivities"]
    return m, layers["rna"]


def cluster(m, res_expr, res_sp, weights, seed):
    mods = ["rna"] + EXTRA
    sl.spatialleiden_multimodal(
        m, resolution={**{n: res_expr for n in mods}, "spatial": res_sp},
        layer_weights={"rna": 1.0, **{n: float(weights.get(n, 0)) for n in EXTRA}, "spatial": float(weights["spatial"])},
        directed={**{n: False for n in mods}, "spatial": False},
        n_iterations=-1, random_state=seed, key_added="tmp")
    return m.obs["tmp"].astype(str).to_numpy()


def run_job(job):
    data, run, k, seed = job
    t0 = time.time()
    m, rna = build(data, run)
    ei, ej = e0.edges(rna)
    ct = rna.obs.celltype_niche.to_numpy()
    base_w = {"spatial": 0}
    r_expr, _ = e0.search_resolution_steps(lambda r: len(set(cluster(m, r, 1.0, base_w, seed))), k)
    rows, labels = [], {}
    for name, extra in CONFIGS.items():
        w = {**extra, "spatial": SPATIAL_W}
        r_sp, ok = e0.search_resolution_steps(lambda r: len(set(cluster(m, r_expr, r, w, seed))), k, start=0.4)
        r_e, fallback = r_expr, False
        if not ok:  # extra layers dominate: adjust the shared expression resolution instead
            r_e, ok = e0.search_resolution_steps(lambda r: len(set(cluster(m, r, r_sp, w, seed))), k, start=r_expr)
            fallback = True
        lab = cluster(m, r_e, r_sp, w, seed)
        labels[name] = lab
        obs = e0.coherence(lab, ei, ej)
        shuf = e0.shuffled_coherence(lab, ei, ej)
        rows.append({"slide": run, "seed": seed, "config": name, **{f"w_{n}": w.get(n, 0) for n in EXTRA},
                     "w_spatial": SPATIAL_W, "target_k": k, "n_clusters": len(set(lab)),
                     "reached_k": len(set(lab)) == k, "expression_fallback": fallback,
                     "res_expression": r_e, "res_spatial": r_sp, "coherence": obs, "coherence_shuffled": shuf,
                     "coherence_gap": obs - shuf, "ARI_celltype_niche": adjusted_rand_score(ct, lab)})
    return run, seed, rows, labels, list(rna.obs_names), f"{run} seed {seed}: {time.time() - t0:.0f}s"


def main():
    data, out = Path(sys.argv[1]), Path(sys.argv[2])
    out.mkdir(parents=True, exist_ok=True)
    k = {run: int(pd.read_csv(data / "visium_tables" / f"{run}_spots.csv", usecols=["celltype_niche"])
                  .celltype_niche.nunique()) for run in SLIDES}
    jobs = [(data, run, k[run], s) for run in SLIDES for s in SEEDS]
    res = []
    with get_context("spawn").Pool(int(sys.argv[3]) if len(sys.argv) > 3 else 8) as pool:
        for r in pool.imap_unordered(run_job, jobs):
            print(r[5], flush=True)
            res.append(r)
    df = pd.DataFrame([row for r in res for row in r[2]])
    stab = []
    for run in SLIDES:
        per = {r[1]: r[3] for r in res if r[0] == run}
        names = next(r[4] for r in res if r[0] == run)
        tab = {}
        for name in CONFIGS:
            ok = [s for s in SEEDS if len(set(per[s][name])) == k[run]]
            aris = [adjusted_rand_score(per[a][name], per[b][name]) for a, b in combinations(ok, 2)]
            stab.append({"slide": run, "config": name, "seed_ARI_mean": np.mean(aris) if aris else np.nan})
            for s in SEEDS:
                tab[f"{name}_s{s}"] = per[s][name]
        pd.DataFrame(tab, index=names).to_csv(out / f"labels_{run}.csv.gz")
    df = df.merge(pd.DataFrame(stab), on=["slide", "config"])
    df.to_csv(out / "mm_configs.csv", index=False)
    print(f"done: {len(df)} runs, {int((~df.reached_k).sum())} missed K, "
          f"{int(df.expression_fallback.sum())} used the expression-resolution fallback")


if __name__ == "__main__":
    main()
