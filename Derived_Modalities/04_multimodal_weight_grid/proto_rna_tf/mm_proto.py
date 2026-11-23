"""First multimodal SpatialLeiden prototype: RNA + TF activity (CollecTRI) + space.

Uses spatialleiden.spatialleiden_multimodal on a MuData with two modalities, "rna" and
"tf", plus the spatial grid graph. All loading, graphs, coherence and the resolution
search are imported from experiment 0 (exp0_sweep.py), so the settings match exactly:
    rna graph: genes >=10 spots, 2,000 HVGs (seurat_v3), normalize_total, log1p,
               scale w/o centring, PCA 30, 15-NN
    tf graph:  CollecTRI ULM scores, z-score per TF, PCA 30, 15-NN
    spatial:   Visium grid, 6 neighbours; all graphs undirected
Grid: RNA weight 1; spatial weight in {0.5, 0.75}; TF weight in {0, 0.25, 0.5, 1, 2};
5 seeds. K = number of the authors' cell-type niches on the slide (8 / 8 / 9).
Holding K: per seed, the expression resolution (shared by rna and tf) is searched with
TF and spatial weights 0 (= the experiment-0 ratio-0 baseline); then, for each setting,
the spatial resolution is searched to hit K. Baseline for "ARI to baseline": TF weight 0
at the same spatial weight and seed (= RNA + space, experiment 0).

Usage (env spatialleiden, from Derived_Modalities/):
    python 04_multimodal_weight_grid/proto_rna_tf/mm_proto.py data/kuppe_2022 results/kuppe_2022/mm_proto_rna_tf
"""

import importlib.util
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
SPATIAL_W = [0.5, 0.75]
TF_W = [0, 0.25, 0.5, 1.0, 2.0]
SEEDS = [0, 1, 2, 3, 4]
DIRECTED = {"rna": False, "tf": False, "spatial": False}


def build_mdata(data: Path, run: str) -> md.MuData:
    rna = e0.load_slide(data, run)
    e0.preprocess(rna)
    tf = e0.load_tf_layer(data, run)
    e0.preprocess_tf(tf)
    assert list(rna.obs_names) == list(tf.obs_names), "spot order differs between layers"
    keep = ["array_row", "array_col", "celltype_niche", "molecular_niche"]
    r = ad.AnnData(obs=rna.obs[keep].copy(), obsp={"connectivities": rna.obsp["connectivities"]})
    t = ad.AnnData(obs=pd.DataFrame(index=tf.obs_names), obsp={"connectivities": tf.obsp["connectivities"]})
    m = md.MuData({"rna": r, "tf": t})
    m.obsp["spatial_connectivities"] = rna.obsp["spatial_connectivities"]
    m.obs[keep] = rna.obs[keep].to_numpy()
    return m, rna


def cluster(m, res_expr, res_sp, w_tf, w_sp, seed, key="tmp"):
    sl.spatialleiden_multimodal(
        m, resolution={"rna": res_expr, "tf": res_expr, "spatial": res_sp},
        layer_weights={"rna": 1.0, "tf": float(w_tf), "spatial": float(w_sp)},
        directed=DIRECTED, n_iterations=-1, random_state=seed, key_added=key)
    return m.obs[key].astype(str).to_numpy()


def run_job(job):
    data, run, k, seed = job
    t0 = time.time()
    m, rna = build_mdata(data, run)
    ei, ej = e0.edges(rna)
    ct = m.obs.celltype_niche.to_numpy()
    # expression resolution: TF and spatial weight 0 (the experiment-0 ratio-0 baseline)
    r_expr, ok0 = e0.search_resolution_steps(
        lambda r: len(set(cluster(m, r, 1.0, 0, 0, seed))), k)
    rows, labels = [], {}
    for w_sp in SPATIAL_W:
        for w_tf in TF_W:
            r_sp, ok = e0.search_resolution_steps(
                lambda r: len(set(cluster(m, r_expr, r, w_tf, w_sp, seed))), k, start=0.4)
            lab = cluster(m, r_expr, r_sp, w_tf, w_sp, seed)
            labels[(w_sp, w_tf)] = lab
            obs = e0.coherence(lab, ei, ej)
            shuf = e0.shuffled_coherence(lab, ei, ej)
            rows.append({
                "slide": run, "seed": seed, "spatial_weight": w_sp, "tf_weight": w_tf,
                "target_k": k, "n_clusters": len(set(lab)), "reached_k": len(set(lab)) == k,
                "res_expression": r_expr, "res_spatial": r_sp,
                "coherence": obs, "coherence_shuffled": shuf, "coherence_gap": obs - shuf,
                "ARI_to_rna_space": adjusted_rand_score(labels[(w_sp, 0)], lab),
                "ARI_celltype_niche": adjusted_rand_score(ct, lab),
            })
    msg = f"{run} seed {seed}: {len(rows)} runs in {time.time() - t0:.0f}s"
    return run, seed, rows, labels, list(m.obs_names), msg


def main():
    data, out = Path(sys.argv[1]), Path(sys.argv[2])
    out.mkdir(parents=True, exist_ok=True)
    k = {run: int(pd.read_csv(data / "visium_tables" / f"{run}_spots.csv", usecols=["celltype_niche"])
                  .celltype_niche.nunique()) for run in SLIDES}
    jobs = [(data, run, k[run], seed) for run in SLIDES for seed in SEEDS]
    results = []
    t0 = time.time()
    with get_context("spawn").Pool(int(sys.argv[3]) if len(sys.argv) > 3 else 8) as pool:
        for res in pool.imap_unordered(run_job, jobs):
            print(res[5], flush=True)
            results.append(res)
    df = pd.DataFrame([r for res in results for r in res[2]])
    stab = []
    for run in SLIDES:
        per = {res[1]: res[3] for res in results if res[0] == run}
        names = next(res[4] for res in results if res[0] == run)
        lab_tab = {}
        for (w_sp, w_tf) in per[SEEDS[0]]:
            okseeds = [s for s in SEEDS if len(set(per[s][(w_sp, w_tf)])) == k[run]]
            aris = [adjusted_rand_score(per[a][(w_sp, w_tf)], per[b][(w_sp, w_tf)])
                    for a, b in combinations(okseeds, 2)]
            stab.append({"slide": run, "spatial_weight": w_sp, "tf_weight": w_tf,
                         "seed_ARI_mean": np.mean(aris) if aris else np.nan, "n_seeds_reached_k": len(okseeds)})
            for s in SEEDS:
                lab_tab[f"sp{w_sp}_tf{w_tf}_s{s}"] = per[s][(w_sp, w_tf)]
        pd.DataFrame(lab_tab, index=names).to_csv(out / f"labels_{run}.csv.gz")
    df = df.merge(pd.DataFrame(stab), on=["slide", "spatial_weight", "tf_weight"])
    df.to_csv(out / "mm_proto_sweep.csv", index=False)
    print(f"done: {len(df)} runs, {int((~df.reached_k).sum())} missed K, {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
