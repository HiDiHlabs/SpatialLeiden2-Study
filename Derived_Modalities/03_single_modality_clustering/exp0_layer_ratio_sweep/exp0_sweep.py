"""Experiment 0: SpatialLeiden layer_ratio sweep on single Kuppe 2022 Visium slides.

One slide at a time, RNA only, `spatialleiden.spatialleiden()` (two layers: the
expression kNN graph and the spatial grid graph). Varies `layer_ratio` and the Leiden
seed; holds the number of clusters fixed. See README.md in this folder for the design.

Steps per slide
  1. Build AnnData: raw counts from the Zenodo h5ad (.raw), joined by barcode to the
     authors' per-spot table (true coordinates, array grid position, published niches).
  2. Preprocess as in the SpatialLeiden 2 preprint (Methods): genes in >=10 spots,
     2,000 HVGs (seurat_v3, counts), normalize_total, log1p, scale without centring,
     PCA 30, 15-NN graph. Spatial graph: squidpy grid, 6 neighbours. Both graphs are
     computed ONCE per slide, so the latent k never changes within the sweep.
  3. Target cluster number K: from the ratio-0 run (plain Leiden at resolution 1.0,
     seed 0) of the reference slide; the same K is used on every slide (--k overrides).
  4. Per seed: layer_ratio 0 -> search the latent resolution giving K with
     sl.spatialleiden(layer_ratio=0) itself, so the baseline uses the same optimiser as
     every other ratio. Sanity check (finding, 23 Nov 2026): sl.leiden and
     sl.multiplex_leiden on the latent layer alone agree exactly (ARI 1), but adding the
     spatial layer at weight 0 changes the optimiser's search path, so
     sl.spatialleiden(layer_ratio=0) lands in a different local optimum of the same
     objective (ARI ~0.8 to plain Leiden). We record that ARI next to the seed-to-seed
     ARI of plain Leiden, which is the noise it should be compared with.
     layer_ratio > 0 -> latent resolution fixed at that seed's ratio-0 value, spatial
     resolution searched to give K (search_resolution_spatial).
  5. Record per (slide, layer_ratio, seed): cluster count, resolutions, search outcome,
     spatial coherence (share of neighbouring spot pairs with the same label) minus the
     same after shuffling labels within the slide (20 shuffles), ARI to the ratio-0
     partition of the same seed, ARI/NMI to the published niches (context only), cluster
     sizes, runtime. Seed stability is computed afterwards (pairwise ARI among seeds).

Usage (from Derived_Modalities/):
    python 03_single_modality_clustering/exp0_layer_ratio_sweep/exp0_sweep.py \
        data/kuppe_2022 results/kuppe_2022/exp0 [--probe] [--k K]
"""

import argparse
import os

for _v in (
    "OMP_NUM_THREADS",
    "MKL_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "NUMBA_NUM_THREADS",
):
    os.environ.setdefault(
        _v, "1"
    )  # one thread per job; parallelism comes from the job pool
# The imports below come after the thread limits on purpose (numpy reads them at import).
import json  # noqa: E402
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
from sklearn.metrics import (  # noqa: E402
    adjusted_rand_score,
    normalized_mutual_info_score,
)
from spatialleiden._resolution_search import _search_resolution  # noqa: E402

SLIDES = {  # run id -> (h5ad name, label)
    "ACH0013": ("Visium_GT_IZ_P9_rep2.h5ad", "P9 IZ rep2 (ACH0013)"),
    "ACH0012": ("Visium_GT_IZ_P9.h5ad", "P9 IZ (ACH0012)"),
    "10X001": ("Visium_control_P1.h5ad", "P1 control (10X001)"),
}
REFERENCE = "ACH0013"
RATIOS = [0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5]
SEEDS = [0, 1, 2, 3, 4]
N_SEARCH = 30  # search iterations (package default 15)
N_SHUFFLE = 20
GRAPH_SEED = 0  # PCA / neighbour graph, fixed for the whole sweep
DIRECTED = (False, False)  # as in the preprint code: symmetric kNN graph, regular grid


def load_slide(data: Path, run: str) -> ad.AnnData:
    h5 = ad.read_h5ad(data / "processed1" / SLIDES[run][0])
    spots = pd.read_csv(data / "visium_tables" / f"{run}_spots.csv").set_index(
        "barcode"
    )
    assert set(spots.index) == set(h5.obs_names), f"{run}: barcodes differ"
    spots = spots.loc[h5.obs_names]  # explicit join by barcode
    a = ad.AnnData(
        X=h5.raw.X.copy(),
        obs=pd.DataFrame(index=h5.obs_names),
        var=pd.DataFrame(index=h5.raw.var_names),
    )
    assert np.allclose(
        np.asarray(a.X.sum(1)).ravel(), spots.nCount_Spatial.to_numpy()
    ), "counts mismatch"
    a.layers["counts"] = a.X.copy()
    a.obsm["spatial"] = spots[["pxl_col_fullres", "pxl_row_fullres"]].to_numpy(float)
    for c in ["array_row", "array_col", "celltype_niche", "molecular_niche"]:
        a.obs[c] = spots[c].to_numpy()
    return a


def preprocess(a: ad.AnnData) -> None:
    sc.pp.filter_genes(a, min_cells=10)
    sc.pp.highly_variable_genes(a, n_top_genes=2000, flavor="seurat_v3", layer="counts")
    sc.pp.normalize_total(a)
    sc.pp.log1p(a)
    sc.pp.scale(a, zero_center=False)
    sc.pp.pca(a, n_comps=30, mask_var="highly_variable", random_state=GRAPH_SEED)
    sc.pp.neighbors(a, n_neighbors=15, n_pcs=30, random_state=GRAPH_SEED)
    sq.gr.spatial_neighbors(a, coord_type="grid", n_neighs=6)


def edges(a: ad.AnnData):
    i, j = a.obsp["spatial_connectivities"].nonzero()
    return i, j


def coherence(labels, i, j) -> float:
    return float((labels[i] == labels[j]).mean())


def run_seed(job):
    """One (slide, seed): every layer_ratio. Runs in its own process."""
    data, run, k, seed = job
    t0 = time.time()
    a = load_slide(data, run)
    preprocess(a)  # deterministic: identical graphs in every job of this slide
    ei, ej = edges(a)
    ct = a.obs.celltype_niche.to_numpy()
    mol = a.obs.molecular_niche.to_numpy()
    rng = np.random.default_rng(seed)  # label shuffles
    deg = np.bincount(np.bincount(ei, minlength=a.n_obs)).tolist()
    msgs = [
        f"{run} seed {seed}: {a.n_obs} spots, {a.n_vars} genes, spatial degree {deg}, "
        f"preprocessing {time.time() - t0:.0f}s"
    ]
    rows, labels = [], {}

    def n_ratio0(res):
        sl.spatialleiden(
            a,
            resolution=(res, 1.0),
            layer_ratio=0,
            directed=DIRECTED,
            n_iterations=-1,
            random_state=seed,
            key_added="tmp",
        )
        return a.obs["tmp"].nunique()

    # ratio 0: search the latent resolution for K with spatialleiden(layer_ratio=0) itself
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        r_lat = _search_resolution(n_ratio0, k, start=1.0, step=0.1, n_iter=N_SEARCH)
    failed0 = any("Correct resolution not found" in str(x.message) for x in w)
    # sanity check: plain Leiden and latent-only multiplex at the same resolution and seed
    sl.leiden(
        a, resolution=r_lat, directed=False, random_state=seed, key_added="leiden_plain"
    )
    plain = a.obs["leiden_plain"].astype(str).to_numpy()
    latent_only = sl.multiplex_leiden(
        a.obsp["connectivities"],
        directed=False,
        layer_weights=1.0,
        partition_kwargs={"resolution_parameter": r_lat},
        random_state=seed,
    )
    sanity = None
    for ratio in RATIOS:
        t1 = time.time()
        key = f"r{ratio}_s{seed}"
        failed = failed0
        r_sp = np.nan
        if ratio == 0:
            sl.spatialleiden(
                a,
                resolution=(r_lat, 1.0),
                layer_ratio=0,
                directed=DIRECTED,
                n_iterations=-1,
                random_state=seed,
                key_added=key,
            )
            sanity = {
                "slide": run,
                "seed": seed,
                "res_latent": r_lat,
                "K_plain_leiden": int(len(set(plain))),
                "K_spatialleiden_ratio0": int(a.obs[key].nunique()),
                "ARI_leiden_vs_latent_only_multiplex": adjusted_rand_score(
                    plain, latent_only
                ),
                "ARI_leiden_vs_spatialleiden_ratio0": adjusted_rand_score(
                    plain, a.obs[key]
                ),
            }
        else:
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                r_sp = sl.search_resolution_spatial(
                    a,
                    k,
                    n_iter=N_SEARCH,
                    resolution=(r_lat, 1.0),
                    layer_ratio=ratio,
                    directed=DIRECTED,
                    n_iterations=-1,
                    random_state=seed,
                    key_added=key,
                )
            failed = any("Correct resolution not found" in str(x.message) for x in w)
            # re-run at the found resolution so the stored labels are exactly that run
            sl.spatialleiden(
                a,
                resolution=(r_lat, r_sp),
                layer_ratio=ratio,
                directed=DIRECTED,
                n_iterations=-1,
                random_state=seed,
                key_added=key,
            )
        lab = a.obs[key].astype(str).to_numpy()
        labels[key] = lab
        sizes = np.sort(np.bincount(pd.factorize(lab)[0]))[::-1]
        obs = coherence(lab, ei, ej)
        shuf = float(
            np.mean([coherence(rng.permutation(lab), ei, ej) for _ in range(N_SHUFFLE)])
        )
        base = labels[f"r0_s{seed}"]
        rows.append(
            {
                "slide": run,
                "label": SLIDES[run][1],
                "layer_ratio": ratio,
                "seed": seed,
                "n_clusters": len(sizes),
                "target_k": k,
                "search_failed": failed,
                "res_latent": r_lat,
                "res_spatial": r_sp,
                "coherence": obs,
                "coherence_shuffled": shuf,
                "coherence_gap": obs - shuf,
                "ARI_to_ratio0": adjusted_rand_score(base, lab),
                "ARI_celltype_niche": adjusted_rand_score(ct, lab),
                "NMI_celltype_niche": normalized_mutual_info_score(ct, lab),
                "ARI_molecular_niche": adjusted_rand_score(mol, lab),
                "NMI_molecular_niche": normalized_mutual_info_score(mol, lab),
                "largest_cluster_frac": sizes[0] / a.n_obs,
                "smallest_cluster_frac": sizes[-1] / a.n_obs,
                "n_clusters_below_1pct": int((sizes < 0.01 * a.n_obs).sum()),
                "cluster_sizes": ";".join(map(str, sizes)),
                "runtime_s": time.time() - t1,
            }
        )
        msgs.append(
            f"  {run} ratio {ratio:<4} seed {seed}: K={len(sizes)} res=({r_lat:.3f}, {r_sp:.3f}) "
            f"gap={obs - shuf:.3f} ARI0={rows[-1]['ARI_to_ratio0']:.3f} failed={failed}"
        )
    meta = pd.DataFrame(
        {
            "x": a.obsm["spatial"][:, 0],
            "y": a.obsm["spatial"][:, 1],
            "array_row": a.obs.array_row.to_numpy(),
            "array_col": a.obs.array_col.to_numpy(),
            "celltype_niche": ct,
            "molecular_niche": mol,
        },
        index=a.obs_names,
    )
    return run, seed, rows, labels, sanity, plain, meta, msgs


def assemble_slide(run, results, out_dir):
    """Combine the per-seed jobs of one slide: seed stability, sanity noise level, label table."""
    results = sorted(results, key=lambda r: r[1])
    rows = [row for r in results for row in r[2]]
    labels = {k: v for r in results for k, v in r[3].items()}
    plain = {r[1]: r[5] for r in results}
    df = pd.DataFrame(rows)
    stab = []
    for ratio in RATIOS:
        aris = [
            adjusted_rand_score(labels[f"r{ratio}_s{s1}"], labels[f"r{ratio}_s{s2}"])
            for s1, s2 in combinations(SEEDS, 2)
        ]
        stab.append(
            {
                "slide": run,
                "layer_ratio": ratio,
                "seed_ARI_mean": np.mean(aris),
                "seed_ARI_min": np.min(aris),
                "seed_ARI_max": np.max(aris),
            }
        )
    df = df.merge(pd.DataFrame(stab), on=["slide", "layer_ratio"])
    plain_aris = [
        adjusted_rand_score(plain[s1], plain[s2]) for s1, s2 in combinations(SEEDS, 2)
    ]
    sanity = pd.DataFrame([r[4] for r in results])
    sanity["plain_leiden_seed_ARI_mean"] = np.mean(plain_aris)
    sanity["plain_leiden_seed_ARI_min"] = np.min(plain_aris)
    labs = pd.concat(
        [pd.DataFrame(labels, index=results[0][6].index), results[0][6]], axis=1
    )
    labs.to_csv(out_dir / f"labels_{run}.csv.gz")
    return df, sanity


def main():
    p = argparse.ArgumentParser()
    p.add_argument("data", type=Path)
    p.add_argument("out", type=Path)
    p.add_argument(
        "--probe", action="store_true", help="only report K at resolution 1.0 per slide"
    )
    p.add_argument("--k", type=int, default=None)
    p.add_argument("--jobs", type=int, default=15)
    p.add_argument(
        "--k-from",
        choices=["celltype_niche", "molecular_niche"],
        default=None,
        help="set K per slide to the number of the authors' niches of this kind on that slide",
    )
    args = p.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    logf = open(args.out / "exp0_log.txt", "a")

    def log(msg):
        print(msg, flush=True)
        logf.write(msg + "\n")
        logf.flush()

    log(
        f"spatialleiden {sl.__version__}, scanpy {sc.__version__}, squidpy {sq.__version__}"
    )
    k_per_slide = None
    if args.k_from:
        k_per_slide = {
            run: int(
                pd.read_csv(
                    args.data / "visium_tables" / f"{run}_spots.csv",
                    usecols=[args.k_from],
                )[args.k_from].nunique()
            )
            for run in SLIDES
        }
        log(
            f"K per slide = number of authors' {args.k_from} labels on the slide: {k_per_slide}"
        )
    elif args.probe or args.k is None:
        ks = {}
        for run in SLIDES:
            a = load_slide(args.data, run)
            preprocess(a)
            sl.leiden(
                a, resolution=1.0, directed=False, random_state=0, key_added="probe"
            )
            ks[run] = int(a.obs["probe"].nunique())
            log(
                f"probe {run}: Leiden resolution 1.0, seed 0 -> {ks[run]} clusters; "
                f"published cell-type niches on slide: {a.obs.celltype_niche.nunique()}"
            )
        (args.out / "probe.json").write_text(json.dumps(ks, indent=1))
        if args.probe:
            return
        args.k = ks[REFERENCE]
    if k_per_slide is None:
        log(
            f"target K = {args.k} (from the ratio-0 run of {REFERENCE} unless overridden)"
        )
        k_per_slide = {run: args.k for run in SLIDES}

    # one process per (slide, seed); Leiden is single-threaded, so this is ~15x faster
    jobs = [
        (args.data, run, k_per_slide[run], seed) for run in SLIDES for seed in SEEDS
    ]
    t0 = time.time()
    per_slide = {run: [] for run in SLIDES}
    with get_context("spawn").Pool(min(args.jobs, len(jobs))) as pool:
        for res in pool.imap_unordered(run_seed, jobs):
            per_slide[res[0]].append(res)
            for m in res[7]:
                log(m)
            log(f"finished {res[0]} seed {res[1]} ({time.time() - t0:.0f}s)")
    all_df, all_san = [], []
    for run in SLIDES:
        df, san = assemble_slide(run, per_slide[run], args.out)
        all_df.append(df)
        all_san.append(san)
    pd.concat(all_df).to_csv(args.out / "exp0_sweep.csv", index=False)
    pd.concat(all_san).to_csv(args.out / "exp0_sanity_ratio0.csv", index=False)
    log(f"done in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
