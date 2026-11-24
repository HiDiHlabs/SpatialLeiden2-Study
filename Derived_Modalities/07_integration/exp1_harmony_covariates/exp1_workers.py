"""Worker functions for exp1_covariates.py, run in spawned processes (one BLAS thread each).

leiden_k: step-halving search for a resolution giving exactly K clusters (exp0_sweep's
search_resolution_steps, plain spatialleiden.leiden, undirected) on one variant's 15-NN graph,
plus the silhouette of the labels on that variant's embedding (stratified subsample).
umap: scanpy UMAP of one variant from its precomputed neighbour graph.
"""
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "03_single_modality_clustering" / "exp0_layer_ratio_sweep"))

_G = {}


def init(graphs, n_obs, sil_idx, sil_emb):
    _G.update(graphs=graphs, n_obs=n_obs, sil_idx=sil_idx, sil_emb=sil_emb)


def leiden_k(job):
    import spatialleiden as sl
    from exp0_sweep import search_resolution_steps
    from sklearn.metrics import silhouette_samples

    v, k, seed = job
    with threadpool_limits(1):
        a = ad.AnnData(obs=pd.DataFrame(index=np.arange(_G["n_obs"]).astype(str)))
        a.obsp["c"] = _G["graphs"][v]["connectivities"]

        def n_at(r):
            sl.leiden(a, resolution=r, directed=False, random_state=seed, key_added="_t", neighbors_key="c")
            return a.obs["_t"].nunique()

        r, ok = search_resolution_steps(n_at, k)
        n_at(r)
        lab = a.obs["_t"].astype(str).to_numpy()
        sil = float(silhouette_samples(_G["sil_emb"][v], lab[_G["sil_idx"]]).mean())
    return v, k, r, ok, lab, sil


def umap(job):
    import scanpy as sc

    v, seed = job
    g = _G["graphs"][v]
    with threadpool_limits(1):
        a = ad.AnnData(obs=pd.DataFrame(index=np.arange(_G["n_obs"]).astype(str)))
        a.obsm["X_v"] = g["emb"]
        a.obsp[f"{v}_connectivities"] = g["connectivities"]
        a.obsp[f"{v}_distances"] = g["distances"]
        a.uns[v] = g["uns"]
        sc.tl.umap(a, neighbors_key=v, random_state=seed)
    return v, a.obsm["X_umap"]
