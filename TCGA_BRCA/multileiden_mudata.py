"""MuData helpers for SpatialLeiden2-style multimodal MultiLeiden.

Matches HiDiHlabs/SpatialLeiden2-Study VisiumHD-CRC / SpatialGlue pattern:
  MuData({"rna": ..., "HE": ...}) with per-modality connectivities and
  spatial neighbors on mdata, then spatialleiden.spatialleiden_multimodal.
"""

from __future__ import annotations

from typing import Mapping

import mudata as md
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
from anndata import AnnData
from mudata import MuData
from scipy import sparse


def _as_float32_csr(X) -> sparse.csr_matrix:
    if sparse.issparse(X):
        return X.tocsr().astype(np.float32)
    return sparse.csr_matrix(np.asarray(X, dtype=np.float32))


def build_gene_he_mudata(
    *,
    gene_X,
    gene_var_names,
    he_X,
    spatial: np.ndarray,
    obs: pd.DataFrame | None = None,
    obs_names=None,
    gene_obsm: Mapping[str, np.ndarray] | None = None,
    he_obsm: Mapping[str, np.ndarray] | None = None,
    rna_key: str = "rna",
    he_key: str = "HE",
) -> MuData:
    """Build MuData with gene and HE modalities aligned on the same spots."""

    n = spatial.shape[0]
    if obs_names is None:
        obs_names = pd.Index([str(i) for i in range(n)], name="spot")
    if obs is None:
        obs = pd.DataFrame(index=obs_names)
    else:
        obs = obs.copy()
        obs.index = obs_names

    rna = AnnData(
        X=_as_float32_csr(gene_X),
        obs=obs.copy(),
        var=pd.DataFrame(index=pd.Index(list(gene_var_names), name="gene")),
    )
    if gene_obsm:
        for k, v in gene_obsm.items():
            rna.obsm[k] = np.asarray(v)

    he = AnnData(
        X=np.asarray(he_X, dtype=np.float32),
        obs=obs.copy(),
        var=pd.DataFrame(
            index=pd.Index([f"he_{i}" for i in range(np.asarray(he_X).shape[1])], name="he_feature")
        ),
    )
    if he_obsm:
        for k, v in he_obsm.items():
            he.obsm[k] = np.asarray(v)

    mdata = md.MuData({rna_key: rna, he_key: he})
    mdata.obs = obs.copy()
    mdata.obsm["spatial"] = np.asarray(spatial, dtype=np.float64)
    return mdata


def mudata_from_pipeline_anndata(
    adata: AnnData,
    *,
    morph_key: str = "pca",
    rna_key: str = "rna",
    he_key: str = "HE",
) -> MuData:
    """Split a gene+morphology AnnData into SpatialLeiden2 MuData."""
    if morph_key not in adata.obsm and "pca" in adata.obsm:
        morph_key = "pca"
    if morph_key not in adata.obsm:
        raise ValueError(f"Expected morphology in adata.obsm['{morph_key}'] (has {list(adata.obsm)})")
    if "spatial" not in adata.obsm:
        raise ValueError("Expected adata.obsm['spatial']")

    gene_X = adata.layers["counts"] if "counts" in adata.layers else adata.X
    return build_gene_he_mudata(
        gene_X=gene_X,
        gene_var_names=adata.var_names,
        he_X=adata.obsm[morph_key],
        spatial=adata.obsm["spatial"],
        obs=adata.obs,
        obs_names=adata.obs_names,
        gene_obsm={k: adata.obsm[k] for k in ("X_pca", "X_umap") if k in adata.obsm},
        rna_key=rna_key,
        he_key=he_key,
    )


def _obsm_c_contiguous(adata: AnnData, key: str) -> None:
    """Force C-contiguous float64 for UMAP/numba knn (avoids TypingError on float32 A-order)."""
    adata.obsm[key] = np.ascontiguousarray(adata.obsm[key], dtype=np.float64)


def prepare_modality_neighbors(
    mdata: MuData,
    *,
    rna_key: str = "rna",
    he_key: str = "HE",
    n_pcs_rna: int = 50,
    n_pcs_he: int | None = None,
    n_neighs_spatial: int = 10,
    coord_type: str = "generic",
    random_state: int = 42,
    use_distance2connectivity: bool = False,
    log1p_rna: bool = False,
) -> None:
    """PCA + knn per modality; spatial neighbors on mdata (SpatialLeiden2 style).

    Set log1p_rna=True only for raw count matrices. DeepSpot-M atlas X is already
    log1p-CP10k — leave log1p_rna=False.

    Spatial graph: keep squidpy ``spatial_connectivities`` by default (VisiumHD-CRC
    multimodal). ``distance2connectivity`` + directed spatial at equal layer
    weights makes the spatial layer a no-op in multiplex Leiden; only enable it
    for classic two-layer SpatialLeiden (``layer_ratio``) paths.
    """
    import spatialleiden as sl

    rna = mdata.mod[rna_key]
    if "counts" in rna.layers:
        rna.X = rna.layers["counts"]
        log1p_rna = True
    if log1p_rna:
        sc.pp.log1p(rna)
    n_comps = min(n_pcs_rna, max(1, rna.n_vars - 1), max(1, rna.n_obs - 1))
    sc.pp.pca(rna, n_comps=n_comps, random_state=random_state)
    _obsm_c_contiguous(rna, "X_pca")
    # method='gauss' avoids umap fuzzy_simplicial_set (numba TypingError on float32 A-order knn_dists)
    sc.pp.neighbors(rna, use_rep="X_pca", method="gauss", random_state=random_state)

    he = mdata.mod[he_key]
    if n_pcs_he is None:
        # Use raw HE embedding as neighbor space (FM features).
        he.obsm["X_he"] = np.ascontiguousarray(np.asarray(he.X), dtype=np.float64)
        he_rep = "X_he"
    else:
        sc.pp.pca(he, n_comps=min(n_pcs_he, max(1, he.n_vars - 1), max(1, he.n_obs - 1)), random_state=random_state)
        he_rep = "X_pca"
        _obsm_c_contiguous(he, "X_pca")
    sc.pp.neighbors(he, use_rep=he_rep, method="gauss", random_state=random_state)

    if "spatial" not in mdata.obsm:
        raise ValueError("mdata.obsm['spatial'] is required")
    sq.gr.spatial_neighbors(mdata, coord_type=coord_type, n_neighs=n_neighs_spatial)
    if use_distance2connectivity and "spatial_distances" in mdata.obsp:
        mdata.obsp["spatial_connectivities"] = sl.distance2connectivity(mdata.obsp["spatial_distances"])


def run_spatialleiden_multimodal(
    mdata: MuData,
    *,
    resolution: float | Mapping[str, float] = 1.0,
    layer_weights: float | Mapping[str, float] = 1.0,
    neighbors_key: str | Mapping[str, str] = "connectivities",
    spatial_neighbors_key: str = "spatial_connectivities",
    directed: bool | Mapping[str, bool] = False,
    key_added: str = "spatialleiden",
    random_state: int = 42,
    n_iterations: int = -1,
) -> None:
    """Wrapper around spatialleiden.spatialleiden_multimodal."""
    from spatialleiden import spatialleiden_multimodal

    spatialleiden_multimodal(
        mdata,
        resolution=resolution,
        layer_weights=layer_weights,
        neighbors_key=neighbors_key,
        spatial_neighbors_key=spatial_neighbors_key,
        directed=directed,
        key_added=key_added,
        random_state=random_state,
        n_iterations=n_iterations,
    )


def copy_mudata_clusters_to_anndata(mdata: MuData, adata: AnnData, columns: list[str] | None = None) -> None:
    """Copy clustering columns from mdata.obs onto a parallel AnnData."""
    cols = columns or [c for c in mdata.obs.columns if c not in adata.obs.columns or c.endswith("leiden") or "morph" in c or "gene" in c or "spatial" in c]
    for c in mdata.obs.columns:
        if columns is not None and c not in columns:
            continue
        if c in mdata.obs:
            adata.obs[c] = mdata.obs[c].values
