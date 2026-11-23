"""Helpers to store FM embeddings without cloning gene counts or WSI images."""

from __future__ import annotations

import numpy as np
import scanpy as sc


def squeeze_2d(x, n_obs: int) -> np.ndarray:
    x = np.asarray(x)
    x = np.squeeze(x)
    if x.ndim == 1:
        x = np.reshape(x, (n_obs, -1))
    elif x.ndim > 2:
        x = np.reshape(x, (x.shape[0], -1))
    if x.shape[0] != n_obs:
        raise ValueError(f"Embedding has {x.shape[0]} rows, expected {n_obs}")
    return np.asarray(x, dtype=np.float32)


def drop_spatial_images(adata) -> None:
    spatial = adata.uns.get("spatial")
    if not isinstance(spatial, dict):
        return
    for spec in spatial.values():
        if isinstance(spec, dict):
            spec.pop("images", None)


def write_thin_obsm(src, features, key: str, path: str) -> None:
    out = sc.AnnData(obs=src.obs.copy())
    out.obsm[key] = squeeze_2d(features, src.n_obs)
    out.write(path)
