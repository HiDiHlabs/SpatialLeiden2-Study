"""Embedding loader shared by exp5_unfixed_k_multimodal.py (same logic as exp4_unfixed_k.py).

src formats:
  npz:<path>:<key>        obs_names "<run>_<barcode>" (exp1 RNA, tf_harmony)
  h5ad:<path>:<obsm key>  read with h5py (files from newer anndata); spots mapped to Visium runs
                          by identical barcode sets per obs "sample", barcodes from obs
                          "original_spot" (the colleague's CCC file)
Returns a DataFrame (spots x dims) indexed "<run>_<barcode>".
"""
from pathlib import Path

import h5py
import numpy as np
import pandas as pd


def _dec(a):
    return np.array([x.decode() if isinstance(x, bytes) else x for x in a], dtype=object)


def _cat(f, col):
    return _dec(f[f"obs/{col}/categories"][:])[f[f"obs/{col}/codes"][:]]


def load(src: str, data: Path, log=print) -> pd.DataFrame:
    kind, path, key = src.split(":", 2)
    if kind == "npz":
        e = np.load(path, allow_pickle=True)
        return pd.DataFrame(e[key], index=e["obs_names"].astype(str))
    with h5py.File(path, "r") as f:
        X = f[f"obsm/{key}"][:]
        smp, bc = _cat(f, "sample"), _cat(f, "original_spot")
    tables = {p.name.removesuffix("_spots.csv"): set(pd.read_csv(p, usecols=["barcode"]).barcode)
              for p in (Path(data) / "visium_tables").glob("*_spots.csv")}
    run_of = {}
    for s in np.unique(smp):
        b = set(bc[smp == s])
        hits = [r for r, t in tables.items() if t == b]
        assert len(hits) == 1, f"{s}: spot tables with the same barcodes: {hits}"
        run_of[s] = hits[0]
    log(f"{Path(path).name}: {len(run_of)} samples mapped to Visium runs")
    return pd.DataFrame(X, index=[f"{run_of[s]}_{b}" for s, b in zip(smp, bc)])
