#!/usr/bin/env python3
"""Run SpatialLeiden2 multimodal MultiLeiden on a MuData (rna + HE + spatial).

Input: .h5mu from retrieve_deepspotm_tcga_brca.py (or SpatialLeiden2-style MuData).
Writes clustering labels onto mdata.obs.

By default runs one weighted multimodal clustering (``--layer-weights`` /
``--key-added``). Pass ``--all-submethods`` to also write the six unweighted
SpatialLeiden columns written to obs:
  gene, gene_location, morphology, morphology_location, morphology_gene, gene_morph_loc

Example
-------
  python TCGA_BRCA/spatialleiden_mudata_run.py \\
    --input ./data/TCGA_BRCA/mudata/SLIDE.h5mu \\
    --output ./data/TCGA_BRCA/spatialleiden/SLIDE.h5mu \\
    --all-submethods --n-clusters 8
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import mudata as md
import numpy as np
import scanpy as sc
import spatialleiden as sl

sys.path.insert(0, str(Path(__file__).resolve().parent))
from multileiden_mudata import (  # noqa: E402
    prepare_modality_neighbors,
    run_spatialleiden_multimodal,
)

# Same order as plot_spatial_clusters.SLEIDEN_SPATIAL_COL_ORDER / aggregate scores.
UNWEIGHTED_SUBMETHODS = [
    "gene",
    "gene_location",
    "morphology",
    "morphology_location",
    "morphology_gene",
    "gene_morph_loc",
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--rna-key", default="rna")
    p.add_argument("--he-key", default="HE")
    p.add_argument("--n-clusters", type=int, default=None, help="If set, grid-search a shared resolution")
    p.add_argument("--resolution", type=float, default=0.55)
    p.add_argument("--n-neighs-spatial", type=int, default=8)
    p.add_argument("--coord-type", choices=["generic", "grid"], default="grid")
    p.add_argument("--layer-weights", default="rna:1.5,HE:0.5,spatial:2")
    p.add_argument("--random-state", type=int, default=42)
    p.add_argument("--key-added", default="spatialleiden")
    p.add_argument(
        "--all-submethods",
        action="store_true",
        help="Write all six unweighted Multileiden submethod columns (see UNWEIGHTED_SUBMETHODS)",
    )
    return p.parse_args()


def _parse_weights(text: str) -> dict[str, float]:
    out: dict[str, float] = {}
    for part in text.split(","):
        k, v = part.split(":")
        out[k.strip()] = float(v)
    return out


def _neighbors_key_for_weights(weights: dict[str, float], rna_key: str, he_key: str) -> dict[str, str]:
    nk: dict[str, str] = {}
    if rna_key in weights:
        nk[rna_key] = "connectivities"
    if he_key in weights:
        nk[he_key] = "connectivities"
    return nk


def _search_shared_resolution(mdata, n_clusters: int, weights: dict, args) -> float:
    """Coarse search for a single resolution yielding ~n_clusters."""
    best_res = args.resolution
    best_diff = 1e9
    nk = _neighbors_key_for_weights(weights, args.rna_key, args.he_key)
    for res in np.linspace(0.1, 2.0, 20):
        run_spatialleiden_multimodal(
            mdata,
            resolution=float(res),
            layer_weights=weights,
            neighbors_key=nk or "connectivities",
            directed=False,
            key_added="_res_search",
            random_state=args.random_state,
            n_iterations=5,
        )
        n = int(mdata.obs["_res_search"].nunique())
        diff = abs(n - n_clusters)
        if diff < best_diff:
            best_diff = diff
            best_res = float(res)
        if diff == 0:
            break
    if "_res_search" in mdata.obs:
        del mdata.obs["_res_search"]
    return best_res


def _run_leiden_one_graph(conn, *, resolution: float, random_state: int) -> np.ndarray:
    """Leiden on a single connectivity matrix via multiplex_leiden (one layer)."""
    return np.asarray(
        sl.multiplex_leiden(
            conn,
            directed=False,
            layer_weights=1.0,
            partition_kwargs=[{"resolution_parameter": float(resolution)}],
            random_state=random_state,
        )
    )


def _run_multiplex(
    *conns,
    resolutions: tuple[float, ...],
    directed: tuple[bool, ...],
    random_state: int,
    layer_weights: float | tuple[float, ...] = 1.0,
) -> np.ndarray:
    return np.asarray(
        sl.multiplex_leiden(
            *conns,
            directed=directed,
            layer_weights=layer_weights,
            partition_kwargs=[{"resolution_parameter": float(r)} for r in resolutions],
            random_state=random_state,
        )
    )


def run_all_unweighted_submethods(mdata, *, rna_key: str, he_key: str, resolution: float, random_state: int) -> None:
    """Fill UNWEIGHTED_SUBMETHODS on mdata.obs using prepared neighbor graphs.

    Spatial layer matches SpatialLeiden2 VisiumHD multimodal: squidpy connectivities
    (no distance2connectivity) and directed=False. Using distance2connectivity with
    directed spatial at equal weights leaves gene_location / gene_morph_loc identical
    to the non-spatial counterparts.
    """
    rna = mdata.mod[rna_key]
    he = mdata.mod[he_key]
    rna_conn = rna.obsp["connectivities"]
    he_conn = he.obsp["connectivities"]
    sp_conn = mdata.obsp["spatial_connectivities"]
    res = float(resolution)
    # VisiumHD-CRC multimodal uses directed=False for every layer (incl. spatial).
    undirected = False

    mdata.obs["gene"] = _run_leiden_one_graph(rna_conn, resolution=res, random_state=random_state)
    mdata.obs["morphology"] = _run_leiden_one_graph(he_conn, resolution=res, random_state=random_state)

    mdata.obs["gene_location"] = _run_multiplex(
        rna_conn,
        sp_conn,
        resolutions=(res, res),
        directed=(undirected, undirected),
        random_state=random_state,
    )
    mdata.obs["morphology_location"] = _run_multiplex(
        he_conn,
        sp_conn,
        resolutions=(res, res),
        directed=(undirected, undirected),
        random_state=random_state,
    )
    mdata.obs["morphology_gene"] = _run_multiplex(
        he_conn,
        rna_conn,
        resolutions=(res, res),
        directed=(undirected, undirected),
        random_state=random_state,
    )
    mdata.obs["gene_morph_loc"] = _run_multiplex(
        rna_conn,
        he_conn,
        sp_conn,
        resolutions=(res, res, res),
        directed=(undirected, undirected, undirected),
        random_state=random_state,
    )

    for col in UNWEIGHTED_SUBMETHODS:
        mdata.obs[col] = mdata.obs[col].astype("category")
        print(f"  {col}: n_clusters={mdata.obs[col].nunique()}")


def main() -> None:
    args = parse_args()
    mdata = md.read(args.input)
    if args.rna_key not in mdata.mod or args.he_key not in mdata.mod:
        raise SystemExit(f"Need modalities {args.rna_key!r} and {args.he_key!r}; has {list(mdata.mod)}")

    prepare_modality_neighbors(
        mdata,
        rna_key=args.rna_key,
        he_key=args.he_key,
        n_neighs_spatial=args.n_neighs_spatial,
        coord_type=args.coord_type,
        random_state=args.random_state,
    )
    weights = _parse_weights(args.layer_weights)
    resolution = args.resolution
    if args.n_clusters is not None:
        resolution = _search_shared_resolution(mdata, args.n_clusters, weights, args)
        print(f"Selected resolution={resolution:.3f} for ~{args.n_clusters} clusters")

    if args.all_submethods:
        print(f"Running unweighted Multileiden submethods at resolution={resolution:.3f}")
        run_all_unweighted_submethods(
            mdata,
            rna_key=args.rna_key,
            he_key=args.he_key,
            resolution=resolution,
            random_state=args.random_state,
        )
        mdata.obs["multileiden"] = mdata.obs["gene_morph_loc"]
        # Keep default multimodal key as an alias of gene_morph_loc unless weights differ.
        if args.key_added not in mdata.obs:
            mdata.obs[args.key_added] = mdata.obs["gene_morph_loc"]
    else:
        nk = _neighbors_key_for_weights(weights, args.rna_key, args.he_key)
        run_spatialleiden_multimodal(
            mdata,
            resolution=resolution,
            layer_weights=weights,
            neighbors_key=nk or "connectivities",
            directed=False,
            key_added=args.key_added,
            random_state=args.random_state,
        )
        mdata.obs["gene_morph_loc"] = mdata.obs[args.key_added].astype("category")
        mdata.obs["multileiden"] = mdata.obs["gene_morph_loc"]

    # Cache HE UMAP for plotting (scanpy only; no squidpy).
    he = mdata.mod[args.he_key]
    if "X_umap" not in he.obsm:
        if "X_he" not in he.obsm:
            he.obsm["X_he"] = np.ascontiguousarray(np.asarray(he.X), dtype=np.float64)
        else:
            he.obsm["X_he"] = np.ascontiguousarray(he.obsm["X_he"], dtype=np.float64)
        if "connectivities" not in he.obsp:
            sc.pp.neighbors(he, use_rep="X_he", method="gauss", random_state=args.random_state)
        sc.tl.umap(he, random_state=args.random_state)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    mdata.write(args.output)
    cols = [c for c in UNWEIGHTED_SUBMETHODS if c in mdata.obs.columns]
    nmsg = ", ".join(f"{c}={mdata.obs[c].nunique()}" for c in cols) if cols else f"{args.key_added}={mdata.obs[args.key_added].nunique()}"
    print(f"Wrote {args.output}; {nmsg}")


if __name__ == "__main__":
    main()
