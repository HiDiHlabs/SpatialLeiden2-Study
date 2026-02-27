#!/usr/bin/env python

import argparse
import warnings
from pathlib import Path

warnings.filterwarnings("ignore", category=FutureWarning)


def main():
    parser = argparse.ArgumentParser(
        description="Run SpatialLeiden on a batch integrated sample"
    )
    parser.add_argument("input", help="Path to anndata file.", type=Path)
    parser.add_argument("neighbors", help="Path to npz file.", type=Path)
    parser.add_argument("output", help="Path to parquet file.", type=Path)
    # parser.add_argument(
    #     "-n", "--n_cluster", help="Number of clusters", required=True, type=int
    # )
    parser.add_argument(
        "-w",
        "--weight",
        help="Weight ratio for gene expression:spatial layer",
        required=True,
        type=float,
    )
    parser.add_argument(
        "--use_batch",
        help="Whether to use the first batch key for building disjoint spatial graphs",
        action="store_true",
    )
    parser.add_argument(
        "--coord_type",
        help="Coordinate type of the spatial coordinates.",
        choices=["grid", "generic"],
    )
    parser.add_argument(
        "--n_neighs",
        help="Number of neighbors to use.",
        required=False,
        type=int,
        default=6,
    )
    parser.add_argument(
        "--seed", help="Seed for randomness.", required=False, type=int, default=42
    )

    args = parser.parse_args()

    import anndata as ad
    import squidpy as sq
    from scipy.sparse import load_npz
    from spatialleiden import search_resolution

    adata = ad.read_h5ad(args.input)
    adata.obsp["connectivities"] = load_npz(args.neighbors)

    library_key = adata.uns["batch_keys"][0] if args.use_batch else None
    sq.gr.spatial_neighbors(
        adata,
        library_key=library_key,
        coord_type=args.coord_type,
        n_neighs=args.n_neighs,
    )

    directed = False if args.coord_type == "grid" else True
    n_cluster = adata.uns["n_clusters"]

    _ = search_resolution(
        adata,
        n_cluster,
        latent_kwargs={
            "directed": False,
            "random_state": args.seed,
        },
        spatial_kwargs={
            "directed": (False, directed),
            "layer_ratio": args.weight,
            "random_state": args.seed,
        },
    )

    columns = ["label", "spatialleiden"]
    if args.use_batch:
        columns += list(adata.uns["batch_keys"])

    labels = adata.obs[columns]
    labels.loc[:, ["x", "y"]] = adata.obsm["spatial"]
    labels.to_parquet(args.output)


if __name__ == "__main__":
    main()
