#!/usr/bin/env python

import argparse
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description="Batch integration using scanorama")
    parser.add_argument("input", help="Path to anndata file.", type=Path)
    parser.add_argument(
        "output", help="Path to npz file for neighbor graph.", type=Path
    )
    parser.add_argument(
        "--n_neighs",
        help="Number of neighbors to use per batch.",
        required=False,
        type=int,
        default=15,
    )
    parser.add_argument(
        "--seed", help="Seed for randomness.", required=False, type=int, default=42
    )

    args = parser.parse_args()

    import sys

    from multisample_utils import preprocess_neighbors

    # prevent importing from cwd
    sys.path.remove(str(Path(__file__).parent.resolve()))

    import anndata as ad
    import scanpy.external as sce
    from scipy.sparse import save_npz

    adata = ad.read_h5ad(args.input)
    sce.pp.scanorama_integrate(
        adata,
        key=adata.uns["batch_keys"][0],
        basis="X_pca",
        adjusted_basis="X_pca_scanorama",
        verbose=False,
    )
    preprocess_neighbors(
        adata, obsm_key="X_pca_scanorama", n_neighbors=args.n_neighs, seed=args.seed
    )
    save_npz(args.output, adata.obsp["connectivities"])


if __name__ == "__main__":
    main()
