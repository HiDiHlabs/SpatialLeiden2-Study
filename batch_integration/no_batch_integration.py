#!/usr/bin/env python

import argparse
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description="No batch integration, default kNN")
    parser.add_argument("input", help="Path to anndata file.", type=Path)
    parser.add_argument(
        "output", help="Path to npz file for neighbor graph.", type=Path
    )
    parser.add_argument(
        "--n_neighs",
        help="Number of neighbors.",
        required=False,
        type=int,
        default=15,
    )
    parser.add_argument(
        "--seed", help="Seed for randomness.", required=False, type=int, default=42
    )

    args = parser.parse_args()

    import anndata as ad
    from multisample_utils import preprocess_neighbors
    from scipy.sparse import save_npz

    adata = ad.read_h5ad(args.input)
    preprocess_neighbors(adata, "X_pca", n_neighbors=args.n_neighs, seed=args.seed)
    save_npz(args.output, adata.obsp["connectivities"])


if __name__ == "__main__":
    main()
