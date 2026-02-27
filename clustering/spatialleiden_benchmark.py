#!/usr/bin/env python

import argparse
import warnings
from pathlib import Path

warnings.filterwarnings("ignore", category=FutureWarning)


def main():
    parser = argparse.ArgumentParser(
        description="Run SpatialLeiden for n iterations (no output)"
    )
    parser.add_argument("input", help="Path to anndata file.", type=Path)
    parser.add_argument(
        "--weight",
        help="Weight ratio for gene expression:spatial layer",
        required=True,
        type=float,
    )
    parser.add_argument(
        "--resolution",
        help="Resolution parameter for the partition.",
        required=False,
        type=float,
        default=1.0,
    )
    parser.add_argument(
        "--n_iter", help="Number of iterations", required=False, type=int, default=5
    )
    parser.add_argument(
        "--directed", help="Whether spatial graph is directed", action="store_true"
    )
    parser.add_argument("--seed", help="Seed for randomness.", required=True, type=int)

    args = parser.parse_args()

    import anndata as ad
    from spatialleiden import spatialleiden

    adata = ad.read_h5ad(args.input)

    spatialleiden(
        adata,
        resolution=(args.resolution, args.resolution),
        directed=(False, args.directed),
        layer_ratio=args.weight,
        n_iterations=args.n_iter,
        random_state=args.seed,
    )


if __name__ == "__main__":
    main()
