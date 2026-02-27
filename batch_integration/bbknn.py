#!/usr/bin/env python

import argparse
import warnings
from pathlib import Path

warnings.filterwarnings("ignore", category=FutureWarning)


def main():
    parser = argparse.ArgumentParser(description="Batch integration using bbknn")
    parser.add_argument("input", help="Path to anndata file.", type=Path)
    parser.add_argument(
        "output", help="Path to npz file for neighbor graph.", type=Path
    )
    parser.add_argument(
        "--n_neighs",
        help="Number of neighbors to use per batch.",
        required=False,
        type=int,
        default=5,
    )
    parser.add_argument(
        "--seed", help="Seed for randomness.", required=False, type=int, default=42
    )

    args = parser.parse_args()

    import sys

    # prevent importing from cwd
    sys.path.remove(str(Path(__file__).parent.resolve()))

    import anndata as ad
    import scanpy.external as sce
    from scipy.sparse import save_npz

    adata = ad.read_h5ad(args.input)
    sce.pp.bbknn(
        adata,
        batch_key=adata.uns["batch_keys"][0],  # args.batch_key,
        neighbors_within_batch=args.n_neighs,
        pynndescent_random_state=args.seed,
    )
    save_npz(args.output, adata.obsp["connectivities"])


if __name__ == "__main__":
    main()
