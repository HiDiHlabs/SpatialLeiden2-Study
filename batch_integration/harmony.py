#!/usr/bin/env python

import argparse
import warnings
from pathlib import Path

warnings.filterwarnings("ignore", category=FutureWarning)


def main():
    parser = argparse.ArgumentParser(description="Batch integration using harmony")
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

    import anndata as ad
    import harmonypy as hm
    from multisample_utils import preprocess_neighbors
    from scipy.sparse import save_npz

    adata = ad.read_h5ad(args.input)

    harmony_result = hm.run_harmony(
        adata.obsm["X_pca"], adata.obs, adata.uns["batch_keys"], random_state=args.seed
    )
    adata.obsm["X_pca_harmony"] = harmony_result.Z_corr

    preprocess_neighbors(
        adata, obsm_key="X_pca_harmony", n_neighbors=args.n_neighs, seed=args.seed
    )
    save_npz(args.output, adata.obsp["connectivities"])


if __name__ == "__main__":
    main()
