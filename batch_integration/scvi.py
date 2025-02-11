#!/usr/bin/env python

import argparse
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description="Batch integration using scVI")
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
        "--n_comp",
        help="Number of components in latent space.",
        required=False,
        type=int,
        default=30,
    )
    parser.add_argument(
        "--worker",
        help="Number of workers for data loading.",
        required=False,
        type=int,
        default=7,
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
    import scvi
    import torch
    from scipy.sparse import save_npz

    train_kwargs = (
        dict(accelerator="gpu", devices=1) if torch.cuda.is_available() else {}
    )

    scvi.settings.dl_num_workers = args.worker

    adata = ad.read_h5ad(args.input)
    adata.layers["counts"] = adata.X

    other_batch_keys = (
        adata.uns["batch_keys"][1:] if len(adata.uns["batch_keys"]) > 1 else None
    )
    scvi.model.SCVI.setup_anndata(
        adata,
        layer="counts",
        batch_key=adata.uns["batch_keys"][0],
        categorical_covariate_keys=other_batch_keys,
    )
    model = scvi.model.SCVI(
        adata, n_layers=2, n_latent=args.n_comp, gene_likelihood="nb"
    )
    model.train(enable_progress_bar=False, **train_kwargs)

    adata.obsm["X_scVI"] = model.get_latent_representation()
    preprocess_neighbors(
        adata, obsm_key="X_scVI", n_neighbors=args.n_neighs, seed=args.seed
    )
    save_npz(args.output, adata.obsp["connectivities"])


if __name__ == "__main__":
    main()
