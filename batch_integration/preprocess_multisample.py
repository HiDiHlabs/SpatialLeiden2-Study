#!/usr/bin/env python

import argparse
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(
        description="Preprocess multiple samples from one experiment"
    )
    parser.add_argument("input", help="Path to anndata file(s).", nargs="+")
    parser.add_argument("output", help="Path to anndata file.", type=Path)
    parser.add_argument(
        "-b",
        "--batch",
        help="Batch key(s). Must be present in obs. Only first will be used for selecting HVGs and should be unique across all batches",
        required=False,
        nargs="*",
        default=None,
    )
    parser.add_argument("--hvg", help="Whether to detect HVGs.", action="store_true")
    parser.add_argument(
        "--seed", help="Seed for randomness.", required=False, type=int, default=42
    )

    args = parser.parse_args()

    import anndata as ad
    from multisample_utils import preprocess_genes, preprocess_scaling

    if len(args.input) == 1:
        adata = ad.read_h5ad(args.input[0])
        batch_key = None
    elif len(args.input) > 1 and args.bacth is not None:
        adata = ad.concat(
            [ad.read_h5ad(file) for file in args.input],
            label="sample_file",
            index_unique="_",
        )
        adata.uns["batch_keys"] = args.batch
        batch_key = adata.uns["batch_keys"][0]
    else:
        raise Exception("If more than 1 file is supported batch must be set")

    preprocess_genes(adata, batch_key=batch_key)
    preprocess_scaling(adata, seed=args.seed)

    adata.write(args.output)


if __name__ == "__main__":
    main()
