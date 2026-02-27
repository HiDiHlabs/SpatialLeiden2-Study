#!/usr/bin/env python

import argparse
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description="Run SpatialGlue")
    parser.add_argument("input", help="Path to mudata file.", type=Path)
    parser.add_argument("output", help="Path to output parquet-file.", type=Path)
    parser.add_argument(
        "--seed", help="Seed for randomness.", required=False, type=int, default=42
    )

    args = parser.parse_args()

    import mudata as md
    import scanpy as sc
    import torch
    from SpatialGlue.preprocess import (
        clr_normalize_each_cell,
        construct_neighbor_graph,
        fix_seed,
        pca,
    )
    from SpatialGlue.SpatialGlue_pyG import Train_SpatialGlue
    from SpatialGlue.utils import clustering

    fix_seed(args.seed)

    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    mdata = md.read_h5mu(args.input)

    mdata.mod["RNA"].obsm["spatial"] = mdata.obsm["spatial"]
    mdata.mod["ADT"].obsm["spatial"] = mdata.obsm["spatial"]

    adata1 = mdata.mod["RNA"]
    adata2 = mdata.mod["ADT"]

    # process according to tutorial
    # https://github.com/JinmiaoChenLab/SpatialGlue_tutorials/blob/main/docs/source/Tutorial%205_data%20integration%20for%20triplet%20omics%20data%20(Simulation).ipynb
    n_protein = adata2.n_vars

    # RNA
    sc.pp.highly_variable_genes(adata1, flavor="seurat_v3", n_top_genes=3_000)
    sc.pp.normalize_total(adata1, target_sum=1e4)
    sc.pp.log1p(adata1)
    adata1.obsm["feat"] = pca(adata1, n_comps=n_protein)

    # Protein
    adata2 = clr_normalize_each_cell(adata2)
    adata2.obsm["feat"] = pca(adata2, n_comps=n_protein)

    data = construct_neighbor_graph(adata1, adata2)
    model = Train_SpatialGlue(data, device=device, random_seed=args.seed, epochs=200)

    # train model
    output = model.train()

    adata = adata1.copy()
    adata.obsm["emb_latent_omics1"] = output["emb_latent_omics1"].copy()
    adata.obsm["emb_latent_omics2"] = output["emb_latent_omics2"].copy()
    adata.obsm["SpatialGlue"] = output["SpatialGlue"].copy()
    adata.obsm["alpha"] = output["alpha"]
    adata.obsm["alpha_omics1"] = output["alpha_omics1"]
    adata.obsm["alpha_omics2"] = output["alpha_omics2"]

    tool = "mclust"  # mclust, leiden, and louvain
    clustering(
        adata, key="SpatialGlue", add_key="SpatialGlue", n_clusters=5, method=tool
    )

    adata.obs[["SpatialGlue"]].to_parquet(args.output)


if __name__ == "__main__":
    main()
