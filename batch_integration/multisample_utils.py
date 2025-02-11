from anndata import AnnData

SEED: int = 42


def preprocess_genes(
    adata: AnnData, batch_key: str | None, min_cells: int = 10, genes: int = 2_000
):
    import scanpy as sc

    sc.pp.filter_genes(adata, min_cells=min_cells)
    sc.pp.highly_variable_genes(
        adata, flavor="seurat_v3", n_top_genes=genes, batch_key=batch_key
    )


def preprocess_scaling(adata: AnnData, n_pcs: int = 30, seed: int = SEED):
    import scanpy as sc

    layer = "scaled"
    adata.layers[layer] = adata.X.copy()
    sc.pp.normalize_total(adata, layer=layer)
    sc.pp.log1p(adata, layer=layer)
    # sc.pp.scale(adata, zero_center=False)
    sc.pp.pca(adata, n_comps=n_pcs, zero_center=True, random_state=seed, layer=layer)


def preprocess_neighbors(
    adata: AnnData, obsm_key: str, *, n_neighbors: int = 15, seed: int = SEED
):
    import scanpy as sc

    sc.pp.neighbors(
        adata, n_neighbors=n_neighbors, n_pcs=None, use_rep=obsm_key, random_state=seed
    )
