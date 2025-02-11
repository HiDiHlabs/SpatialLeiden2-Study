#!/usr/bin/env python

from pathlib import Path

import anndata as ad
import pandas as pd
from abc_atlas_access.abc_atlas_cache.abc_project_cache import AbcProjectCache
from scipy.sparse import csr_matrix

download_base = Path("../data/abc_atlas")
data_path = Path("../data/ABC_atlas")

# download data
abc_cache = AbcProjectCache.from_cache_dir(download_base)

expression_file = abc_cache.get_data_path(
    directory="Zhuang-ABCA-1", file_name="Zhuang-ABCA-1/raw"
)

cell_file = abc_cache.get_metadata_path(
    directory="Zhuang-ABCA-1", file_name="cell_metadata_with_cluster_annotation"
)

ccf_file = abc_cache.get_metadata_path(
    directory="Zhuang-ABCA-1-CCF", file_name="ccf_coordinates"
)


# prepare anndata object
coord_columns = ["x", "y", "z"]
cell_meta_columns = ["neurotransmitter", "class", "subclass"]

adata = ad.read_h5ad(expression_file)
adata.X = csr_matrix(adata.X)

cells = pd.read_csv(
    cell_file,
    index_col="cell_label",
    usecols=["cell_label"] + coord_columns + cell_meta_columns,
)

adata = adata[adata.obs_names.isin(cells.index), :]

adata.obs = adata.obs.join(cells)
adata.obsm["spatial"] = adata.obs[coord_columns].to_numpy()
adata.obs = adata.obs.drop(columns=coord_columns).astype("category")

adata.uns["batch_keys"] = ["brain_section_label"]

adata.write_h5ad(data_path / "Zhuang-ABCA-1.h5ad")
