"""Save the experiment-0 RNA PCA (30 components) per slide for scCellFie's smoothing graph.

Runs in env spatialleiden (seurat_v3 HVGs need scikit-misc, which env sccellfie lacks).
Same steps as exp0_sweep.py: genes in >=10 spots, 2,000 HVGs (seurat_v3, counts),
normalize_total, log1p, scale without centring, PCA 30 (random_state 0).
Usage (from Derived_Modalities/): python 02_modalities/sccellfie/save_rna_pca.py data/kuppe_2022
Output: data/kuppe_2022/modalities/rna_pca30_<run>.npy (rows in h5ad barcode order)
"""
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import scanpy as sc

H5AD = {"10X001": "Visium_control_P1.h5ad", "ACH0013": "Visium_GT_IZ_P9_rep2.h5ad", "ACH0012": "Visium_GT_IZ_P9.h5ad"}
DATA = Path(sys.argv[1])
for run, f in H5AD.items():
    h5 = ad.read_h5ad(DATA / "processed1" / f)
    a = ad.AnnData(X=h5.raw.X.copy(), var=h5.raw.var.copy())
    a.obs_names = h5.obs_names
    a.layers["counts"] = a.X.copy()
    sc.pp.filter_genes(a, min_cells=10)
    sc.pp.highly_variable_genes(a, n_top_genes=2000, flavor="seurat_v3", layer="counts")
    sc.pp.normalize_total(a)
    sc.pp.log1p(a)
    sc.pp.scale(a, zero_center=False)
    sc.pp.pca(a, n_comps=30, mask_var="highly_variable", random_state=0)
    np.save(DATA / "modalities" / f"rna_pca30_{run}.npy", a.obsm["X_pca"].astype(np.float32))
    print(run, a.obsm["X_pca"].shape)
