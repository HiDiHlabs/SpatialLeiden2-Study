"""Metabolic task activity per Visium spot with scCellFie (0.6.0) on a Kuppe 2022 slide.

scCellFie scores metabolic genes per spot, combines them through gene-protein-reaction
rules into reaction activities, and those into metabolic task activities (Recon2.2-based
tasks). This is inferred metabolic *activity*, not flux: no flux balance is solved.

Recipe (as in our earlier MASLD scCellFie runs):
  * raw counts in X; n_counts = the spot's total counts (authors' nCount_Spatial)
  * scCellFie database: local copy in data/kuppe_2022/sccellfie_db (Recon2.2 tasks,
    thresholds, gene-reaction tables)
  * smoothing (alpha 0.33, the default) over a 10-NN graph in EXPRESSION space
    (PCA of the RNA, built as in experiment 0, saved by save_rna_pca.py). Not over spatial neighbours: that would
    put spatial structure into the metabolic layer before SpatialLeiden adds space,
    making any later comparison circular.

Usage (env sccellfie, from Derived_Modalities/):
    python 02_modalities/sccellfie/run_sccellfie.py <run_id> data/kuppe_2022
Output: data/kuppe_2022/modalities/sccellfie_tasks_<run>.h5ad (spots x tasks, same spot
order as the RNA h5ad), sccellfie_reactions_<run>.parquet, sccellfie_qc_<run>.json
"""

import json
import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import sccellfie

H5AD = {"10X001": "Visium_control_P1.h5ad", "ACH0013": "Visium_GT_IZ_P9_rep2.h5ad",
        "ACH0012": "Visium_GT_IZ_P9.h5ad"}
run, DATA = sys.argv[1], Path(sys.argv[2])
OUT = DATA / "modalities"
OUT.mkdir(exist_ok=True)
t0 = time.time()

h5 = ad.read_h5ad(DATA / "processed1" / H5AD[run])
spots = pd.read_csv(DATA / "visium_tables" / f"{run}_spots.csv").set_index("barcode")
assert set(spots.index) == set(h5.obs_names), "barcodes differ"
spots = spots.loc[h5.obs_names]
A = ad.AnnData(X=sp.csr_matrix(h5.raw.X), obs=spots.copy(), var=pd.DataFrame(index=h5.raw.var_names))
A.obs["n_counts"] = spots["nCount_Spatial"].to_numpy(np.float32)
assert np.allclose(np.asarray(A.X.sum(1)).ravel(), A.obs["n_counts"]), "counts mismatch"
A.obsm["spatial"] = spots[["pxl_col_fullres", "pxl_row_fullres"]].to_numpy(float)

# expression kNN graph for smoothing: the experiment-0 RNA PCA, saved by save_rna_pca.py (env spatialleiden)
A.obsm["X_pca_rna"] = np.load(DATA / "modalities" / f"rna_pca30_{run}.npy")
assert A.obsm["X_pca_rna"].shape == (A.n_obs, 30)
sc.pp.neighbors(A, n_neighbors=10, use_rep="X_pca_rna", random_state=0)
print(run, A.shape, f"graph {time.time() - t0:.0f}s", flush=True)

res = sccellfie.run_sccellfie_pipeline(
    A, organism="human", sccellfie_data_folder=str(DATA / "sccellfie_db"), n_counts_col="n_counts",
    process_by_group=False, groupby=None, neighbors_key="neighbors", n_neighbors=10, batch_key=None,
    threshold_key="sccellfie_threshold", smooth_cells=True, alpha=0.33, chunk_size=5000,
    disable_pbar=True, save_folder=None, save_filename=None, verbose=True)
R = res["adata"]
T = R.metabolic_tasks
X = T.X.toarray() if sp.issparse(T.X) else np.asarray(T.X)
assert np.isfinite(X).all()
assert list(T.obs_names) == list(h5.obs_names), "spot order changed"
out = ad.AnnData(X=X.astype(np.float32), obs=spots.copy(), var=pd.DataFrame(index=T.var_names))
out.obsm["spatial"] = A.obsm["spatial"]
out.uns["sccellfie"] = {"version": sccellfie.__version__, "smoothing": "alpha 0.33, 10-NN expression PCA graph",
                        "n_counts": "nCount_Spatial", "database": "data/kuppe_2022/sccellfie_db"}
out.write_h5ad(OUT / f"sccellfie_tasks_{run}.h5ad", compression="gzip")
Rx = R.reactions
RX = Rx.X.toarray() if sp.issparse(Rx.X) else np.asarray(Rx.X)
pd.DataFrame(RX.astype(np.float32), index=Rx.obs_names, columns=Rx.var_names).to_parquet(
    OUT / f"sccellfie_reactions_{run}.parquet")
qc = {"run": run, "spots": int(out.n_obs), "tasks": int(out.n_vars),
      "tasks_zero_variance": int((X.std(0) == 0).sum()), "reactions": int(RX.shape[1]),
      "metabolic_genes_kept": int(R.n_vars), "runtime_s": round(time.time() - t0)}
(OUT / f"sccellfie_qc_{run}.json").write_text(json.dumps(qc, indent=1))
print(qc, flush=True)
