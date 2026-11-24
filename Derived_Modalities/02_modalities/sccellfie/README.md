# scCellFie metabolic task activity

scCellFie 0.6.0 (env `sccellfie`, spec `../../envs/sccellfie.yml`). Metabolic **task activity** (~210 tasks),
not flux. Smoothing uses the RNA expression graph, not spatial neighbours.

1. `save_rna_pca.py` (env spatialleiden): saves the experiment-0 RNA PCA (`rna_pca30_<run>.npy`), because the
   sccellfie env lacks scikit-misc for seurat_v3 HVGs.
2. `run_sccellfie.py <run> data/kuppe_2022` (env sccellfie): writes `sccellfie_tasks_<run>.h5ad`,
   `sccellfie_reactions_<run>.parquet`, `sccellfie_qc_<run>.json` to `data/kuppe_2022/modalities/`.
   The database is in `data/kuppe_2022/sccellfie_db`.
