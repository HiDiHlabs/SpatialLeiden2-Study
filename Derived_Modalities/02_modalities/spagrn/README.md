# SpaGRN spatial regulons

SpaGRN 1.1.4 (env `spagrn`, spec and the one-line pysal patch in `../../envs/`). Regulons are inferred **per
slide** (spatially variable genes, TF-target co-expression, cisTarget motif pruning; activating regulons only).

- `run_spagrn.py <run> data/kuppe_2022 results/kuppe_2022/spagrn --workers 8`: needs `CISTARGET_DIR` (hg38 v10
  rankings, motif table, allTFs_hg38.txt from resources.aertslab.org/cistarget). About 1-2 h per slide.
  Output: `results/kuppe_2022/spagrn/<run>/<run>_spagrn.h5ad` (AUCell in `obsm["auc_mtx"]`, regulons in
  `uns["regulon_dict"]`) + `tmp_files/`. Done for 10X001, ACH0012, ACH0013 (copy on OneDrive).
- `compare_spagrn_collectri.py` (env with decoupler, `../../envs/decoupler.yml`): regulon counts, target overlap,
  spot-wise agreement with CollecTRI ULM, P9 rep2 artery ranking -> `results/kuppe_2022/spagrn/compare_collectri/`.
  Findings in `handoff/HANDOFF_2026-11-24_start_here.md`, section 3.
