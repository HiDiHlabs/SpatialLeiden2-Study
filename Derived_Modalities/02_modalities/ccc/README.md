# CCC (cell-cell communication)

Input: spot x ligand-receptor local scores (LIANA+-style, 0-1, `LIGAND^RECEPTOR`) computed by a team member
for P1 control (10X001), P9 RZ (ACH0023), P9 IZ (ACH0012) and P9 IZ rep2 (ACH0013). Stored as
`data/kuppe_2022/modalities/ccc_<run>.csv` (not in git; copy on OneDrive and in the data bundle, see
`handoff/HANDOFF_2026-11-24_start_here.md`). Method, kernel bandwidth and pair filtering are not documented yet.

**The scores are spatially smoothed upstream.** About a quarter of each spot's CCC feature-space neighbours are
physical neighbours (vs 4-8 % for RNA and TF), so with SpatialLeiden we use spatial weight 0.25 for CCC.

- `ccc_maps.py`: maps of selected pairs (SPP1->CD44, FN1->ITGAV/ITGB1, ...) on the four slides.
- The CCC layer is loaded and graphed by `03_.../exp0_sweep.py` (`load_ccc_layer`, `preprocess_ccc`).
