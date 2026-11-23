# Experiment 0 panels: SpatialLeiden `layer_ratio` sweep

Panels from `03_single_modality_clustering/exp0_layer_ratio_sweep/`. Design: `README.md` there. Answers: `RESULT.md` there.

**Slides:** P9 IZ rep2 (`ACH0013`), P9 IZ (`ACH0012`), P1 control (`10X001`).
**Sweep:** `layer_ratio` 0–2.5, 5 seeds, RNA only.

One subfolder per choice of K, the number of domains each map must have:

| Subfolder | K |
|---|---|
| `K15_brief/` | 15 on all slides (the brief's rule: plain Leiden at resolution 1.0 on `ACH0013`) |
| `K9/` | 9 on all slides |
| `K_authors_celltype_niches/` | Per slide, the number of the authors' cell-type niches on that slide: 8 / 8 / 9 |
| `K_authors_molecular_niches/` | Per slide, the number of the authors' molecular niches: 8 / 11 / 12 |

**Files in each subfolder:**

- `fig1-3_sweep`: coherence above shuffled, ARI to the ratio-0 partition, and seed stability, each against `layer_ratio`.
- `fig5_resolution_sizes_niches`: the spatial resolution found, the largest cluster's share, and ARI to the authors' niches.
- `fig4_domains_<run>.png`: H&E, cell2location maps, and domains at `layer_ratio` 0 / 0.5 / 1 / 2 beside the authors' niches.
- `exp0_sweep.csv`: one row per slide × ratio × seed.
- `exp0_sanity_ratio0.csv`: the ratio-0 vs plain Leiden check.

To rebuild, run `exp0_sweep.py`, then `exp0_figures.py` (about 2 minutes per K setting).
