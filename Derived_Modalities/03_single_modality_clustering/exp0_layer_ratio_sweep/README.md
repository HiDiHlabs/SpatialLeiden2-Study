# Experiment 0: SpatialLeiden `layer_ratio` sweep on single slides

**Question.** On a single slide, with RNA only, how does the partition change as a function of the spatial layer's weight alone?

SpatialLeiden always runs on two graphs: the expression (latent) graph and the spatial graph. This makes it multiplex, not multimodal: space is geometry, not a modality. This experiment uses `spatialleiden.spatialleiden()`. It does not use `spatialleiden_multimodal()`, a composition layer, derived modalities or several samples.

## Design

| | Value |
|---|---|
| Slides | `GT/IZ_P9_rep2` (`ACH0013`) and `GT/IZ_P9` (`ACH0012`): same patient and lesion, about 2.4× difference in depth. Plus `control_P1` (`10X001`), for comparison with a colleague's run |
| Modality | RNA only (raw counts from the Zenodo h5ad, joined by barcode to the authors' per-spot table) |
| Expression graph | Genes in ≥ 10 spots → 2,000 HVGs (`seurat_v3`, counts) → `normalize_total` → `log1p` → `scale(zero_center=False)` → PCA 30 → `sc.pp.neighbors(k=15)`. Built once per slide, so k is fixed |
| Spatial graph | `sq.gr.spatial_neighbors(coord_type="grid", n_neighs=6)` |
| `directed` | `(False, False)`, as in the preprint's code (`clustering/spatialleiden_multisample.py`). The scanpy graph is symmetric and the grid is regular |
| `n_iterations` | `-1` (Leiden runs to convergence) |
| `layer_ratio` | 0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5 |
| Seeds (`random_state`) | 0–4 at every `layer_ratio` |
| Number of clusters K | **Primary: 15** on all three slides. This is the ratio-0 run (plain Leiden, resolution 1.0, seed 0) of the reference slide `ACH0013`. At resolution 1.0, `ACH0012` gives 11 and `10X001` gives 12 |
| Holding K fixed | At ratio 0, the expression resolution is searched with `spatialleiden(layer_ratio=0)` itself, separately for each seed. At ratio > 0, the expression resolution stays at that seed's ratio-0 value, and `search_resolution_spatial` finds the spatial resolution that gives K. Both searches run up to 30 iterations |
| Sanity check | `sl.leiden` against `sl.spatialleiden(layer_ratio=0)` at the same resolution and seed. **It failed in the first attempt; see below** |
| Baseline | The ratio-0 run with the same seed. The published niches are context only: they were pooled over 27 slides, built from composition rather than expression, with Louvain on scran graphs, and without coordinates |

**Sensitivity run at K = 9.** K = 15 is more clusters than the authors' 9 cell-type niches (pooled over 27 slides), and more than plain Leiden gives on the other two slides at resolution 1.0 (11 and 12). The whole sweep is therefore repeated at K = 9, the number of published cell-type niches and the value used in the intro notebook. Outputs go to `results/kuppe_2022/exp0_k9/`.

## Sanity check: `layer_ratio = 0` is not plain Leiden (finding, 23 Nov 2026)

On `ACH0013`, for seeds 0–2:

- **Same objective, same result without the spatial graph.** `sl.leiden` and `sl.multiplex_leiden` on the expression layer alone give identical partitions (ARI = 1.000).
- **Adding the spatial graph at weight 0 changes the result anyway.** `sl.spatialleiden(layer_ratio=0)` agrees with plain Leiden at only ARI 0.76–0.84, and can give a different cluster count (13 instead of 15 for seed 0).

The zero-weight layer contributes nothing to the objective. However, leidenalg's multiplex optimiser still considers moves along its edges, so it searches differently and ends in a different local optimum of the same objective.

**Consequence: the baseline is `spatialleiden(layer_ratio=0)` itself,** with its own resolution search, so every ratio, including 0, uses the same optimiser and the same K. `exp0_sanity_ratio0.csv` reports, per seed, the ARI of plain Leiden against both latent-only multiplex and ratio 0, beside the seed-to-seed ARI of plain Leiden. That last number is the noise level the ratio-0 difference should be compared with.

The first attempt, which used plain-Leiden resolutions, was stopped and discarded.

A separate small bug: `multiplex_leiden` with its default `layer_weights=1` (an int) raises `TypeError`, because the length check expects a float. Pass `1.0`. Worth reporting upstream.

## Measures, for every (slide, layer_ratio, seed)

| Measure | What it shows |
|---|---|
| Spatial coherence gap | Share of neighbouring spot pairs with the same label, minus the same after shuffling labels within the slide (mean of 20 shuffles) |
| ARI to ratio 0 | How far the partition moved from the expression-only baseline (same seed) |
| ARI and NMI to the published niches | Context only, not a score |
| Seed stability | Mean, min and max pairwise ARI among the 5 seeds |
| Cluster sizes | Largest and smallest cluster, and the number of clusters below 1 % of spots |
| Resolutions | The resolutions found, and whether the search reached K |

## Pre-specified interpretation

Written on 23 November 2026, before the sweep was run.

- **Space helped:** the coherence gap rises well above its ratio-0 value, the partition stays related to the expression baseline, and domains line up with the infarct core, the vessel and the fibroblast-rich region on the H&E and cell2location maps.
  - Operationally, at some `layer_ratio`: gap ≥ ratio-0 gap + 0.15, ARI to ratio 0 ≥ 0.3, and seed ARI ≥ 0.6.
- **Space did nothing:** the partition barely moves. Operationally, ARI to ratio 0 stays ≥ 0.8, and the gap changes by < 0.05 across the whole sweep.
- **Space over-smoothed:** coherence rises, but ARI to ratio 0 collapses and domains become large blobs unrelated to visible tissue. Operationally, ARI to ratio 0 < 0.2 with domains that do not follow the H&E or cell2location structure.
- **P1 (control, little visible structure):** we expect a smaller coherence gain than on P9 at the same ratio, and no tissue-like domains at moderate ratios. Clear domains on P1 would need checking against its H&E before being believed.
- **P9 depth:** if the curves of `ACH0013` (deep) and `ACH0012` (typical depth) differ in shape, the effect of space depends on depth.

These operational thresholds are our own, set before seeing any results. They are not from the brief or the paper.

## Files

- `exp0_sweep.py` runs the sweep. Outputs go to `results/kuppe_2022/exp0/`: `exp0_sweep.csv`, `exp0_sanity_ratio0.csv`, per-slide label tables and the log.
- `exp0_figures.py` makes figures 1–4.
- `RESULT.md` answers the experiment's questions after the run.
