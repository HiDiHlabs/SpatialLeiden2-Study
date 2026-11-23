# Experiment 0: result

Run 23 November 2026: 3 slides × 8 `layer_ratio` values × 5 seeds, at K = 15 (primary, the brief's rule) and K = 9 (sensitivity). Outputs are in `results/kuppe_2022/exp0/` (K = 15) and `results/kuppe_2022/exp0_k9/` (K = 9): `exp0_sweep.csv`, `exp0_sanity_ratio0.csv` and `figures/`.

- **Target hit:** K was reached in 113 of 120 runs at K = 15 and 119 of 120 at K = 9. The misses are flagged in `search_failed`.
- **Runtime:** about 2 minutes per K, with one process per slide × seed.

Numbers below are means over 5 seeds, at K = 9 unless stated. K = 15 gives the same picture.

## Answers to the brief's questions

**1. Does `layer_ratio = 0` reproduce plain Leiden?**

Not exactly, but the difference is no larger than ordinary seed noise.

- **Same objective:** plain Leiden and latent-only `multiplex_leiden` are identical (ARI = 1).
- **Different run:** `spatialleiden(layer_ratio=0)` agrees with plain Leiden at ARI 0.55–0.84 (per slide: P1 0.55, P9 rep2 0.68, P9 0.84). Plain Leiden agrees with *itself* across seeds at only 0.47–0.77 (P1 0.47, P9 rep2 0.64, P9 0.77).
- **Why:** the zero-weight spatial layer changes the optimiser's search path, not the objective (README, "Sanity check").
- **Consequence:** the baseline is `spatialleiden(layer_ratio=0)` itself, so every ratio uses the same optimiser. With that, ratio 0 is a valid no-space baseline.

**2. Over what range does coherence rise meaningfully?**

Between `layer_ratio` 0.25 and 1.0. Above 1.5 it saturates: the gap above shuffled reaches about 0.79–0.83 at K = 9 (0.84–0.86 at K = 15) on **all three slides**, including the control.

**3. Where do domains stop tracking expression?**

ARI to the ratio-0 partition falls steeply:

| Slide | ARI at 0.5 | ARI at 1.0 | Below 0.3 from |
|---|---|---|---|
| P1 control | 0.31 | 0.11 | 0.75 (0.15) |
| P9 rep2 (ACH0013) | 0.46 | 0.27 | 1.0 |
| P9 (ACH0012) | 0.54 | 0.43 | Not reached; it stays at 0.36–0.43 |

**4. Is it stable across seeds?**

Moderately, and the baseline itself is no more stable.

- **P9:** mean pairwise ARI 0.54–0.83 at every ratio.
- **P1:** low without space (0.39–0.50 up to ratio 0.75). It rises to about 0.64–0.69 at ratio ≥ 1, where space imposes the same tiling regardless of the seed.

**5. Same shape on both P9 sections, or depth-dependent?**

Same shape. The shallower section (`ACH0012`, about 3,500 UMIs per spot) starts more coherent (gap 0.47 vs 0.36 at ratio 0) and stays closer to its baseline than the deeper section (`ACH0013`, about 8,600 UMIs). So greater depth did not mean more spatial structure in expression here.

**6. At the working range, do domains follow the H&E?** (`figures/fig4_domains_*.png`, seed 0)

- **P9 rep2 at 0.5:** domains follow the tissue.
  - A narrow domain runs along the vSMC-lined vessel.
  - Cardiomyocyte-rich tissue at the top forms one domain.
  - The fibroblast- and myeloid-rich lower right forms its own domains.
  - At 1.0 the vessel domain survives, but the rest coarsens into blocks. At 2.0 the slide is cut into compact, tile-like patches, a Voronoi-like tessellation that ignores most of the visible structure.
- **P1 control:** without space, the clusters are speckled, with one broad region. From 1.0 upwards the slide is tiled into similar-sized patches that correspond to nothing in the H&E or the cell2location maps.

## Verdict against the pre-specified interpretation

| Slide | Verdict | Evidence |
|---|---|---|
| P9 (`ACH0012`) | **Space helped**, at 0.5–1.0 | At 0.5: gap +0.19 over ratio 0, ARI to ratio 0 = 0.54, seed ARI 0.80 |
| P9 rep2 (`ACH0013`) | **Space helped**, borderline, at about 0.5 | At 0.5: gap +0.14 (threshold 0.15), ARI 0.46, seed ARI 0.66. Domains visibly follow the vessel and the lesion |
| Both P9 sections | **Over-smoothed** above about 1.0–1.5 | Tile-like domains |
| P1 control | **Over-smoothed** from 0.75 | ARI to ratio 0 < 0.2; domains unrelated to the tissue |

- **The disease slide has more spatial structure than the control,** as expected. That shows **without any spatial term**: the coherence gap at ratio 0 is 0.36–0.47 on P9 against 0.13 on P1.
- **Coherence alone cannot tell real structure from imposed structure.** At high ratios, P1 reaches the same coherence as P9. The useful signal is coherence gained *while* the partition stays related to the expression baseline and to visible tissue.

**Working range for the next experiments:** `layer_ratio` ≈ 0.5, up to about 0.75. The preprint used 0.8 for single Visium slides.

## Notes

- **The authors' niches are context, not truth.** ARI to their cell-type niches is about 0.13–0.21 on P9 and only 0.02–0.04 on P1. On P1 their niches are themselves close to spatial noise, with coherence near chance, as in the exploration report.
- **The cluster count matters.** At K = 15 the same patterns appear, with more fragmentation, and P1's baseline is even less stable (seed ARI 0.44).
- **One slide each, five seeds:** these are observations, not a benchmark.

## Cluster number matched to the authors (added the same day)

To check whether the low agreement with the published niches is just a matter of cluster number, the sweep was rerun with K set per slide to the number of the authors' niches on that slide:

- **Matched to cell-type niches** (`--k-from celltype_niche`): K = 8 on P1, 8 on P9, 9 on P9 rep2. Output: `results/kuppe_2022/exp0_k_celltype_niche/`.
- **Matched to molecular niches** (`--k-from molecular_niche`): K = 8, 11 and 12. Output: `results/kuppe_2022/exp0_k_molecular_niche/`.

ARI to the published niches at `layer_ratio` 0.5 (cell-type / molecular):

| Slide | K = 15 | K = 9 | K = cell-type niches | K = molecular niches |
|---|---|---|---|---|
| P1 control | 0.02 / 0.02 | 0.03 / 0.03 | 0.04 / 0.04 | 0.04 / 0.04 |
| P9 (ACH0012) | 0.13 / 0.14 | 0.14 / 0.16 | 0.16 / 0.18 | 0.16 / 0.16 |
| P9 rep2 (ACH0013) | 0.16 / 0.10 | 0.21 / 0.14 | 0.21 / 0.14 | 0.17 / 0.11 |

**Matching the cluster number barely changes the agreement.**

- At best it reaches about 0.2 on P9, and it stays at about 0.04 on P1. So the gap to the published niches is not caused by the cluster number. It comes from what is clustered and how: their niches group spots by cell-type composition (or pooled expression) across 27 slides, without coordinates. Ours use one slide's expression plus space.
- Coherence, distance from the baseline and seed stability show the same pattern as at K = 15 and K = 9. The conclusions above do not depend on K.

**Reproducibility check:** where two runs share the same K on the same slide, the numbers are identical. P9 rep2 has K = 9 in both the "K = 9" and the "cell-type niches" runs, and P1 has K = 8 in both matched runs.
