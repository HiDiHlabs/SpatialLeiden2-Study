# TF layer integration result: Harmony on the TF PCs

**The TF activity layer is integrated with Harmony on its PCs, using the covariates fixed in experiment 1 (`patient` + `sample` + `visium_batch`). The two alternatives (no correction, per-slide z-score) are reported as comparison.**

Decided 24 Nov 2026. Fixed for everything downstream; cite this file. Covariates: `exp1_harmony_covariates/RESULT.md`.

## Recipe

| Step | Setting |
|---|---|
| Input | CollecTRI TF activity, decoupler 2.2.0 ULM (tmin = 5), all 27 sections (`data/kuppe_2022/modalities/collectri_all27/`, `02_modalities/collectri/collectri_tf_activity_all.py`) |
| Features | the 622 TFs scored in all 27 sections; no highly variable feature selection |
| Scaling | z-score per TF over all 88,704 spots pooled (centred, unit variance, no clipping) |
| PCA | 30 PCs, scanpy `sc.pp.pca`, seed 0 |
| Harmony | harmonypy 0.2.0 on the 30 PCs; covariates `patient`, `sample`, `visium_batch`; package defaults: theta 2 and lambda 1 per covariate, sigma 0.1, 100 clusters (min(N/30, 100)), block size 0.05, `max_iter_harmony` 10, `max_iter_kmeans` 20, seed 0. Converged after 4 iterations |
| Graph / UMAP | 15-NN on the 30 Harmony dimensions (euclidean); scanpy UMAP (min_dist 0.5), seed 0 |

Same Harmony settings as the RNA layer (experiment 1). Code: `07_integration/tf_harmony.py`; outputs in `results/kuppe_2022/tf_harmony/` (`tf_embeddings.npz` holds `X_pca` and `X_harmony` in the spot order of `tf_obs.csv.gz`).

## Numbers (pooled; per zone in `tf_variants_compare.csv`)

| Variant | LISI slide | LISI patient | LISI batch | spots with LISI slide < 1.5 | cLISI zone | silhouette zone |
|---|---|---|---|---|---|---|
| none (pooled z-score) | 1.47 | 1.30 | 1.02 | 51 % | 0.91 | 0.009 |
| **Harmony on TF PCs** | 6.47 | 5.27 | 1.73 | 0.1 % | 0.51 | −0.035 |
| per-slide z-score | 6.64 | 5.60 | 1.62 | 1.3 % | 0.50 | −0.023 |

LISI = median effective number of slides / patients / batches in a spot's neighbourhood (1 = all neighbours from one). cLISI zone scaled to [0, 1] (1 = pure). Silhouette on a stratified subsample of 3,000 spots per zone.

## Why Harmony

- **The uncorrected TF layer needs correction.** Neighbourhoods are almost entirely single-batch (LISI batch 1.02), although batch is not confounded with zone (experiment 1, §0); IZ and FZ sections form patient-specific lobes, and half the spots sit in single-slide neighbourhoods.
- **Harmony and per-slide z-scoring tie on the numbers.** Both clear slide structure (LISI slide 6.5 / 6.6) and both cost about the same zone separation (cLISI 0.51 / 0.50). Per zone it is mixed: Harmony keeps more of CTRL (silhouette 0.033 vs 0.005), per-slide slightly more of IZ (cLISI 0.62 vs 0.59). The decision rule of experiment 1 does not separate them.
- **Harmony is chosen for consistency and mechanism, not because it scores better.** It is the recipe of the RNA layer, so differences between layers are not confounded with differences in correction; and it aligns structure shared across slides instead of recentring every section, which removes real between-section differences in mean TF activity along with the slide effect.
- **Zone survives as a gradient.** After Harmony, IZ concentrates on one side of a single continuous manifold and the myogenic zones (control, RZ, BZ) on the other, as for RNA in experiment 1.

## Caveats

- harmonypy is not bit-reproducible across runs with the same seed (multithreaded BLAS in its k-means): reruns differ by a median of 1e-4 (max 0.3) per coordinate, with every PC correlating at 1.000. Use the saved `X_harmony`, not a recomputation, when labels must match.
- The TF scores of the three 23 Nov slides (10X001, ACH0012, ACH0013) are reproduced exactly by the all-27 run (max difference 6e-8).
