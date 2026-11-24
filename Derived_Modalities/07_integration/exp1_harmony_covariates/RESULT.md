# Experiment 1 result: Harmony covariates for the integrated expression layer

**The chosen covariates are the authors' (`patient` + `sample` + `visium_batch`), and the other two variants (`none`, `batch`) are reported as robustness.**

Fixed for everything downstream. Cite this file rather than re-deciding; use the same covariate set for the derived layers (PROGENy, TF activity) unless a layer's own diagnostic says otherwise.

Run 24 Nov 2026: 27 sections, 88,704 spots, 11,681 shared genes, one PCA (30 PCs) corrected three ways, Leiden at exactly K = 9 and K = 12 for every variant (seed 0). Column mapping and pipeline details: `NOTES.md`. Outputs: `results/kuppe_2022/exp1_harmony_covariates/`.

## Numbers (pooled; per zone in `exp1_covariate_comparison.csv`)

| Variant | K | ARI batch | ARI sample | ARI zone | ARI cell-type niche | ARI molecular niche | iLISI batch | cLISI zone | silhouette zone |
|---|---|---|---|---|---|---|---|---|---|
| none | 9 | 0.079 | 0.368 | 0.437 | 0.113 | 0.071 | 0.009 | 0.995 | 0.039 |
| batch | 9 | 0.008 | 0.307 | 0.276 | 0.135 | 0.073 | 0.264 | 0.964 | 0.038 |
| **authors** | 9 | 0.017 | 0.043 | 0.087 | 0.118 | 0.138 | 0.634 | 0.633 | −0.056 |
| none | 12 | 0.096 | 0.491 | 0.419 | 0.123 | 0.075 | 0.009 | 0.995 | 0.039 |
| batch | 12 | 0.009 | 0.390 | 0.319 | 0.137 | 0.071 | 0.264 | 0.964 | 0.038 |
| **authors** | 12 | 0.012 | 0.047 | 0.076 | 0.113 | 0.121 | 0.634 | 0.633 | −0.056 |

iLISI and cLISI are scaled to [0, 1] (iLISI 1 = fully mixed; cLISI 1 = neighbourhoods pure in zone). Embedding metrics do not depend on K.

## What depends on the integration choice, and what does not

- **Robust:** agreement with the authors' cell-type niches is the same for all three variants (ARI 0.11–0.14 at both K). The spot-level composition structure does not depend on the correction.
- **Depends strongly:** sample structure. Without sample in the covariates, the clusters largely follow sections: ARI against sample 0.37–0.49 (`none`) and 0.31–0.39 (`batch`), up to 0.98–1.00 within BZ, where each of the three sections forms its own clusters. The authors' covariates remove it (ARI sample ≈ 0.04; iLISI 0.63 vs 0.01 / 0.26).
- **Correcting on `visium_batch` alone mixes the two batches but leaves the per-section islands.** Within BZ (70 % 10X spots), ARI against batch stays at 0.54 at K = 9 for `batch` and drops to 0.004 only with the authors' covariates.
- **ARI against zone falls from 0.44 (`none`) to 0.09 (`authors`).** This is not read as lost biology: zone is a section-level label nested in sample (§0), so a spot clustering can only score high on it by clustering sections, which is what `none` does (see the UMAPs: one island per section). After correction, zone persists as a gradient across a continuous manifold (IZ at one side), not as separate clusters.
- **Molecular niches agree best with `authors`** (0.12–0.14 vs ≈ 0.07). This is expected and partly circular: the authors defined their molecular niches on an expression space integrated this way.

## Descriptive: silhouette-optimal K on the authors' embedding

Silhouette of the Leiden labels is highest at K = 4 (0.170, the lower edge of the K = 4–20 range searched), then 0.151 at K = 5, and falls to 0.065–0.08 from K = 6 on, with no second peak at 9 or 12. The criterion favours the coarsest split; it is not used for the comparison (`exp1_silhouette_optimal_k.csv`).
