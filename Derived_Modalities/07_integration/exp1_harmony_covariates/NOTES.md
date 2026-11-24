# Experiment 1 notes: Harmony covariates for the integrated expression layer

## Data facts that constrain the pipeline

- **No mitochondrial genes in any section.** None of the 27 processed Visium h5ads (Zenodo 6578047) contains `MT-` genes, in `X` or in `.raw`. `percent.mt` therefore comes from the authors' metadata (spot tables) and cannot be recomputed, regressed out, or used as a filter.
- **Gene sets differ per section** (13,144 to 17,709 genes in `.raw`). The joint object keeps the 11,681 genes shared by all 27 sections and renormalises their raw counts.
- **27 runs, not 28.** ACH0021 (the second RZ_P6 run) has no processed h5ad; the paper's sample-level analyses also use n = 27. The §0 section-level crosstab counts all 28 runs. Without ACH0021, ACH RZ goes from 5 to 4 sections, which does not change the non-confounding conclusion.

## Column mapping (spec name → authors' column)

| Spec | Authors | Values |
|---|---|---|
| `visium_batch` | `batch` | 10X (8 runs), ACH (19 of the 27) |
| `zone` | `major_labl` | CTRL, RZ, BZ, IZ, FZ |
| `sample` | `hca_sample_id` | one per Visium run |
| `patient`, `patient_group` | same | |

## §0 confounding audit (see `exp1_confounding.csv`, `exp1_confounding_summary.csv`)

- Batch × zone is not near-diagonal (Cramér's V 0.32 over sections, χ² p = 0.59). The gate passes.
- In spots, 70 % of BZ spots are 10X (2 of the 3 BZ sections), so metrics are stratified by zone.
- Patient is nearly nested in batch (V = 0.94), so a `batch + patient` variant was dropped.
- `patient_group` nearly duplicates zone (group_2 = all IZ; group_3 = five of the six FZ sections), so ARI against `patient_group` is not reported.

## Design decisions (settled 24 Nov 2026)

- Primary pipeline = the authors' covariates (patient + sample + batch). `none` and `batch` are sensitivity variants: the question is whether the result depends on the integration choice.
- One PCA (genes above; HVG seurat_v3 2,000 on counts with `batch_key="sample"`; normalise, log1p, scale without centring, 30 PCs), corrected three ways with harmonypy 0.2.0 (default parameters, seed 0).
- K fixed: 9 (published cell-type niches) and 12 (molecular niches). Variants are compared within K only.
- The silhouette-optimal K on the integrated embedding is reported as a descriptive result only.
- Metrics: ARI (sklearn), iLISI / cLISI from `harmonypy.compute_lisi` (perplexity 30), sklearn silhouette on a stratified subsample of 3,000 spots per zone. No scib-metrics (it would upgrade numpy to 2.x in this env). Every metric is reported pooled and per zone.
- Figure 2 colours the UMAPs by batch, zone and the authors' cell-type niche (not patient group, which duplicates zone).
