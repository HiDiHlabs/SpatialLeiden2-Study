# Handoff: TF activity, scCellFie and the first multimodal SpatialLeiden prototype

23 November 2026, late afternoon. This continues `HANDOFF_2026-11-23_cohort_and_first_slides.md` (section 5). Paths are relative to `Derived_Modalities/`. **Nothing in this handoff is pushed yet.**

## 1. Slides

All work below uses the same three Visium sections:

| Label | Run | What it is |
|---|---|---|
| P1 control | `10X001` | Healthy donor myocardium. Shallow (3,264 UMIs/spot), 32 % mitochondrial reads. Negative control |
| P9 IZ | `ACH0012` | Ischaemic zone, 2 days after onset. **Normal myocardium + granulation tissue** (28 % cardiomyocytes, 24 % fibroblasts). Typical depth (3,540 UMIs/spot) |
| P9 IZ rep2 | `ACH0013` | Same patient, same day, **mostly granulation tissue** (12 % cardiomyocytes, 40 % fibroblasts). Very deep (8,582 UMIs/spot). The slide in the paper's Fig. 2b |

The two P9 sections differ in tissue composition and depth, **not in time point or degree of necrosis.** The pathologist notes for both say most necrotic cardiomyocytes are already cleared. They are separate pieces of tissue, not replicates or consecutive sections.

## 2. Modalities derived from expression

| Modality | Tool | Output | Status |
|---|---|---|---|
| TF activity | decoupler 2.2.0 ULM + CollecTRI (env `masld_atlas`, read-only) | `data/kuppe_2022/modalities/collectri_ulm_<run>.h5ad` (~680 TFs) | Done |
| Metabolic task activity | scCellFie 0.6.0 (env `sccellfie`, read-only) | `data/kuppe_2022/modalities/sccellfie_tasks_<run>.h5ad` (~210 tasks; reactions in `.parquet`) | Done, not yet used |
| Spatial regulons | SpaGRN 1.1.4 (env `spagrn`) | `results/kuppe_2022/spagrn/<run>/` | **Running** (screens `spagrn_<run>`; at the Geary's C step at 17:33 CET) |

**Notes:**

- **scCellFie infers metabolic *task activity*, not flux.** Smoothing uses the RNA expression graph, not spatial neighbours, to avoid building spatial structure into the layer. That graph is the experiment-0 RNA PCA, saved by `02_modalities/sccellfie/save_rna_pca.py`, because env `sccellfie` lacks scikit-misc.
- **SpaGRN needed two fixes:**
  - The default `preprocess()` drops spots with ≥ 4,000 genes; the limit is raised to keep all spots.
  - A pysal 2.0 / esda mismatch broke the Geary's C step. It is fixed by a one-line patch to the installed `spagrn/autocor.py`, documented in `envs/spagrn.yml` and `envs/README.md`.
- **Depth confound on P1:** a typical TF has ≥ 5 targets detected in 16 % of P1 spots, against 70 % on P9 rep2. P1's derived layers are therefore exploratory.

## 3. TF activity checks against the biology

**Figure:** `results/kuppe_2022/collectri/figures/collectri_expected_tfs_p1_p9.png`, script `02_modalities/collectri/collectri_expected_tfs.py`.

- **HIF1A (hypoxia) and RELA (NF-κB):** higher on both P9 sections than on P1, as expected for an acute infarct.
- **MEF2C (cardiomyocyte identity):** uniformly high on P1, and on P9 high only in surviving myocardium, a negative image of the lesion.
- **Depth control:** MEF2C goes opposite to the depth difference, so the P1 vs P9 contrast is not a depth artefact.

## 4. SpatialLeiden on TF activity instead of RNA

**Code:** `exp0_sweep.py --latent collectri`, which builds the TF graph as z-score per TF → PCA 30 → 15-NN, with the same spatial graph, grid and K as experiment 0.
**Results:** `results/kuppe_2022/exp0tf_*` (four K settings), compared by `exp0_compare_rna_tf.py` in `exp0tf_k_celltype_niche/compare_rna_tf/`.

- **P9:** TF activity alone carries spatial structure: a coherence gap of 0.22–0.31 without space, against 0.36–0.50 for RNA. With space, TF domains follow the vessel and the lesion.
- **P1:** TF clusters are spatially random (gap 0.03) and unrelated to the RNA clusters (ARI 0.04).
  - They are **not** depth strata: cluster membership explains 13 % of log-UMI variance, the same as for RNA.
  - The agreement between RNA and TF at high spatial weight (about 0.5) comes from the shared spatial tiling, not from shared biology.

## 5. Multimodal prototype: RNA + TF activity + space

**Code:** `04_multimodal_weight_grid/proto_rna_tf/mm_proto.py` (reuses the experiment-0 functions) and `mm_proto_figures.py`.
**Results:** `results/kuppe_2022/mm_proto_rna_tf/`.

**Design:**

| | Setting |
|---|---|
| Layers | RNA graph (exactly as in experiment 0), TF graph, spatial grid; `spatialleiden_multimodal` |
| Weights | RNA 1; spatial 0.5 or 0.75; TF 0, 0.25, 0.5, 1, 2 |
| K | The authors' cell-type niche count per slide (8 / 8 / 9) |
| Seeds | 5 |
| Resolutions | RNA and TF share the expression resolution (searched once per seed); the spatial resolution is searched to hit K |
| Baseline | TF weight 0 (= RNA + space) |

**Where the weights come from:**

- **RNA = 1** is the reference.
- **Spatial 0.5–0.75** is experiment 0's working range. It is consistent with the preprint's 0.8 for single Visium slides, but chosen from our data.
- **The TF grid is our own.** The preprint gives no rule; it names manual weighting as an open problem. Its only multimodal example (Visium HD CRC) used RNA 1.5, H&E 0.5, spatial 2.

**Findings:**

- **P9:** TF weight 0.25–0.5 keeps 60–80 % agreement with RNA + space. Coherence is slightly lower and seed stability equal or better. TF-driven subdivisions appear inside the lesion.
- **P1:** TF weight ≥ 1 breaks the structure (coherence gap 0.61 → 0.12). The noisy TF layer must stay at low weight.
- **Agreement with the authors' niches** does not change at any TF weight.
- **Open design flaw:** 21 of 150 runs missed K, almost all at TF weight 1–2 (fewer domains than the target). Searching only the spatial resolution cannot restore K when the TF layer dominates. The next version should also search the shared expression/TF resolution.

## 6. Summary figure

**Figure:** `results/kuppe_2022/mm_proto_rna_tf/figures/fig_rna_tf_combo.png/.pdf`, script `04_multimodal_weight_grid/proto_rna_tf/fig_rna_tf_combo.py`.

**Layout:**

- Rows: P1 control (`10X001`), P9 IZ (`ACH0012`), P9 IZ rep2 (`ACH0013`).
- Columns: H&E | RNA only | RNA + space | TF activity + space | RNA + TF + space | authors' cell-type niches | SPP1 expression | NPPB expression.

**Settings for every clustering column:** K = the authors' cell-type niche count, seed 0, colours matched to RNA + space. No footnote and no ARI on the maps; settings and coherence are in `fig_rna_tf_combo_settings.csv`.

| Column | RNA | TF | Spatial |
|---|---|---|---|
| RNA only | 1 | — | 0 |
| RNA + space | 1 | — | 0.75 |
| TF activity + space | — | 1 | 0.75 |
| RNA + TF + space | 1 | 0.5 | 0.75 |

**Markers, chosen from the paper:**

- **SPP1:** the paper reports SPP1+ macrophages expanded in the ischaemic group. It is detected in 43–46 % of P9 spots against 8 % on P1. On P9 (`ACH0012`), its high patch coincides with the domain every clustering finds.
- **NPPB:** the paper's stressed-cardiomyocyte (vCM3) marker. It is high in P9's surviving myocardium bordering the lesion, and absent inside the lesion.

**Coherence (edge-weighted):**

| Slide | RNA only | RNA + space | TF + space | RNA + TF + space |
|---|---|---|---|---|
| P1 | 0.38 | 0.78 | 0.72 | 0.67 |
| P9 rep2 | 0.59 | 0.73 | 0.73 | 0.68 |
| P9 | 0.75 | 0.91 | 0.80 | 0.88 |

## 7. Other changes since the last push

- **Code review fixes to experiment 0** (one coherence definition, misses excluded from means, fresh fixed-seed shuffles, no private API, explicit resolution search, gene-order assert). Partitions are unchanged and numbers move by at most 0.003–0.027. Details are in `RESULT.md`.
- **The intro notebook** reruns cleanly with these fixes.
- **Experiment-0 panels** in `figures/exp0_layer_ratio/` are refreshed, with `runs_missing_k_per_ratio.csv`.

## 8. Next steps

1. **SpaGRN:** when it finishes, map the regulon activities (AUCell) on the three slides and compare its regulons with the CollecTRI TF activities. SpaGRN infers regulons per slide, so TF sets will differ between slides.
2. **Multimodal v2:**
   - search the shared expression/TF resolution as well, so high TF weights reach K;
   - add the scCellFie layer;
   - optionally one run at the preprint's CRC proportions (RNA 1.5 : TF 0.5 : spatial 2), for reference.
3. **Overnight robustness grid** (experiment 1; design in the previous handoff, section 5.7). Estimated at 1–3 h staged, 4–10 h fuller.
4. **Push**, when you say so: everything above (scripts, READMEs, notebook, panels). Not in git: data, the slide atlas PDF, and crops of published figures.
