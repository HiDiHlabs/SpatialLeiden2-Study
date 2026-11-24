# Work of 24 November 2026 (track C, day 2)

This folder is an index of the work done on 24 November 2026. The scripts themselves stay in their usual folders, because many of them use relative paths and imports that break if they are moved. Paths below are relative to `Derived_Modalities/`; run the scripts from `Derived_Modalities/` as each one's docstring says.

Where to start: the two handoffs, [fig1hij_spatialleiden](../handoff/HANDOFF_2026-11-24_fig1hij_spatialleiden.md) (Kuppe Fig. 1h–j with SpatialLeiden as panel j) and [ccc_panel_j_for_colleague](../handoff/HANDOFF_2026-11-24_ccc_panel_j_for_colleague.md) (panel j on the CCC layer, for another VM).

## Integration across the 27 slides (`07_integration/`)

| File | What it does |
|---|---|
| [exp1_harmony_covariates/exp1_confounding.py](../07_integration/exp1_harmony_covariates/exp1_confounding.py) | Experiment 1, step 0: confounding audit of the candidate Harmony covariates (batch, zone, patient, sample) |
| [exp1_harmony_covariates/exp1_covariates.py](../07_integration/exp1_harmony_covariates/exp1_covariates.py) | Experiment 1: Harmony on RNA with three covariate sets (none, batch, authors'), compared at fixed K |
| [exp1_harmony_covariates/exp1_workers.py](../07_integration/exp1_harmony_covariates/exp1_workers.py) | Worker functions for `exp1_covariates.py` (Leiden K search, UMAP) |
| [exp1_harmony_covariates/NOTES.md](../07_integration/exp1_harmony_covariates/NOTES.md), [RESULT.md](../07_integration/exp1_harmony_covariates/RESULT.md) | Notes and result: the authors' covariates (patient + sample + visium_batch) are used |
| [tf_unintegrated_umap.py](../07_integration/tf_unintegrated_umap.py) | Unintegrated UMAPs of the CollecTRI TF layer: is there slide / patient structure? |
| [tf_harmony.py](../07_integration/tf_harmony.py) | Harmony on the TF layer PCs, compared with the unintegrated layer |
| [tf_variants_compare.py](../07_integration/tf_variants_compare.py) | Three TF-layer variants (none, harmony, per-slide z-score) on batch mixing vs zone preservation |
| [RESULT_tf_harmony.md](../07_integration/RESULT_tf_harmony.md) | Result: the TF layer is integrated with Harmony on its PCs |
| [ccc_harmony.py](../07_integration/ccc_harmony.py) | CCC layer embedding across all 27 slides, built like the TF layer |
| [umaps_rna_tf_ccc.py](../07_integration/umaps_rna_tf_ccc.py) | UMAPs of RNA / TF / CCC by patient and zone, with joint SpatialLeiden domains (Kuppe ED Fig. 1a–b style) |

## Joint SpatialLeiden over all 27 slides (`03_single_modality_clustering/exp2_joint_harmony/`)

| File | What it does |
|---|---|
| [exp2_joint_sweep.py](../03_single_modality_clustering/exp2_joint_harmony/exp2_joint_sweep.py) | Experiment 2: joint SpatialLeiden on Harmony RNA, spatial weight sweep at fixed K = 9 |
| [exp3_joint_multimodal_sweep.py](../03_single_modality_clustering/exp2_joint_harmony/exp3_joint_multimodal_sweep.py) | Experiment 3: RNA + CCC or TF + space, spatial weight sweep at K = 9 |
| [exp3_fast_spatial.py](../03_single_modality_clustering/exp2_joint_harmony/exp3_fast_spatial.py) | Faster finish of experiment 3 (narrower resolution search, one process per weight) |
| [exp4_unfixed_k.py](../03_single_modality_clustering/exp2_joint_harmony/exp4_unfixed_k.py) | Experiment 4: joint SpatialLeiden with K left free across spatial weights |
| [exp4_loader.py](../03_single_modality_clustering/exp2_joint_harmony/exp4_loader.py) | Shared embedding loader (npz / h5ad) for experiments 4 and 5 |
| [exp5_unfixed_k_multimodal.py](../03_single_modality_clustering/exp2_joint_harmony/exp5_unfixed_k_multimodal.py) | Experiment 5: multimodal counterpart of experiment 4 |
| [inspect_panel_j_runs.ipynb](../03_single_modality_clustering/exp2_joint_harmony/inspect_panel_j_runs.ipynb) | Notebook comparing the panel j variants |

## Modalities (`02_modalities/`)

| File | What it does |
|---|---|
| [collectri/collectri_tf_activity_all.py](../02_modalities/collectri/collectri_tf_activity_all.py) | CollecTRI TF activity (decoupler ULM) on all 27 sections |
| [progeny/progeny_kuppe.py](../02_modalities/progeny/progeny_kuppe.py) | PROGENy pathway activity per spot, computed as in Kuppe et al. 2022 |
| [liana/plot_liana_all_samples.py](../02_modalities/liana/plot_liana_all_samples.py) | Plots of LIANA+ local ligand-receptor scores across all samples |

## Kuppe et al. 2022 figure reproductions (`01_data_preprocess/`)

| File | What it does |
|---|---|
| [kuppe_fig1/fig1hij.py](../01_data_preprocess/kuppe_fig1/fig1hij.py) | Fig. 1h–i reproduced, panel j replaced by joint SpatialLeiden domains |
| [kuppe_fig2/fig2b_spatialleiden.py](../01_data_preprocess/kuppe_fig2/fig2b_spatialleiden.py) | Fig. 2b: the authors' niches next to SpatialLeiden domains at increasing spatial weight |
| [kuppe_ed_fig3/ed_fig3abc.py](../01_data_preprocess/kuppe_ed_fig3/ed_fig3abc.py) | Extended Data Fig. 3a–c: composition UMAP by patient / region, niche proportions |
| [kuppe_ed_fig3/ed_fig3d.py](../01_data_preprocess/kuppe_ed_fig3/ed_fig3d.py) | Extended Data Fig. 3d: cell-type and marker-gene feature plots on the composition UMAP |
| [../common/plot_patient_region_umap.py](../common/plot_patient_region_umap.py) | Shared plotting: any embedding by patient and region, in the authors' palettes |

## Scripts changed today

- [exp0_sweep.py](../03_single_modality_clustering/exp0_layer_ratio_sweep/exp0_sweep.py), [mm_configs.py](../04_multimodal_weight_grid/proto_rna_tf/mm_configs.py), [vessel_band_across_modalities.py](../04_multimodal_weight_grid/proto_rna_tf/vessel_band_across_modalities.py), [downloads/kuppe_2022.py](../downloads/kuppe_2022.py), [figures/README.md](../figures/README.md)
- `spatialleiden.yml` (repo root): pins `pynndescent=0.5.13` and `setuptools<81`, adds `datashader==0.19.1`, all needed for the env to run these scripts

Figure summaries (`figures/kuppe_fig1h/*.json`) and the [figure README](../figures/kuppe_fig1h/README.md) are committed; the figure PNGs are not.
