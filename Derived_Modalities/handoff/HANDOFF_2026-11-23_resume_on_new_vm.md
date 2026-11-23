# Handoff: resume track C on a new VM (24 November 2026)

State at the end of 23 November 2026, branch `hack/c-derived-modalities`. Paths are relative to `Derived_Modalities/` unless stated. Earlier handoffs, in order: `HANDOFF_2026-11-23.md`, `HANDOFF_2026-11-23_cohort_and_first_slides.md`, `HANDOFF_2026-11-23_multimodal_prototype.md`. This one covers what came after them and how to pick up.

## 1. What is in git and what is not

| In git (this push) | Not in git |
|---|---|
| All scripts (`01_`–`05_`), READMEs, handoffs, the intro notebook | Kuppe data (`data/`), tens of GB, from Zenodo |
| Env specs: `../spatialleiden.yml`, `envs/spagrn.yml`, `envs/sccellfie.yml` | Derived layers (`data/kuppe_2022/modalities/`, 154 MB) |
| **Cluster labels** of every run: `results/kuppe_2022/**/labels_*.csv.gz`, and the summary CSVs | The colleague's original CCC tables |
| Key figures (PNG): the multipanels, Kuppe Fig. 3 comparison, vessel panels | SpaGRN outputs (2.4 GB), PDFs, `panels_2026-11-23/` (a local copy folder) |

The `results/` folder is git-ignored. The labels, CSVs and PNGs listed above were added explicitly. With them, every figure script runs without re-clustering, provided the data are in place.

## 2. Getting the data onto the new VM

**Option A (fastest): copy the bundle made on the old VM.** It is in the folder `transfer_2026-11-23/`, next to the repo clone (outside the repo), with `SHA256SUMS`.

| File | Size | Contents |
|---|---|---|
| `derived_modalities_data_2026-11-23.tar` | 1.1 GB | `data/kuppe_2022/` subset: the 4 Visium h5ads (P1 control, P9 IZ, P9 IZ rep2, P9 RZ), `visium_tables/` (spot tables with niches, cell2location, coordinates, H&E low-res), `modalities/` (TF, CCC, scCellFie, RNA PCA), `sccellfie_db/`, `metadata/` |
| `derived_modalities_spagrn_2026-11-23.tar` | 2.5 GB | `results/kuppe_2022/spagrn/` (optional) |

Unpack both from inside `Derived_Modalities/` with `tar -xf <file>`, then check with `sha256sum -c SHA256SUMS`. **This is the only copy of the colleague's CCC tables outside their own machine** (`data/kuppe_2022/modalities/ccc_<run>.csv`).

**Option B: rebuild from the source.**

1. `downloads/kuppe_2022.py --sets processed1 visium_raw metadata` (about 43 GB).
2. `01_data_preprocess/kuppe_visium_exploration/02_extract_spaceranger_files.py` to make `visium_tables/`.
3. `02_modalities/collectri/collectri_tf_activity.py` for TF activity (decoupler 2.2, a few minutes).
4. `02_modalities/sccellfie/save_rna_pca.py`, then `run_sccellfie.py <run> data/kuppe_2022` for each slide (env `sccellfie`).
5. **CCC cannot be rebuilt:** ask the colleague for the four `*_CCC_scores.csv` files and copy them to `data/kuppe_2022/modalities/ccc_<run>.csv` (10X001, ACH0012, ACH0013, ACH0023).

## 3. Environments

| Env | Spec | Used for |
|---|---|---|
| `spatialleiden` | `../spatialleiden.yml` plus SpatialLeiden 0.4 (the package, editable install) | All clustering, figures, SIS |
| decoupler | any env with decoupler 2.2.0 and scanpy (we used an existing one) | `02_modalities/collectri/collectri_tf_activity.py` only |
| `sccellfie` | `envs/sccellfie.yml` (scCellFie 0.6.0; export of the working env) | `02_modalities/sccellfie/run_sccellfie.py` |
| `spagrn` | `envs/spagrn.yml` plus the one-line `autocor.py` patch in `envs/README.md` | `02_modalities/spagrn/run_spagrn.py`. Needs `CISTARGET_DIR` pointing at the cisTarget hg38 v10 files |

Run everything from `Derived_Modalities/`. The scripts set one BLAS thread per worker and use process pools; `--jobs` / the last argument sets the number of workers.

## 4. What was done today after the multimodal-prototype handoff

**Single-layer sweeps on the new layers** (same design as experiment 0: K = authors' cell-type niches, layer_ratio 0–2.5, 5 seeds):

- `results/kuppe_2022/exp0ccc_k_celltype_niche`: CCC (`exp0_sweep.py --latent ccc`)
- `results/kuppe_2022/exp0met_k_celltype_niche`: scCellFie tasks (`--latent sccellfie`)

**Named multimodal configurations.** `04_multimodal_weight_grid/proto_rna_tf/mm_configs.py` runs four configs: RNA + TF, RNA + CCC, RNA + metabolic, and RNA + all three. Weights are RNA 1, each extra layer 0.5, spatial 0.75.

- **Results:** `results/kuppe_2022/mm_configs`. 60 runs, all reached K.
- **K-search fallback:** if the spatial resolution cannot reach K, the shared expression resolution is searched instead. Which runs used it is recorded per run.
- **Settings through environment variables:** `MM_SPATIAL_W` and `MM_CONFIGS`. Example: `MM_SPATIAL_W=0.25 MM_CONFIGS=rna_ccc` → `results/kuppe_2022/mm_configs_sp0.25`.

**CCC is pre-smoothed, so it gets spatial weight 0.25.**

- About a quarter of each spot's CCC feature-space neighbours are also physical neighbours, against 4–8 % for RNA and TF.
- At spatial weight 0.75, CCC domains turn into large blobs.
- At 0.25, CCC + space is about as smooth as RNA + space at 0.75. RNA + CCC + space agrees with RNA + space at least as well as at 0.75 (ARI 0.70 on P9 IZ, 0.57 on rep2), and is less blobby.
- The comparison table is `results/kuppe_2022/mm_configs_sp0.25/compare_ccc_spatial_weight.csv`.

**Figures** (`fig_modality_combo.py`, output in `results/kuppe_2022/mm_configs/figures/`):

| Panel argument | Figure | Shows |
|---|---|---|
| `ccc` | `fig_ccc_combo` | RNA vs CCC |
| `ccc_sp025` | `fig_ccc_sp025_combo` | CCC at spatial weight 0.25 and 0.75 |
| `met` | `fig_met_combo` | RNA vs scCellFie metabolic tasks (fatty-acid oxidation, anaerobic glycolysis) |
| `all` | `fig_all_combo` | All modalities side by side |

The methods-slide brief for Claude in PowerPoint is `handoff/HANDOFF_2026-11-23_methods_slide_rna_tf_ccc.md`.

**The artery on P9 IZ rep2 (ACH0013).**

- **The vessel domain:** RNA + TF + space finds a diagonal band (domain 6, 266 spots). It carries MYH11, ACTG2 and GJA5 (an arterial marker), MYOCD/SRF TF activity, and JAG1 → NOTCH3 signalling; 22 % of it is vascular smooth-muscle cells. Figure: `fig_p9rep2_vessel_band`.
- **The authors' molecular niches already contain it.** 51 of the 55 artery-wall spots are their molecular niche 11, which their Fig. 3b points out. Only their cell-type niches split the vessel. Our domain adds continuity: one band of 266 spots, against 86 scattered niche-11 spots.
- **Nested SpatialLeiden inside the band** (`vessel_band_subclusters.py`) splits it reproducibly (5/5 seeds) into artery wall vs a mixed margin.
  - The margin has the **same** fibroblast markers as the surrounding tissue (COL1A1, COL3A1, DCN, LUM, POSTN). It differs only by fading smooth-muscle signal (MYH11: wall 2.9 → margin 1.0 → just outside 0.3 → rest of slide 0.1).
  - **It is not a fibrotic sheath.** At K = 3 the extra cluster is the small islands off the band.
  - Figure: `fig_p9rep2_vessel_subclusters`.
- **Which modalities find the vessel** (`vessel_band_across_modalities.py`):
  - TF alone finds the core cleanly (97 % pure); CCC alone finds it less cleanly.
  - CCC + space and the metabolic layer lose it.
  - All RNA combinations keep it; RNA + TF is the tightest.

**Kuppe Fig. 3 on our results** (`fig_kuppe_fig3_compare.py` → `results/kuppe_2022/kuppe_fig3_compare/`):

- Authors' molecular niches next to our RNA + TF + space domains. Our domains are coloured by the best-matching molecular niche.
- Also: UMAPs, composition heatmaps with one-sided tests, niche shares, and the myogenic niches with MYBPC3 / ANKRD2 on P1 and P9 IZ.
- Our domains are contiguous where the authors' niches are speckled.

**Synergistic Information Score (SIS)** (Richter et al. 2026, Nat. Mach. Intell.; the PDF is in `../../refs/` on the old VM, not in git):

- **Method** (`05_synergy_sis/sis_layers.py` → `results/kuppe_2022/sis/sis_layers.csv`): linear probes on RNA PCA, the other layer, and the two concatenated, with 5-fold spatial-block cross-validation and a shuffled-layer null.
- **P9: TF, CCC and metabolic tasks are redundant with RNA** (SIS −0.06 to +0.02, within the null range).
- **RNA + CCC is positive (+0.15 to +0.17) only for neighbourhood composition.** That reflects CCC's built-in spatial smoothing, not synergy.
- **P1's positive values come from a weak RNA baseline.**
- SpatialLeiden fuses graphs, not embeddings, so SIS here tests the layers, not our clusterings.

**SpaGRN finished on all three slides** (`results/kuppe_2022/spagrn/<run>/`; in the optional bundle). Not analysed yet.

## 5. Next steps

1. **SpaGRN:** map the regulon AUCell scores and compare them with the CollecTRI TF activities. Does MYOCD (or SRF) come out as a regulon in the P9 rep2 artery?
2. **The overnight robustness grid (experiment 1).** Feature-selection options worth adding:
   - number of HVGs: 1,000 / 2,000 / 5,000;
   - TF filter (number of targets detected), or unscaled ULM scores;
   - CCC: top N pairs, or collapse pairs by ligand;
   - number of PCs per layer;
   - extra-layer weights 0.25 / 0.5 / 1 / 2.
3. **Ask the colleague** how the CCC scores were made: method, kernel bandwidth, pair filtering.
4. **Optional:** the P9 remote-zone slide (ACH0023). It has CCC scores but no RNA or TF runs yet.
