# Handoff, 24 November 2026 (evening): the 23 Nov panels with joint, integrated clustering

Track C, branch `hack/c-derived-modalities`. Paths are relative to `Derived_Modalities/`. Builds on `HANDOFF_2026-11-24_fig1hij_spatialleiden.md` (joint pipeline, other VM) and `HANDOFF_2026-11-23_resume_on_new_vm.md` (per-slide work).

## 1. Goal

Redraw the 23 Nov multipanels for **P1 control (10X001), P9 IZ (ACH0012) and P9 IZ rep2 (ACH0013)**, with domains from **one joint SpatialLeiden clustering of all 27 Visium slides** on Harmony-integrated layers, instead of clustering each slide on its own. The two panels:

- RNA / TF / RNA + TF;
- RNA / CCC / RNA + CCC.

The authors' Extended Data Fig. 3 (composition UMAP by patient and region) is the planned sanity check for the integration.

## 2. Strategy

This follows the preprint's multi-sample recipe (Methods: integrate first, build the kNN graph on the corrected space, one spatial graph per sample) and its `clustering/spatialleiden_multisample.py`.

| Step | Per slide (23 Nov) | Joint (this run) |
|---|---|---|
| Spots | one slide (~4,300) | all 27 slides (88,704) |
| RNA embedding | slide PCA 30 | shared genes, HVG 2,000 (batch-aware), PCA 30, **Harmony** (patient + sample + Visium batch) |
| TF layer | CollecTRI ULM, z-score per slide, PCA 30 | CollecTRI ULM on 27 slides, z-score pooled, PCA 30, Harmony |
| CCC layer | colleague's LIANA+ scores (4 slides) | **our LIANA+ reconstruction** on 27 slides (section 3); 286 pairs scored on all slides; z-score pooled, PCA 30, Harmony |
| Expression graphs | 15-NN within the slide | 15-NN on the Harmony embedding (neighbours can be on other slides) |
| Spatial graph | Visium grid, 6 neighbours | the same, built per slide (`library_key="sample"`), so no edge crosses slides |
| Clustering | `spatialleiden` / `spatialleiden_multimodal` | the same functions; the MuData holds the integrated layers |
| Weights | RNA 1, second layer 0.5, spatial 0 / 0.25 / 0.75 | the same, plus spatial 1.2 (the preprint's value for integrated multi-sample data) |
| K | authors' niche count on the slide (8 / 8 / 9) | 9, shared by all 27 slides |
| Domain identity | per slide | shared across slides: a colour is the same domain on every slide |

**How to read the joint panels:**

- **The same spatial weight pulls less jointly.** Part of each spot's expression neighbourhood lies on other slides, so 0.75 jointly is weaker than 0.75 per slide. That is why the 1.2 run is included.
- **Joint domains are coarser on a single slide.** With 9 domains spread over 27 slides, uniform P1 may fall into one or two domains, and fine structures (the P9 rep2 artery band) may merge into larger ones. The figure tables give the ARI against the 23 Nov per-slide clusters to quantify this.

## 3. Our CCC for all 27 slides (`02_modalities/liana/liana_local_all27.py`)

The colleague's scores exist for 4 slides only, and their settings are not documented.

**Calibration on P9 IZ.** We ran LIANA+ 1.8.1 `li.mt.bivariate` with the local cosine score, a Gaussian kernel, cutoff 0.1 and the h5ad's log-normalised X, varying the bandwidth (in multiples of the spot spacing):

| Bandwidth (× spot spacing) | 0.4 | 0.6 | **0.75** | 0.855 | 1 | 2 | 3 |
|---|---|---|---|---|---|---|---|
| Median per-pair r with their scores | 0.20 | 0.77 | **0.84** | 0.83 | 0.82 | 0.67 | 0.58 |

**The chosen bandwidth, 0.75, is reproduced on all 4 shared slides:** median per-pair r = 0.79 (P1), 0.84 (P9 IZ), 0.81 (P9 rep2), 0.79 (P9 RZ). See `results/kuppe_2022/ccc_all27_ours/compare_colleague.csv`.

**Treat it as a reconstruction, not the colleague's layer.** Some setting still differs, and we score more pairs per slide (473–1,298, against their 154–620).

**Outputs:**

- per-slide scores: `data/kuppe_2022/modalities/ccc_all27_ours/ccc_<run>.csv`;
- run summary: `results/kuppe_2022/ccc_all27_ours/liana_qc.csv`.

## 4. Environment

`spatialleiden_joint` is a clone of `spatialleiden` with three pins:

- **`pynndescent=0.5.13`:** 0.7 returns `csr_array` kNN graphs that scanpy 1.10 cannot read on large data (already in `spatialleiden.yml`).
- **`setuptools<81`:** also already in `spatialleiden.yml`.
- **`umap-learn=0.5.12`:** **new trap.** Installing pynndescent 0.5.13 upgraded umap-learn to 0.6.0, which fails in `sc.pp.neighbors` with a numba error ("No conversion from array(float32, 2d, A) to array(float32, 2d, C) for knn_dists"). This pin should be added to `spatialleiden.yml`.

The spec also adds `datashader==0.19.1`.

## 5. Commands (quick mode: seed 0, spatial weights 0, 0.25, 0.75, 1.2)

Run from `Derived_Modalities/` with env `spatialleiden_joint` (CollecTRI and LIANA+ with an env that has decoupler 2.2 and liana 1.8).

```
# layers
python 02_modalities/collectri/collectri_tf_activity_all.py data/kuppe_2022 results/kuppe_2022/collectri_all27
python 02_modalities/liana/liana_local_all27.py data/kuppe_2022 results/kuppe_2022/ccc_all27_ours 6
# integration
python 07_integration/exp1_harmony_covariates/exp1_covariates.py data/kuppe_2022 results/kuppe_2022/exp1_harmony_covariates
python 07_integration/tf_harmony.py data/kuppe_2022 results/kuppe_2022/tf_harmony
python 07_integration/ccc_harmony.py harmony data/kuppe_2022 data/kuppe_2022/modalities/ccc_all27_ours results/kuppe_2022/ccc_harmony_ours
# joint SpatialLeiden (each with EXP_SEEDS=0 EXP_RATIOS=0,0.25,0.75,1.2; drop them for the full 5 seeds x 8 weights)
S=03_single_modality_clustering/exp2_joint_harmony
python $S/exp2_joint_sweep.py data/kuppe_2022 results/kuppe_2022/exp1_harmony_covariates/embeddings.npz results/kuppe_2022/exp2_joint_harmony authors 6
python $S/exp2_joint_sweep.py data/kuppe_2022 results/kuppe_2022/tf_harmony/tf_embeddings.npz results/kuppe_2022/exp2_joint_harmony_tf X_harmony 6
python $S/exp2_joint_sweep.py data/kuppe_2022 results/kuppe_2022/ccc_harmony_ours/ccc_embeddings.npz results/kuppe_2022/exp2_joint_harmony_ccc_ours X_harmony 6
python $S/exp3_joint_multimodal_sweep.py data/kuppe_2022 results/kuppe_2022/exp1_harmony_covariates/embeddings.npz:authors \
    tf=results/kuppe_2022/tf_harmony/tf_embeddings.npz:X_harmony results/kuppe_2022/exp3_joint_rna_tf 0.5 6
python $S/exp3_joint_multimodal_sweep.py data/kuppe_2022 results/kuppe_2022/exp1_harmony_covariates/embeddings.npz:authors \
    ccc=results/kuppe_2022/ccc_harmony_ours/ccc_embeddings.npz:X_harmony results/kuppe_2022/exp3_joint_rna_ccc_ours 0.5 6
# panels
python $S/fig_joint_three_slides.py data/kuppe_2022 results/kuppe_2022 rna_tf
python $S/fig_joint_three_slides.py data/kuppe_2022 results/kuppe_2022 ccc
```

**Runner:** on this VM all of it is chained in `results/kuppe_2022/run_quick.sh` (screen `joint_quick`). The exit code is in `QUICK_EXIT`, and the timeline in `run_quick.log`.

**Speed:** one Leiden run on the joint graph takes about 30 s, and holding K needs about 10–20 runs per setting. The full sweeps (5 seeds × 8–9 weights) take roughly 3 h for all five layers. Quick mode takes about 30 min.

## 6. Panels (`03_single_modality_clustering/exp2_joint_harmony/fig_joint_three_slides.py`)

| Panel | Columns |
|---|---|
| `rna_tf` | H&E · RNA only · RNA + space · TF only · TF + space · RNA + TF + space · authors' niches · SPP1 · NPPB |
| `ccc` | H&E · RNA only · RNA + space · CCC only · CCC + space · RNA + CCC + space · authors' niches · SPP1 → CD44 (our CCC) |

**Settings:**

- **Spatial weight:** 0.75 for RNA and TF, 0.25 for CCC, as on 23 Nov. Change with `FIG_SPACE_W` / `FIG_CCC_W`, for example `FIG_SPACE_W=1.2`.
- **Colours:** each run's domains are matched over all 88,704 spots to RNA + space, and coloured by domain size.

**Outputs** (in `results/kuppe_2022/joint_three_slides/`): `fig_joint_<panel>.png/.pdf` and `_settings.csv`. The table lists, per slide and column: domains present, coherence, ARI against the 23 Nov per-slide clustering of the same layers, and ARI against the authors' cell-type niches.

## 7. State at handoff (24 Nov, ~17:30 CET)

| Step | Status |
|---|---|
| CCC on 27 slides (ours) | done |
| CollecTRI on 27 slides | done; identical to the 23 Nov files on the 3 shared slides |
| RNA Harmony embedding (`embeddings.npz`, keys none / batch / authors) | saved |
| TF Harmony embedding (`tf_embeddings.npz`) | saved |
| CCC Harmony embedding (`ccc_harmony_ours/`, 286 pairs) | saved |
| Integration diagnostics of experiment 1 (LISI, silhouette, UMAPs, `RESULT.md` numbers) and of `tf_harmony.py` | **stopped for speed**, embeddings kept. Rerun both scripts in full for the sanity check (section 8) |
| Five joint sweeps (quick mode) | running, from 17:21 CET; panels chained after them; expected ~17:45–18:00 CET |

**harmonypy is not bit-reproducible.** These embeddings differ slightly from the other VM's. Do not mix label files from the two VMs.

## 8. Next steps

1. **Check the quick results:**
   - every run reached K = 9 (`exp2_sweep.csv` / `exp3_sweep.csv`: `n_clusters`, search flags);
   - the domains per slide and the ARI against 23 Nov in the panel tables.
2. **Integration sanity check against the authors' Extended Data Fig. 3:**
   - reproduce ED Fig. 3a–d as published (`01_data_preprocess/kuppe_ed_fig3/ed_fig3abc.py`, `ed_fig3d.py`; composition UMAP, no integration involved);
   - draw the same patient / region UMAPs on the RNA, TF and CCC Harmony embeddings (`common/plot_patient_region_umap.py`);
   - rerun the experiment 1 and `tf_harmony.py` diagnostics in full (iLISI on patient, cLISI on zone).

   **What good integration looks like:** patients mix, while regions stay as separated as on the composition UMAP.
3. **Full sweeps overnight:** 5 seeds × all weights, for seed stability.
4. **Optional:**
   - an scVI embedding for RNA as a robustness check (the preprint's preferred integration);
   - `umap-learn=0.5.12` into `spatialleiden.yml`;
   - the colleague's own 27-slide CCC, if it becomes available, to replace our reconstruction.
5. **Commit:**
   - code: `02_modalities/liana/liana_local_all27.py`, `03_.../exp2_joint_harmony/fig_joint_three_slides.py`, this handoff;
   - small results with `git add -f`: label CSVs, sweep CSVs, panel PNGs.

   Nothing from this evening is committed yet.
