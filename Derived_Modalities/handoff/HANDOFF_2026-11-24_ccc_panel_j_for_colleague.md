# Handoff, 24 November 2026: panel j on CCC, for running on another VM

Track C, branch `hack/c-derived-modalities`. Paths are relative to `Derived_Modalities/`. Background: `HANDOFF_2026-11-24_fig1hij_spatialleiden.md`.

## 0. Deadline version: about 45 min, run in this order

The full sweeps (5 seeds × 8 weights) take 45–60 min each. For today's deadline, run **quick mode**: seed 0 and only the four plotted weights. Each sweep then takes about 20–30 min.

1. **Unpack** `colleague_bundle_ccc_2026-11-24.tar.gz` inside `Derived_Modalities/`. It holds the scripts, the RNA Harmony embedding and the PROGENy scores.
2. **Put the data in place:** `data/kuppe_2022/` with `visium_tables/` and `processed1/` (section 2).
3. **CCC embedding** (about 3 min):

   ```
   python 07_integration/ccc_harmony.py harmony data/kuppe_2022 <your_non_integrated_ccc> results/kuppe_2022/ccc_harmony
   ```

   `<your_non_integrated_ccc>` is a folder of `ccc_<run>.csv` or one `.h5ad`; see section 3.
4. **Start both sweeps at the same time,** each with half your cores (`<n>` = cores / 2):

   ```
   EXP_SEEDS=0 EXP_RATIOS=0,0.5,1.0,1.5 python -W ignore 03_single_modality_clustering/exp2_joint_harmony/exp2_joint_sweep.py \
       data/kuppe_2022 results/kuppe_2022/ccc_harmony/ccc_embeddings.npz results/kuppe_2022/exp2_joint_harmony_ccc X_harmony <n> &
   EXP_SEEDS=0 EXP_RATIOS=0,0.5,1.0,1.5 python -W ignore 03_single_modality_clustering/exp2_joint_harmony/exp3_joint_multimodal_sweep.py \
       data/kuppe_2022 results/kuppe_2022/exp1_harmony_covariates/embeddings.npz:authors \
       ccc=results/kuppe_2022/ccc_harmony/ccc_embeddings.npz:X_harmony \
       results/kuppe_2022/exp3_joint_multimodal_rna_ccc 0.5 <n> &
   ```

   Quick mode leaves only 3 searches at once per sweep, so a few cores each is enough.
5. **Figures,** once each sweep has written `exp2_labels.csv.gz`: the two commands in section 6.

If something fails, send the `exp2.log` / `exp3.log` and the traceback.

## 1. What to make

We are rebuilding Kuppe et al. 2022 Fig. 1h–j:

- **Panels h and i** reproduce the paper (cell2location proportions; PROGENy TGFβ and NFκB). They are fixed and identical in every version.
- **Panel j** shows SpatialLeiden domains at increasing spatial weight.

Each version of panel j:

1. **Joint clustering:** one clustering of all 27 Visium slides (88,704 spots) at **K = 9**, the authors' cell-type niches.
2. **Rows:** the four Fig. 1 slides.

   | Row | Visium run |
   |---|---|
   | Ctrl P7 | ACH003 |
   | BZ P3 | 10X0026 |
   | IZ P3 | 10X0017 |
   | FZ P14 | ACH005 |

3. **Columns:** spatial weight 0, 0.5, 1 and 1.5. Weight 0 is expression only.

| Version | Panel j input | Who |
|---|---|---|
| 1 | RNA (Harmony) | done on our VM |
| 2 | TF activity (CollecTRI, Harmony) | done on our VM |
| **3** | **CCC alone** | **you** |
| **4** | **Multimodal: RNA + CCC + spatial graph** | **you** |

We ran out of compute, so versions 3 and 4 run on your VM.

## 2. What you need on your VM

1. **Code:** this branch. The scripts below must be committed and pushed first; if they are not on the branch yet, they are also in the bundle (step 3).
2. **Data:** `data/kuppe_2022/`, from the OneDrive folder `SpatialLeiden_hackathon/derived_modalities_2026-11-23/`, `bundles/` (see `HANDOFF_2026-11-24_start_here.md`). You need:
   - `processed1/`, the 27 Visium h5ads;
   - `visium_tables/`, the spot tables `<run>_spots.csv`;
   - `visium_raw/metadata-Visium.csv`.
3. **Our bundle** (`results/` and `data/` are not in git). It contains:

   | File | Used for |
   |---|---|
   | `results/kuppe_2022/exp1_harmony_covariates/embeddings.npz` | RNA Harmony embedding (key `authors`); needed for version 4 |
   | `data/kuppe_2022/modalities/progeny/` | PROGENy scores for the four slides; needed for panel i |

   Unpack the bundle inside `Derived_Modalities/`.
4. **Your CCC for all 27 slides.** Section 3 says which version to use.
5. **Env:** `spatialleiden` (from `../spatialleiden.yml`), with `mudata` and `harmonypy` 0.2.0. Check with `python -c "import spatialleiden, mudata, harmonypy, squidpy"`.

## 3. Step 1: CCC embedding

Script: `07_integration/ccc_harmony.py`. It builds the CCC layer the same way as the TF layer (`07_integration/RESULT_tf_harmony.md`), so all layers are comparable:

1. Keep the pairs present on all 27 slides.
2. z-score each pair over all spots pooled.
3. PCA with 30 components (seed 0).
4. Harmony on the PCs with covariates patient + sample + visium_batch (seed 0). These are the covariates used for RNA and TF.

**Use your non-integrated CCC (recommended).** The script then integrates it with Harmony exactly like RNA and TF:

```
python 07_integration/ccc_harmony.py harmony data/kuppe_2022 <your_ccc> results/kuppe_2022/ccc_harmony
```

**Or your integrated CCC.** Use this only if you prefer your own integration. The script then runs PCA only, with no second correction, and the embedding key becomes `X_pca`:

```
python 07_integration/ccc_harmony.py pca data/kuppe_2022 <your_integrated_ccc> results/kuppe_2022/ccc_integrated
```

**Input format for `<your_ccc>`:**

- **Either** a folder of per-slide tables `ccc_<run>.csv`: index = spot barcode, columns `LIGAND^RECEPTOR`.
- **Or** one `.h5ad` (spots × pairs) with an obs column naming the run (`--sample-col`, default `hca_sample_id`). obs_names are barcodes or `<run>_<barcode>`.
- `<run>` is the Visium run id (`10X001`, `ACH0012`, …; column `hca_sample_id` in `metadata-Visium.csv`).
- Every spot of `visium_tables/<run>_spots.csv` must have scores.

**Checks:**

- The script prints how many pairs are shared by all 27 slides. It stops below 50; if that happens, relax your per-slide pair filter.
- About 3 min on our VM (88,704 spots; Harmony takes most of it).

Output: `ccc_embeddings.npz` (obs_names `<run>_<barcode>`, `X_pca`, plus `X_harmony` in harmony mode) and `ccc_pairs.txt`.

## 4. Step 2: version 3, SpatialLeiden on CCC alone

Script: `03_single_modality_clustering/exp2_joint_harmony/exp2_joint_sweep.py`. It is the same one we used for RNA and TF.

```
python -W ignore 03_single_modality_clustering/exp2_joint_harmony/exp2_joint_sweep.py \
    data/kuppe_2022 results/kuppe_2022/ccc_harmony/ccc_embeddings.npz \
    results/kuppe_2022/exp2_joint_harmony_ccc X_harmony 30
```

The last two arguments are the embedding key (`X_pca` if you used pca mode) and the number of worker processes. Use about the number of cores, minus a few.

**What it runs:**

- **Graphs:** 15-NN on the embedding, plus a squidpy hexagonal grid built per slide (no spatial edge crosses slides).
- **Weights:** `layer_ratio` 0, 0.25, 0.5, 0.75, 1, 1.5, 2, 2.5; seeds 0–4; K = 9.
- **Holding K:**
  - At ratio 0, the expression resolution is searched with `spatialleiden(layer_ratio=0)` itself.
  - At ratio > 0, the spatial resolution is searched.
  - This follows experiment 0; it differs from the preprint's `search_resolution` on purpose (see `HANDOFF_2026-11-24_fig1hij_spatialleiden.md`).

**Runtime on our VM** (32 cores):

- One Leiden run on the joint graph takes about 30 s.
- The 5 ratio-0 searches take 4–8 min.
- The 35 searches at ratio > 0 take 5–20 min each, in parallel.
- In total, about 45–60 min with 30 workers.

**Outputs:**

| File | Content |
|---|---|
| `exp2_labels.csv.gz` | Spots × `r<ratio>_s<seed>` |
| `exp2_sweep.csv` | Per run: K, resolutions, search outcome, coherence gap, ARIs to ratio 0, to the authors' niches and to sample |
| `exp2_per_slide.csv` | The same per slide |
| `exp2_seed_stability.csv` | Pairwise ARI among seeds per ratio |
| `exp2.log` | Run log |

**Delete `exp2.log` before a rerun**; it is opened in append mode.

**Note on CCC:** the CCC scores are computed over a spatial kernel, so they already contain space. Adding the spatial graph on top counts space twice; yesterday we used spatial weight 0.25 for CCC for this reason. We screen the same weights as for RNA and TF so the versions are comparable, but expect CCC to smooth sooner.

## 5. Step 3: version 4, multimodal RNA + CCC + spatial

Script: `03_single_modality_clustering/exp2_joint_harmony/exp3_joint_multimodal_sweep.py`.

```
python -W ignore 03_single_modality_clustering/exp2_joint_harmony/exp3_joint_multimodal_sweep.py \
    data/kuppe_2022 \
    results/kuppe_2022/exp1_harmony_covariates/embeddings.npz:authors \
    ccc=results/kuppe_2022/ccc_harmony/ccc_embeddings.npz:X_harmony \
    results/kuppe_2022/exp3_joint_multimodal_rna_ccc 0.5 30
```

**What it runs:**

- **Layers** (`spatialleiden.spatialleiden_multimodal`): `rna` (15-NN on the RNA Harmony embedding), `ccc` (15-NN on the CCC embedding) and `spatial` (per-slide grid).
- **Weights:** RNA 1, CCC 0.5 (the 5th argument; the `rna_ccc` setting of `04_multimodal_weight_grid/proto_rna_tf/mm_configs.py`). The spatial weight is swept over the same values as the layer ratio in step 2.
- **Holding K:**
  - At spatial weight 0, the shared expression resolution (RNA and CCC) is searched with the multimodal call itself.
  - At weight > 0, the spatial resolution is searched.
  - If that cannot reach K, the expression resolution is searched instead (column `expression_fallback`), as in `mm_configs.py`.
- **Outputs:** `exp2_labels.csv.gz` (same name and format as step 2, so the figure script reads it unchanged), `exp3_sweep.csv`, `exp3_per_slide.csv`, `exp3_seed_stability.csv`, `exp3.log`.
- **Runtime:** expect it to be slower than step 2 (two expression graphs).

**Quick prototype first (recommended).** Seed 0 and only the four plotted weights, about 20–30 min:

```
EXP_SEEDS=0 EXP_RATIOS=0,0.5,1.0,1.5 python -W ignore 03_single_modality_clustering/exp2_joint_harmony/exp3_joint_multimodal_sweep.py \
    data/kuppe_2022 \
    results/kuppe_2022/exp1_harmony_covariates/embeddings.npz:authors \
    ccc=results/kuppe_2022/ccc_harmony/ccc_embeddings.npz:X_harmony \
    results/kuppe_2022/exp3_proto_rna_ccc 0.5 4
```

## 6. Step 4: figures

Script: `01_data_preprocess/kuppe_fig1/fig1hij.py`. Panel j shows seed 0 at weights 0 / 0.5 / 1 / 1.5. To change them, set `FIG_SL_RATIOS=0,0.5,1.0,2.0`; the labels file needs those columns.

```
# version 3: CCC alone
python 01_data_preprocess/kuppe_fig1/fig1hij.py data/kuppe_2022 figures/kuppe_fig1h \
    results/kuppe_2022/exp2_joint_harmony_ccc/exp2_labels.csv.gz kuppe_fig1hij_ccc CCC
# version 4: RNA + CCC + spatial
python 01_data_preprocess/kuppe_fig1/fig1hij.py data/kuppe_2022 figures/kuppe_fig1h \
    results/kuppe_2022/exp3_joint_multimodal_rna_ccc/exp2_labels.csv.gz kuppe_fig1hij_rna_ccc "RNA + CCC (0.5)"
```

- **Domain numbering:** domains are numbered D1–D9 by size at weight 0. At the other weights, each domain takes the number of the weight-0 domain it overlaps most (Hungarian matching over all spots), so a colour means the same domain in every row and column.
- **Panel j title:** for version 4 it reads "SpatialLeiden domains (RNA + CCC (0.5), Harmony, all 27 slides jointly; K = 9)".
- **Panel j columns:** labelled "ratio …", which is the spatial weight in version 4.

## 7. Checks before sending results back

1. **K reached:** every run has `n_clusters == 9` and `search_failed == False` (`exp2_sweep.csv` / `exp3_sweep.csv`). For version 4, also count how many runs used `expression_fallback`.
2. **Seed stability:** `*_seed_stability.csv`. If stability is low, panel j's seed 0 may not be representative.
3. **Domains per slide:** `*_per_slide.csv` gives the number of domains present on each slide. With 9 domains shared across 27 slides, a slide can show only a few.

**Send back:** the output folders `results/kuppe_2022/exp2_joint_harmony_ccc/` and `.../exp3_joint_multimodal_rna_ccc/`, plus the figures `figures/kuppe_fig1h/kuppe_fig1hij_ccc.*` and `kuppe_fig1hij_rna_ccc.*`.

## 8. Testing status of this code

| Script | Status |
|---|---|
| `exp2_joint_sweep.py` | In production on our VM (RNA and TF runs) |
| `ccc_harmony.py` | Smoke-tested in both modes on all 27 slides, with CollecTRI TF scores standing in for CCC (h5ad input) |
| `exp3_joint_multimodal_sweep.py` | Prototype running on our VM now with RNA + TF (same code path, `tf=` instead of `ccc=`). Check that run's log before relying on it |
| `fig1hij.py` | Used for the RNA draft |
