# Handoff, 24 November 2026: Kuppe Fig. 1h–j with SpatialLeiden as panel j

Track C, branch `hack/c-derived-modalities`. Paths are relative to `Derived_Modalities/`. Continues `HANDOFF_2026-11-24_start_here.md`.

## 1. Goal

Hackathon deliverable: reproduce Kuppe et al. 2022 Fig. 1h (and i), then add columns that show what SpatialLeiden does.

- **Panels h and i:** reproduced as published.
- **Panel j:** the paper's TF binding panel is replaced by SpatialLeiden domains at increasing spatial weight.
- **Two versions of the figure:** j on RNA, and j on TF activity. h and i are identical in both.

**Rows** (the paper's four slides):

| Row | Slide | Visium run |
|---|---|---|
| Ctrl tissue (P7) | control_P7 | ACH003 |
| BZ tissue (P3) | RZ/BZ_P3 | 10X0026 |
| IZ tissue (P3) | IZ_P3 | 10X0017 |
| FZ tissue (P14) | FZ_P14 | ACH005 |

## 2. Done

### Panel h: cell2location, reproduced

- **Source:** the authors' `c2l_props` values (`prop_<celltype>` in `data/kuppe_2022/visium_tables/<run>_spots.csv`).
- **Proportions, not absolute abundance (`c2l`).** We checked both against the published panel. Abundance does not reproduce the IZ P3 and FZ P14 vCM patterns; proportions do. The authors' repo (saezlab/visium_heart) has no script for Fig. 1h–j; their generic c2l plot script, `st_snRNAseq/04_deconvolution/plot_ct_scores.R`, is plain `SpatialFeaturePlot`.
- **Colour scale:** per map, clipped at the 1st and 99th percentile (Seurat `q1` / `q99`, as in their `st_snRNAseq/07_structure/plot_spatial_ex_panels.R`), inferno palette.
- **Matched by eye:** black background, dashed outline, arrows.
- **Status:** all 12 maps match the published panel.

### Panel i: PROGENy, recomputed with the authors' formula

- **Formula:** the top-500 human PROGENy model × log-normalised RNA, then a z-score per pathway per slide (their `st_snRNAseq/utils/funcomics.R`, `add_path_activities`).
- **Why recomputed:** their own `progeny` assay is in Seurat `.rds` objects that are not on this VM.
- **Code:** `02_modalities/progeny/progeny_kuppe.py` (decoupler env). Outputs in `data/kuppe_2022/modalities/progeny/` (not in git).
- **Status:** TGFβ and NFκB match the published panel closely.

### Dropped: CollecTRI as panel j

- **What it was:** a first version of j used CollecTRI TF activity in place of the paper's TF binding (HINT footprinting on snATAC).
- **How it looked:** SRF in P7 and FOS–JUNB in P14 agreed with the paper; MEF2C was more diffuse.
- **Why it went:** the user replaced j with SpatialLeiden.
- **Leftover files:** the CollecTRI-map version is still on disk as `figures/kuppe_fig1h/kuppe_fig1hij.{pdf,png}` and `fig1hij_summary.json`. **Delete them** once the new versions are rendered.

## 3. Running at handoff time: experiment 2, joint SpatialLeiden sweep

**Script:** `03_single_modality_clustering/exp2_joint_harmony/exp2_joint_sweep.py`. Both runs were started in the background on this VM (32 cores).

| Run | Input embedding | Output | Workers |
|---|---|---|---|
| RNA | `results/kuppe_2022/exp1_harmony_covariates/embeddings.npz`, key `authors` | `results/kuppe_2022/exp2_joint_harmony/` | 30 |
| TF | `results/kuppe_2022/tf_harmony/tf_embeddings.npz`, key `X_harmony` | `results/kuppe_2022/exp2_joint_harmony_tf/` | 14 |

**Design** (decided with the user, 24 Nov):

- **One joint clustering of all 27 slides (88,704 spots).** Figures show the four Fig. 1 slides per patient, but the domains come from the multi-sample run, so a colour means the same domain on every slide.
- **Expression graph:** 15-NN (scanpy, seed 0) on the 30 Harmony dimensions. Covariates are patient + sample + batch, as fixed in `07_integration/exp1_harmony_covariates/RESULT.md` and `07_integration/RESULT_tf_harmony.md`.
- **Spatial graph:** squidpy grid, 6 neighbours, built per slide (`library_key="sample"`), so no spatial edge crosses slides. This follows the preprint's `clustering/spatialleiden_multisample.py`.
- **K = 9:** the authors' cell-type niches, pooled over 27 slides; also experiment 1's K. Yesterday's per-slide rule ("niches present on the slide") would give 4 on IZ P3, which is 99 % niche 2, so it was not used.
- **Settings and K search:** as in experiment 0 (`03_single_modality_clustering/exp0_layer_ratio_sweep/README.md`).
  - Search order: at ratio 0, the expression resolution is searched with `spatialleiden(layer_ratio=0)` itself, not plain Leiden (see exp0's sanity-check finding). At ratio > 0, the spatial resolution is searched.
  - Parameters: `directed=(False, False)`, `n_iterations=-1`, up to 30 search iterations.
  - Ratios 0, 0.25, 0.5, 0.75, 1, 1.5, 2, 2.5; seeds 0–4.
- **Runtime:** one Leiden run on the joint graph takes about 30 s. The searches are sequential, so expect roughly 15–30 min per run, and longer while the two runs share cores.

**Outputs per run:**

| File | Content |
|---|---|
| `exp2_labels.csv.gz` | Spots × `r<ratio>_s<seed>`; index `<sample>_<barcode>` |
| `exp2_sweep.csv` | Per (ratio, seed): K, resolutions, search outcome, pooled coherence gap, ARI to ratio 0, to the cell-type and molecular niches and to sample, cluster sizes |
| `exp2_per_slide.csv` | The same per slide, plus the number of domains present on each slide |
| `exp2_seed_stability.csv` | Pairwise ARI among seeds per ratio |
| `exp2.log` | Run log |

**If a run died** (no `exp2_sweep.csv` in its folder): rerun it from `Derived_Modalities/` with the `spatialleiden` env (`/home/gage11/envs/spatialleiden` on this VM):

```
python -W ignore 03_single_modality_clustering/exp2_joint_harmony/exp2_joint_sweep.py \
    data/kuppe_2022 results/kuppe_2022/exp1_harmony_covariates/embeddings.npz \
    results/kuppe_2022/exp2_joint_harmony authors 30
python -W ignore 03_single_modality_clustering/exp2_joint_harmony/exp2_joint_sweep.py \
    data/kuppe_2022 results/kuppe_2022/tf_harmony/tf_embeddings.npz \
    results/kuppe_2022/exp2_joint_harmony_tf X_harmony 14
```

`exp2.log` is opened in append mode, so delete it before a rerun.

## 4. Next steps

1. **Check each run before plotting.**
   - Every (ratio, seed) has `n_clusters == 9` and `search_failed == False` in `exp2_sweep.csv`.
   - Look at seed stability.
   - Look at how many domains appear on each of the four slides (`exp2_per_slide.csv`). With K = 9 shared across 27 slides, a slide can show only a few domains.
2. **Render both versions** (spatialleiden env). h and i are identical in both. j shows ratios 0 / 0.5 / 1 / 2 for seed 0.

   ```
   python 01_data_preprocess/kuppe_fig1/fig1hij.py data/kuppe_2022 figures/kuppe_fig1h \
       results/kuppe_2022/exp2_joint_harmony/exp2_labels.csv.gz kuppe_fig1hij_rna RNA
   python 01_data_preprocess/kuppe_fig1/fig1hij.py data/kuppe_2022 figures/kuppe_fig1h \
       results/kuppe_2022/exp2_joint_harmony_tf/exp2_labels.csv.gz kuppe_fig1hij_tf "TF activity"
   ```

   Domains are numbered D1–D9 by size at ratio 0. At the other ratios, each domain takes the number of the ratio-0 domain it overlaps most (Hungarian matching over all spots). The script also rewrites `kuppe_fig1h.{pdf,png}` (h alone).
3. **Delete the stale CollecTRI-map version** (`kuppe_fig1hij.{pdf,png}`, `fig1hij_summary.json`).
4. **Write up the results:**
   - Add a `README.md` (design, which is section 3 above) and a `RESULT.md` (numbers per ratio, RNA against TF) to `03_single_modality_clustering/exp2_joint_harmony/`.
   - Update `figures/kuppe_fig1h/README.md`, which still describes j as CollecTRI maps, and its row in `figures/README.md`.
5. **Commit, only when the user asks.** Nothing from this block is committed yet.
   - New: `01_data_preprocess/kuppe_fig1/`, `02_modalities/progeny/`, `03_single_modality_clustering/exp2_joint_harmony/`, `figures/kuppe_fig1h/` (add the PDFs with `git add -f`), this handoff.
   - Not ours: `git status` also lists modified and untracked files from other work (for example `exp0_sweep.py`, `mm_configs.py`, `spatialleiden.yml`, `02_modalities/liana/`). Leave them alone unless their owner says otherwise.

## 5. Caveats to keep straight on slides

- **Panel i is recomputed, not the authors' stored assay.** It uses the same formula and the same normalisation (the Seurat RNA data slot = `processed1` `X`), so differences should be negligible.
- **Ratio 0 is not plain Leiden.** A zero-weight spatial layer still changes leidenalg's search path (experiment 0 finding), so the ratio-0 baseline is `spatialleiden(layer_ratio=0)`.
- **K = 9 is shared across 27 slides,** so a single slide may show only a subset of domains, for example IZ P3 being mostly one domain.
- **harmonypy is not bit-reproducible,** so always use the saved embeddings, never a recomputation.
- **Cosmetic:** Arial is not installed on this VM (DejaVu Sans fallback), and arrow positions in h and i are approximate.
