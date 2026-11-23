# Handoff: Kuppe 2022 cohort structure, and P1 + P9 as the first SpatialLeiden test

23 November 2026. Paths are relative to `Derived_Modalities/`. Companion figures:

- `figures/kuppe_cohort/` (cohort overview),
- `figures/kuppe_fig2b/` (reproduction of the paper's Fig. 2b),
- `results/kuppe_2022/slide_atlas/kuppe_visium_slide_atlas.pdf` (one page per Visium run; local only, 42 MB).

## 1. Cohort structure

**Who.** There are 23 patients:

- 4 control donors: non-transplanted donor hearts.
- 7 acute myocardial infarction (MI) patients, sampled 2–11 days after symptom onset.
- 9 transplant recipients at later stages, 31–166 days after onset.
- 3 validation patients with snRNA-seq only (P21–P23).

Together they give 31 samples. Source: the paper's Supplementary Table 1, cross-checked against the Zenodo sample sheets with no mismatches.

**Two ways to group the samples.** Both are available for every Visium section.

| | Zone (where the sample was taken) | Pathologist group (histology of the whole slide) |
|---|---|---|
| Column | `major_labl` | `patient_group` |
| Levels | control, remote (RZ), border (BZ), ischaemic (IZ), fibrotic (FZ) | `group_1` myogenic, `group_2` ischaemic, `group_3` fibrotic |
| Visium sections (27) | 4 / 5 / 3 / 9 / 6 | 13 / 9 / 5 |

**The three pathologist groups.**

| Group | Tissue | Sections |
|---|---|---|
| **Myogenic** (13 sections, 14 runs) | Mostly intact muscle | 4 controls (P1, P7, P8, P17), 5 remote (RZ_P6 has 2 runs), 3 border, and RZ/FZ_P5 |
| **Ischaemic** (9) | Infarct core with necrosis, inflammation and early granulation tissue | IZ_P3, IZ_P10, IZ_P15, IZ_P16, IZ/BZ_P2, GT/IZ_P9, GT/IZ_P9_rep2, GT/IZ_P13, GT/IZ_P15 |
| **Fibrotic** (5) | Scar | FZ_P14, FZ_P18, FZ_P20, FZ/GT_P4, FZ/GT_P19 |

**Things that constrain any analysis.**

- **Controls sit inside the myogenic group.** There is no separate control group at slide level. To compare with healthy tissue, use the zone, or split myogenic into control and the rest.
- **Group is confounded with time.** Ischaemic slides are acute (2–11 days; P13 at 45 days), fibrotic slides are late (40–166 days).
- **Zone is confounded with patient.** Most patients give one zone. Only P2, P3 and P9 contribute several zones, and P15 gives two ischaemic samples.
- **Labels are per slide, not per spot.** No spot-level zone or histology annotation ships with the data.
- **Two Visium processing batches:** ACH (20 runs) and 10X (8 runs).
- **The published niches are not spatial.** The authors' cell-type niches (clustered from cell2location proportions) and molecular niches (clustered from expression) were made without spot coordinates. They are a reference, not ground truth.

## 2. Proposed first implementation: P1 (control) and P9 (ischaemic), the slides of the paper's Fig. 2b

| | P1 (control) | P9 (IZ) |
|---|---|---|
| Section and run | `control_P1`, `10X001` | `GT/IZ_P9_rep2`, `ACH0013` |
| Group, days after onset | Myogenic, control donor | Ischaemic, 2 days |
| Spots | 4,269 | 4,113 |
| Median UMIs / genes per spot | 3,264 / 1,649 (12th deepest of 28 runs) | 8,582 / 3,458 (**2nd deepest of 28**) |
| Median % mitochondrial UMIs | 32 % | 6 % |
| Cell-type niches present | 8; niche 8 alone covers 44 % of spots | 9; niche 4 covers 24 %, niche 5 14 % |
| Cell-type niche coherence (shuffled)\* | 0.32 (0.26): close to chance | 0.45 (0.16): clearly structured |
| Matched snRNA-seq / snATAC-seq nuclei | 11,416 / 3,183 | 3,166 / 1,469 |
| Visium batch | 10X | ACH |

\*Share of a spot's neighbours with the same niche label; in brackets, the same after shuffling labels within the slide.

### Why this pair is a good first step

1. **There is a published picture to compare against.** These are exactly the slides of Fig. 2b. We can show the authors' non-spatial niches beside SpatialLeiden's domains and discuss the differences visually. The reproduction is in `figures/kuppe_fig2b/`.
1. **They span the range of spatial structure.**
   - P9 has strong structure: an infarct core, a vessel lined by vSMCs, and fibroblast- and myeloid-rich regions. Here SpatialLeiden should produce clear, contiguous domains.
   - P1 has almost none: myogenic niches are intermixed spot by spot, and its coherence is near chance. Here SpatialLeiden should *not* invent domains. P1 is the negative control.
1. **The pipeline already runs on P9.** `03_single_modality_clustering/00_intro_spatialleiden_kuppe.ipynb` does this. As `layer_ratio` goes from 0 to 2, coherence rises from 0.58 to 0.94 at 9 clusters, and agreement with the authors' niches falls (ARI 0.21 → 0.18).
1. **Both slides have matched snRNA-seq and snATAC-seq,** so later modalities can be added without changing the slides.

### What to watch out for

1. **P9 is not a typical slide.** It is the second-deepest run in the cohort, so results there may look better than on an average slide. Its sister section `GT/IZ_P9` (`ACH0012`) has typical depth (3,540 UMIs) and is the natural robustness check.
1. **Do not cluster P1 and P9 together as a comparison.** They differ in patient, zone, time, Visium batch and depth. Treat them as two separate case studies. Joint analysis needs the multi-sample workflow with batch correction, which is a later step.
1. **There is no ground truth.**
   - Agreement with the authors' niches (ARI/NMI) is not accuracy. Those niches ignore space, so a good spatial method is *expected* to diverge from them.
   - Judge results instead by spatial coherence against shuffled labels, by stability across random seeds, and by whether domains follow structures visible on the H&E and in the cell2location maps (vessels, infarct core, scar).
1. **High coherence alone is not success.** A large `layer_ratio` makes any slide look contiguous, P1 included. The real question is whether P9's domains are contiguous **and** biologically interpretable, while P1 is not carved into artificial regions.
1. **The authors' colours are not colourblind-safe.** Keep them only where we compare directly with the paper.

### Verdict

Yes, P1 and P9 make a good first implementation, provided they are run as two single-slide case studies with P1 as the negative control, and P9 is checked against its sister section.

## 3. Suggested next steps, in order

1. **Run the intro notebook's pipeline on P1** with the same settings (same `layer_ratio` grid, number of clusters matched to its 8 niches).
   - Report coherence against shuffled labels, ARI/NMI against both niche sets, and stability over about 5 seeds.
   - Expectation: low coherence gain at small `layer_ratio`. If clear domains appear, check them against the H&E before believing them.
1. **Repeat on GT/IZ_P9 (`ACH0012`)** and compare domain composition with `ACH0013`, using the cell2location proportions per domain.
1. **Add one fibrotic slide,** so each pathologist group has one case. FZ_P14 or FZ_P18 are good candidates. Avoid FZ/GT_P4, the weakest run in the cohort: 954 UMIs and 67 % mitochondrial UMIs.
1. **Only then move to all 27 sections** with the multi-sample workflow, and to the derived modalities (project C step 2).

## 4. Where things are

| What | Where |
|---|---|
| Cohort figure and table | `figures/kuppe_cohort/` (`cohort_table.csv`: one row per sample) |
| Fig. 2b reproduction | `figures/kuppe_fig2b/`, script `01_data_preprocess/kuppe_fig2/fig2b.py` |
| Slide atlas (all 28 runs) | `results/kuppe_2022/slide_atlas/` (local), script `01_data_preprocess/kuppe_slide_atlas/slide_atlas.py` |
| Intro notebook (P9) | `03_single_modality_clustering/00_intro_spatialleiden_kuppe.ipynb`; executed copy in `results/kuppe_2022/intro_notebook/` |
| Data traps | `01_data_preprocess/kuppe_visium_exploration/NOTES.md` |
| Per-spot tables (niches, coordinates, cell2location) | `data/kuppe_2022/visium_tables/<run>_spots.csv` (local, made by `kuppe_visium_exploration/run_all.sh`) |

## 5. Update, evening of 23 November 2026: what was done on P1 and P9

### 5.1 Fig. 2b reproduced

`figures/kuppe_fig2b/`, script `01_data_preprocess/kuppe_fig2/fig2b.py`.

- **The published P9 slide is `GT/IZ_P9_rep2` (`ACH0013`),** identified by its shape. P1 is `10X001`.
- **The shipped `ctniche_1…9` labels are the paper's niches 1–9,** checked against the composition described for Fig. 2c.
- **Colours are the authors' own,** from `saezlab/visium_heart`.
- **Licence:** the paper is under exclusive licence to Springer Nature, so crops of published figures are never committed.

### 5.2 Experiment 0: `layer_ratio` sweep

Folder `03_single_modality_clustering/exp0_layer_ratio_sweep/`: README (design and pre-specified interpretation), `exp0_sweep.py`, `exp0_figures.py`, `RESULT.md`. Panels are in `figures/exp0_layer_ratio/`.

**Design:**

- Slides: P9 rep2, P9 (`ACH0012`) and P1.
- RNA only; `layer_ratio` 0–2.5 (8 values); 5 Leiden seeds.
- The number of domains K is held fixed at every ratio by a resolution search.
- Four K settings: 15 (the brief's rule), 9, and K matched to the authors' cell-type niches (8/8/9) or molecular niches (8/11/12).
- About 2 minutes per K setting, with one process per slide × seed.

**Findings:**

- **The baseline is valid, with a caveat.** `spatialleiden(layer_ratio=0)` is not identical to plain Leiden: the zero-weight spatial layer changes the optimiser's search path. The difference (ARI 0.55–0.84) is no larger than plain Leiden's own seed-to-seed noise. The baseline is therefore `spatialleiden(layer_ratio=0)` itself, with its own resolution search.
- **Coherence rises between 0.25 and 1.0 and saturates above about 1.5 on every slide,** P1 included. Coherence alone cannot separate real structure from structure imposed by the spatial term.
- **P9: space helps at about 0.5.** Domains follow the vSMC-lined vessel, the myocardium and the lesion. Above about 1, the slide is cut into tile-like patches.
- **P1: space over-smooths from about 0.75.**
- **P9 has more spatial structure than P1 even without space:** the gap at ratio 0 is 0.36–0.47 against 0.13.
- **The choice of K does not change the conclusions.** Agreement with the authors' niches stays at most about 0.2 on P9 and about 0.04 on P1, even with K matched to their niche counts. The difference comes from inputs and method, not cluster number.
- **Working range for the next experiments:** `layer_ratio` ≈ 0.5–0.75.

**Code review fixes** (details in `RESULT.md`):

- one coherence definition (edge-weighted) everywhere;
- runs that missed K excluded from every mean and from seed stability;
- a fresh fixed-seed generator for the label shuffles;
- a local copy of the search routine instead of a private API;
- an explicit search for the expression resolution;
- a gene-order assert.

All partitions stayed identical, and numbers moved by at most 0.003–0.027.

### 5.3 Intro notebook completed

`03_single_modality_clustering/00_intro_spatialleiden_kuppe.ipynb` is the single notebook for collaborators.

- **Sections 0–6:** introduction on P9.
- **Sections 7–14:** experiment 0 (K choices, the ratio-0 sanity check, holding K fixed, the full sweep, maps for all slides, conclusions).

It runs end to end in about 1 min 40 s. The committed copy has its outputs cleared; an executed copy is in `results/kuppe_2022/intro_notebook/`.

### 5.4 TF activity (decoupler ULM + CollecTRI) on P1 and P9

`02_modalities/collectri/`, run in env `masld_atlas` (decoupler 2.2.0), used read-only.

- **Matrices:** `data/kuppe_2022/modalities/collectri_ulm_<run>.h5ad` (spots × TFs, in the same barcode order as experiment 0), with 681, 693 and 679 TFs scored. 673 TFs are scored on all three slides.
- **Depth confound:** a typical TF has ≥ 5 targets detected in 16 % of P1 spots, against 70 % on P9 rep2 and 31 % on P9. P1 has fewer UMIs (3,264 against 8,582), 32 % mitochondrial reads, and myocardium-dominated spots.
- **Consequence:** P1's TF layer is exploratory. P9's TF maps follow tissue: MEF2C, GATA4 and TEAD1 in myocardium; SPI1, SMAD3 and STAT3 in lesions.

### 5.5 Environments

- **`spagrn` (new):** Python 3.8, SpaGRN 1.1.4, pySCENIC 0.12.1, hotspotsc 1.1.1, built from `envs/spagrn.yml`. It imports cleanly and reads the cisTarget hg38 v10 databases (folder set by `CISTARGET_DIR`). Not yet run.
- **GPU:** SpatialLeiden 0.4.x is CPU-only (igraph + leidenalg); the two T4 GPUs are unused.
  - An unmerged `origin/gpu` branch (June 2025 draft) reimplements it with RAPIDS cuGraph. It calls `cugraph.leiden` with per-layer resolutions and weights, which released cuGraph versions do not appear to support. Unverified, untested.
  - For single Visium slides, parallel CPU runs are the efficient route.

### 5.6 What is in git

- **Commit `260106e`:** panels (Fig. 2b, experiment 0, TF activity), scripts, the notebook, the experiment-0 folder, CollecTRI and `spagrn.yml`.
- **Not yet pushed:** the code review fixes, the refreshed experiment-0 panels (with `runs_missing_k_per_ratio.csv`), `RESULT.md` and the notebook changes, and this update.
- **Never in git:** the slide atlas PDF (41 MB; rebuilt by `slide_atlas.py`) and all data.

### 5.7 Next steps

1. **Overnight robustness grid, when wanted** (experiment 1; design and pre-specified interpretation to be written first). Two parts:
   - **One factor at a time:** PCs 10–50, variable genes 1,000–5,000, scaling, variable-gene method, `n_iterations`, partition type, how K is reached, graph seeds, and spot filtering on P1.
   - **Full grid for the interacting part:** `layer_ratio` × expression k (10–50) × spatial rings (1–3, uniform or distance weights).

   Estimated at about 5,500 runs, 1–3 h on 28 cores. The fuller version (all 8 ratios × 4 K) is about 20,000 runs, 4–10 h. The question is whether the conclusions survive these choices, not which setting wins.
2. **First multimodal experiment on P9:** RNA + TF activity (+ space) with `spatialleiden_multimodal`, `layer_ratio` ≈ 0.5–0.75, over a grid of TF-layer weights. P1 only as exploratory.
3. **SpaGRN on P9 rep2,** and a comparison of its regulons with the CollecTRI TF activities. SpaGRN infers regulons per slide, so TF sets will differ between slides.
4. **Check the cohort figure's filled markers** against the snRNA-seq and snATAC-seq objects (now downloaded).
