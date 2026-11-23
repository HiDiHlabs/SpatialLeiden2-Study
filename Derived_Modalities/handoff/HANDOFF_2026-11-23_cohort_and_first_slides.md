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
