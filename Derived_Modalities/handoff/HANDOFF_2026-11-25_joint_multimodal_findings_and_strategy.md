# Handoff, 24–25 November 2026: joint multimodal SpatialLeiden on 27 slides, findings and strategy

Track C, branch `hack/c-derived-modalities`. Paths are relative to `Derived_Modalities/`. This continues `HANDOFF_2026-11-24_joint_three_slides.md`, which covers the setup, environment, commands and the per-slide vs joint design. **Nothing from 24 Nov evening onwards is committed yet** (section 6).

## 1. What was run

All runs cluster **all 27 Visium slides jointly** (88,704 spots) after Harmony integration (patient + sample + Visium batch). Common settings:

- **Graphs:** each layer gets a 15-NN graph on its Harmony embedding; the spatial graph is the Visium grid, built per slide.
- **K:** fixed at 9 by a resolution search.
- **Seed:** 0 only.

| Layer(s) | Spatial weights | Output (`results/kuppe_2022/`) |
|---|---|---|
| RNA | 0, 0.25, 0.75, 1.2 | `exp2_joint_harmony` |
| TF (CollecTRI, 27 slides) | 0, 0.25, 0.75, 1.2 | `exp2_joint_harmony_tf` |
| CCC (our LIANA+ reconstruction, 286 pairs) | 0, 0.25, 0.75, 1.2 | `exp2_joint_harmony_ccc_ours` |
| RNA 1 + TF 0.5 | 0, 0.25, 0.75, 1.2 | `exp3_joint_rna_tf` |
| RNA 1 + CCC 0.5 | 0, 0.25, 0.75, 1.2 | `exp3_joint_rna_ccc_ours` |
| RNA 1 + TF 0.5 + CCC 0.5 | 0, 0.75 | `exp3b_joint_rna_tf_ccc` (**running** at handoff, screen `joint_rna_tf_ccc`; new script `exp3b_joint_rna_tf_ccc.py`) |

**Search outcome:** every run reached 9 domains, except RNA at spatial weight 0.25 (8). No resolution search failed.

**Our CCC layer** (`02_modalities/liana/liana_local_all27.py`) is a reconstruction of the colleague's scores, which exist for 4 slides only:

- calibrated bandwidth: 0.75 × spot spacing;
- median per-pair r with the colleague's scores: 0.79–0.84 on the 4 shared slides.

## 2. Figures

All are in `results/kuppe_2022/`, with copies in `quick_look/` next to the repo.

| Figure | Script | Shows |
|---|---|---|
| `joint_three_slides/fig_joint_rna_tf.png`, `fig_joint_ccc.png` | `03_.../exp2_joint_harmony/fig_joint_three_slides.py` | The 23 Nov P1 / P9 / P9 rep2 multipanels, with joint domains |
| `joint_domain_enrichment/fig_joint_rna_vs_rna{tf,ccc}.png` | `fig_joint_rna_vs_rnatf.py [tf\|ccc]` | RNA vs single layer vs combination: maps on 4 slides, where domains come from, mean share per region |
| `joint_domain_enrichment/<run>_fig.png` | `domain_enrichment_by_zone.py` | Control vs patient enrichment per domain (slide = unit) |
| `joint_umaps/fig_joint_umaps.png`, `_v2.png` | `07_integration/fig_joint_umaps.py [v2]` | **Main figure**, see below |

**Main figure layout:**

- **Style:** Kuppe ED Fig. 3a–b, Nature format (183 mm wide, Arial 5 / 6 / 8 pt).
- **Panels a–c:** UMAPs of the RNA, TF and CCC Harmony embeddings.
- **Columns:**
  - region, patient (the authors' palettes);
  - the authors' cell-type and molecular niches;
  - our domains.
- **Panel d:** stacked bars of each domain's share per region.
- **The two versions:** v1 ends with TF-activity and communication columns; v2 replaces them with the RNA-only / TF-only / CCC-only clusterings (no spatial graph).
- **Stressed-CM domain:** drawn in pink with a bold label, plus a dashed ring on the RNA map.
- **New column:** RNA + TF + CCC + space is added automatically once `exp3b` has finished.

## 3. Findings

1. **Integration works.** Patients mix on all three UMAPs, and the domains barely follow slides (ARI against sample ≈ 0.03–0.05). Agreement with the authors' cell-type niches is the same or slightly better than per slide for RNA and RNA + TF (P9 rep2: 0.25).
2. **Control vs patient.** Healthy myocardium (CM + EC) is enriched in controls; inflamed scar (fibroblasts + myeloid, SPP1 / POSTN high) and fibrotic scar are enriched in the ischaemic and fibrotic zones. This holds for every clustering (BH-adjusted p < 0.05 with 4 control slides).
3. **Stressed-cardiomyocyte domain** (RNA + TF D5, RNA + CCC D6):
   - **Markers:** NPPB and ANKRD1 high; cardiomyocytes plus proliferating cells.
   - **Where:** border and ischaemic zones, ≥ 2 % of the slide on 16 of 27 slides from 6 patients, with NPPB higher inside it on 15 of 16 and ANKRD1 on 16 of 16. RNA + TF and RNA + CCC mark the same spots (Jaccard 0.69 median).
   - **On the UMAP** it sits between healthy myocardium and scar.
   - **Not significantly enriched by zone:** slide shares vary a lot, and there are only 3 border-zone slides.
4. **It is not invisible to RNA.** RNA + space at weight 0.25 finds a similar domain (Jaccard 0.67). At 0.75 and above, RNA + space merges it into a larger fibroblast + CM domain. TF or CCC keep it separate at the spatial weight that integrated data needs (preprint: about 1.2). That is the honest claim: **the extra layers stabilise a disease-relevant transition state against spatial smoothing; they do not discover it from nothing.**
5. **No other truly multimodal-only domain** (`joint_domain_enrichment/multimodal_domains_vs_rna.csv`):
   - every multimodal domain has an RNA-only counterpart at some weight (Jaccard ≥ 0.5);
   - the one exception is RNA + CCC D4, a second CM + EC domain in remote, ischaemic and control tissue (best RNA match 0.49). It is borderline.
6. **TF alone and CCC alone are noisier and partly circular:**
   - TF-only domains separate spots by TF activity (η² 0.59); CCC-only domains separate them by communication (η² 0.74);
   - with RNA at weight 1 and the extra layer at 0.5, RNA dominates (`joint_umaps/eta2_activity.csv`);
   - TF-only and CCC-only partitions agree poorly with RNA (ARI 0.12–0.33; `joint_umaps/ari_between_clusterings.csv`).
7. **Authors' niches:**
   - their molecular niches follow the RNA embedding;
   - their cell-type niches follow none of our embeddings, because they cluster cell2location proportions, not expression.

## 4. Strategy: finding further disease-specific domains systematically

The aim is to find spatially coherent, reproducible domains that are enriched in a disease zone and missed by RNA alone or by the authors' niches. Releasing K is part of it, but on its own it mostly produces more splits of the large domains. A systematic search needs three axes and a filter.

### 4.1 Explore the parameter space (overnight, resumable)

| Axis | Values | Why |
|---|---|---|
| Resolution (K free) | expression resolution set to give K ≈ 9 / 15 / 25 / 40 at spatial weight 0, then kept fixed (as in `exp4_unfixed_k.py`) | Small, specific domains appear only at higher K. Fixing the resolution instead of K is the preprint's discovery mode |
| Spatial weight | 0.25, 0.75, 1.2, 2 | The preprint uses 1.2–2.5 for integrated data. Fragile states (such as the stressed CMs) are lost at high weight unless another layer holds them |
| Layer combination and weight | RNA; RNA + TF; RNA + CCC; RNA + TF + CCC; extra-layer weight 0.25 / 0.5 / 1 | Shows which layer protects which domain |
| Seeds | 0–4 | Needed for stability (step 4.2) |

**Cost:** about 4 K × 4 weights × 4 combinations × 2 extra-layer weights × 5 seeds ≈ 640 runs. One joint run takes 1–3 min with fixed resolutions (no search), so about 10–20 h on 30 cores. **Start with seed 0 to prune the grid**, then add the seeds.

### 4.2 Filter: what counts as a disease-specific domain

Apply to every domain of every run (a "domain table"):

1. **Reproducible:** it reappears with Jaccard ≥ 0.5 in ≥ 4 of 5 seeds, and in neighbouring settings (the next K or the next weight).
2. **Spatially coherent:** neighbouring spots share it more than chance (coherence gap > 0.2 for that domain); it is not salt-and-pepper.
3. **Not one slide:** ≥ 2 % of the slide on ≥ 3 slides from ≥ 3 patients.
4. **Disease-enriched:** a higher share in a patient region than in control (slide as unit; with only 4 controls, use effect size, at least 3× the control share, plus a permutation test). Also check patient group (ischaemic vs fibrotic) and time since infarction where the metadata allow.
5. **Hidden from the baselines:** best Jaccard < 0.5 with
   - any RNA-only domain at the same K (any spatial weight);
   - the authors' cell-type and molecular niches.

   A domain that passes 1–4 but fails 5 is still useful (disease biology), but it is not a multimodal gain.

### 4.3 Characterise the survivors

- **Composition:** cell2location proportions and the authors' sn-derived cell states. For example, is it vCM3 (stressed cardiomyocytes) in snRNA-seq?
- **Markers:** marker genes (Wilcoxon within slide, then pooled), top TFs (CollecTRI) and ligand-receptor pairs specific to the domain. Check whether the layer that keeps it separate also explains it biologically.
- **Location:** maps on every slide where it occurs, and where it sits relative to the infarct (border, perivascular, subendocardial).
- **Literature:** match against Kuppe's published states (vCM3, SPP1+ macrophages, myofibroblasts, the border-zone niche).

### 4.4 Statistics

Slides, not spots, are the unit. With 4 controls, 5 RZ, 3 BZ, 9 IZ and 6 FZ slides:

- report shares per slide and effect sizes;
- use BH correction over all candidate domains;
- make clear that discovery and testing happen on the same 27 slides, so the findings are exploratory. **Validation option:** hold out patients (leave-one-patient-out), or check in another MI Visium dataset.

### 4.5 Two caveats to keep in mind

- **More K means more domains that are just finer cuts of a gradient** (for example the border-zone transition). Seed and neighbouring-setting stability (4.2.1) is the guard.
- **Our CCC is a reconstruction.** Any CCC-driven candidate should be rechecked with the colleague's scores, or with their method once documented.

## 5. Immediate next steps

1. **The running `exp3b` run.** When `EXP3B_EXIT` = 0, rerun both versions of the main figure (`07_integration/fig_joint_umaps.py ... [v2]`); the RNA + TF + CCC column and its stressed-CM domain (if Jaccard > 0.5) are added automatically. Also run `domain_enrichment_by_zone.py ... exp3b_joint_rna_tf_ccc r0.75_s0`.
2. **Five-seed versions of the quick runs.** Drop `EXP_SEEDS` in the commands of the previous handoff, section 5.
3. **Strategy step 4.1 with seed 0.** Extend `exp4_unfixed_k.py` / `exp5_unfixed_k_multimodal.py` with a list of target K values; they already run one clustering per weight without a search.
4. **Write the domain-table script for 4.2.**
5. **Commit** (only when asked): the new scripts, the handoffs, and small results with `git add -f` (label CSVs, tables, key PNGs). The new scripts are:
   - `02_modalities/liana/liana_local_all27.py`
   - `03_.../exp2_joint_harmony/`:
     - `fig_joint_three_slides.py`
     - `domain_enrichment_by_zone.py`
     - `fig_joint_rna_vs_rnatf.py`
     - `exp3b_joint_rna_tf_ccc.py`
   - `07_integration/fig_joint_umaps.py`

## 6. Not committed, not backed up

Everything from 24 Nov evening onwards is on this VM only: the joint embeddings, labels, figures and our 27-slide CCC tables. Before switching VMs, copy them to the OneDrive folder `SpatialLeiden_hackathon/` and commit the scripts.

## 7. Grid search results (run 25 Nov, 01:21–09:14 CET)

**Code:** `08_grid_search/`:

- `grid_search.py` (calibrate / grid; resumable)
- `run_grid.sh` (chain)
- `domain_table.py` (filter of section 4.2)
- `candidate_families.py` (collapses candidates into regions)
- `fig_region_tf_ccc.py` (maps + top TFs / LR pairs per region)

**Outputs:** `results/kuppe_2022/grid_search/`. The runs folder (40 MB) is not in git; it is on OneDrive.

**What ran:**

- 10 layer configurations × target K 9 / 15 / 25 / 40 × spatial weights 0 / 0.25 / 0.75 / 1.2 / 2 × seeds 0–4: 1,000 runs, all done.
- Calibration hit the target K in 38 of 40 configurations; RNA + TF + CCC at weight 0.25 and at weight 1 reached 24 instead of 25.
- The expression resolution is fixed per configuration and K. The spatial resolution is fixed at 1.0, so K drifts with spatial weight (for example RNA, K 9 at weight 2 → 31 domains).

**Filter result** (`domain_table.csv`, `candidates.csv`):

- 32,660 domains in total;
- 3,268 disease-enriched (reproducible in at least 3 of 4 other seeds, coherent, on ≥ 3 slides from ≥ 3 patients, ≥ 3× the control share);
- 504 of these also pass "Jaccard < 0.5 with RNA alone and with the authors' niches";
- these 504 instances collapse into **47 regions** (`candidate_families.csv`, consensus spots in `family_consensus_spots.csv.gz`).

**Top regions** (TFs and LR pairs from `fig_region_tf_ccc.py`; within-slide effect, ≥ 75 % of slides positive):

| Family | Found in | Cell types / where | Top TFs | Top LR pairs |
|---|---|---|---|---|
| 3 | 207 runs, 9 configurations | Fib + CM + EC, IZ 94 % | MEF2A, SRF, MYOG, NKX2-5, HAND2 | laminins → ITGA7/ITGB1 |
| 0 | 45 runs, 6 configurations | Fib + Myeloid + Prolif, IZ 94 % | SPI1, MYB, FOS, AP-1, STAT1, REL | APOE → LRP1, CD14 → ITGB2, CCN2 → ITGB2, C1QB → LRP1, LGALS9 → CD44 |
| 13 | 33 runs, 8 configurations | CM + EC + Myeloid, IZ 96 % (3 slides) | MEF2C, NR1H4, PPARG | C3 → CD81, PLA2G2A → ITGB1 |
| 7 | 11 runs, 3 configurations | CM + EC + Fib, **BZ 84 %** | MEF2C/A, PPARA, PPARGC1A | laminin → ITGA7, BSG → SLC16A1, MFGE8 → ITGAV |
| 2 | 11 runs, 2 configurations | Fib + CM + EC, **FZ 66 %** (15 slides) | NKX2-5, GATA4, HOPX | THBS4 → CD36, collagens → CD36 |
| 6 | 8 runs, 4 configurations | Fib + Myeloid + EC, **FZ 100 %** (3 slides) | STAT6, SOX6 | complement C3 / C1QB → LRP1 / CD81 |

**How to read them:**

- Many regions sit just under the 0.5 Jaccard threshold with RNA alone (max 0.46–0.50): they are refinements of RNA domains, not entirely new regions.
- The fibrotic-zone regions (families 2, 5, 6, 9, 15, 22) overlap least with RNA alone (0.21–0.37).
- "Fold vs control" is inflated for regions nearly absent from controls; the permutation p-values are weak with 4 control slides.

**Stressed-CM domain across the grid** (`stressed_cm_in_grid.csv`):

- recovered (Jaccard ≥ 0.5) in 25–40 % of RNA + TF runs at spatial weights 0–1.2, against 15–30 % for RNA alone;
- best for RNA + TF at K = 15 (median 0.52–0.55, against 0.42 for RNA);
- never at spatial weight 2.

**This corrects the impression from the fixed-K runs:** the gain over RNA is real but modest and depends on settings.

**Statistics, next:**

- **scCODA** (pertpy): compositional test of region shares per slide across zones, slide as unit.
- **Milo:** neighbourhood differential abundance on the RNA + TF / RNA + CCC graph, independent of clustering, with patient as a covariate (several patients contribute several slides).

Both are limited by 4 control and 3 border-zone slides.
