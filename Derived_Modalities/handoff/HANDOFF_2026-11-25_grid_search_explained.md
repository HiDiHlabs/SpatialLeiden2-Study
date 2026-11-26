# How the overnight grid search worked

Track C, 24–25 November 2026. This is a plain-language explanation for presenting the logic. Technical details, numbers and file locations are in `HANDOFF_2026-11-25_joint_multimodal_findings_and_strategy.md` (section 7). Code: `Derived_Modalities/08_grid_search/`.

## 1. The question

**Does adding a layer derived from gene expression help SpatialLeiden find disease-specific tissue regions that gene expression alone, or the authors' niches, miss?** The two derived layers are transcription-factor (TF) activity and cell-cell communication (CCC).

## 2. The data and the starting point

- **Data:** 27 Visium slides from Kuppe et al. 2022 (88,704 spots) covering control, remote, border, ischaemic and fibrotic zones after myocardial infarction.
- **Three layers per spot:**
  - **RNA:** gene expression.
  - **TF:** transcription-factor activity, computed with CollecTRI / decoupler.
  - **CCC:** ligand-receptor communication scores, computed with LIANA+.
- **Integration:** each layer was integrated across slides with **Harmony** (correcting for patient, sample and Visium batch). All slides are clustered together, so a domain means the same thing on every slide.
- **Graphs:**
  - each layer gives one graph: each spot linked to its 15 most similar spots;
  - the tissue gives a spatial graph: each spot linked to its physical neighbours on its own slide.

## 3. How SpatialLeiden combines layers: weights

**SpatialLeiden looks for one grouping of spots that fits all graphs at once.** Think of it as a vote:

- **RNA says:** "these spots have similar gene expression, group them."
- **TF says:** "these spots have similar TF activity, group them."
- **Space says:** "these spots are next to each other, group them."

When the layers disagree, **the weights decide whose vote counts more**:

| Setting | Meaning |
|---|---|
| RNA 1 | the reference; always 1, because only the ratios between weights matter |
| TF 0.25 | TF's vote counts a quarter as much as RNA's: a light nudge |
| TF 0.5 | half as much as RNA |
| TF 1 | equal say with RNA |
| spatial 0 | location ignored |
| spatial 2 | strong push towards continuous regions |

**Example:** two spots have similar RNA but different TF activity. With TF at 0.25, RNA wins and they stay together. With TF at 1, TF can pull them apart.

## 4. Why a grid, and what it contains

**There is no rule for choosing these settings.** The SpatialLeiden 2 preprint picks them per dataset. Instead of defending one arbitrary choice, we ran every combination of four "dials":

| Dial | What it controls | Options | Count |
|---|---|---|---|
| 1. Layers and their weights | which layers are combined, and how loud the extra one is | RNA alone; RNA + TF; RNA + CCC; RNA + TF + CCC, with the extra layer(s) at 0.25 / 0.5 / 1 | 10 |
| 2. Granularity | how fine the domains are: resolution set to give about this many domains | 9, 15, 25, 40 | 4 |
| 3. Spatial weight | how much location counts | 0, 0.25, 0.75, 1.2, 2 | 5 |
| 4. Random seed | same settings, different random start, to separate stable results from flukes | 0–4 | 5 |

**Total:** 10 × 4 × 5 × 5 = **1,000 clustering runs**, about 7 hours on 28 cores.

**The 10 layer combinations:**

| # | Layers | Weights (RNA : extra) |
|---|---|---|
| 1 | RNA only | 1 |
| 2–4 | RNA + TF | 1 : 0.25, 1 : 0.5, 1 : 1 |
| 5–7 | RNA + CCC | 1 : 0.25, 1 : 0.5, 1 : 1 |
| 8–10 | RNA + TF + CCC | 1 : 0.25 : 0.25, 1 : 0.5 : 0.5, 1 : 1 : 1 |

**How the granularity was fixed.** For each layer combination and target number of domains, the resolution was first tuned once (seed 0, no spatial weight) to give that number of domains. It was then kept fixed for all spatial weights and seeds. The number of domains is therefore free to change as space is added, which shows whether location merges or splits regions. This is the preprint's "discovery mode".

## 5. From 32,660 domains to 47 regions

**The 1,000 runs produce 32,660 domains.** Each domain was kept only if it passed **all five** criteria:

1. **Reproducible:** it reappears in at least 3 of the other 4 seeds.
2. **Spatially coherent:** neighbouring spots share it well above chance, so it is a region, not scattered spots.
3. **Not one sample:** it covers ≥ 2 % of the slide on at least 3 slides from at least 3 patients.
4. **Disease-enriched:** its share in some patient zone is at least 3× its share in control slides, with slides as the unit.
5. **Not already known:** less than 50 % overlap (Jaccard) with every RNA-only domain at the same granularity, and with the authors' cell-type and molecular niches.

| Step | Domains |
|---|---|
| All domains | 32,660 |
| Pass criteria 1–4 (robust and disease-enriched) | 3,268 |
| Pass all five (also new compared with RNA and with the authors) | 504 |
| Grouped into distinct regions (instances overlapping ≥ 50 %) | **47** |

The 504 collapse to 47 because the same region recurs across many settings. **How often a region recurs, and across how many layer combinations, is a measure of how robust it is.** For example, family 3 appears in 207 runs across 9 combinations.

## 6. Characterising each region

**Method.** For every region, on every slide where it covers at least 2 % of the spots:

- compare TF activity and ligand-receptor scores **inside vs outside** the region, **within that slide**, in units of that slide's standard deviation, so patient and slide differences cannot drive the result;
- average over slides;
- count a feature as characteristic if it is higher inside on **at least 75 % of those slides**.

**Main regions found:**

| Theme | Example families | TF programme | Signalling |
|---|---|---|---|
| Border-zone myocardium | 7 | MEF2, PPARGC1A | laminin → ITGA7/B1 |
| Stressed cardiomyocytes | 33, 41 | heat shock (HSF2/4), ER stress (XBP1, ATF6) | CCN1 / PSAP / GRN → SORT1, integrins |
| Ischaemic inflammatory macrophages | 0 | SPI1, AP-1, STAT1 | APOE / C1QB → LRP1, CD14 → ITGB2 |
| Fibrotic complement stroma | 6 | STAT6, SOX6 | C3 → LRP1 / CD81 |
| Fibrotic matrix | 2 | GATA4, NKX2-5, HOPX | THBS4 / collagens → CD36 |

**The fibrotic regions overlap least with RNA alone**, so they are the clearest cases where the extra layer adds structure.

## 7. Honest limitations

- **The stressed-CM domain found the day before** is recovered in only 25–40 % of RNA + TF runs, against 15–30 % for RNA alone. The gain is real but modest and depends on the settings. Exposing this is exactly what the grid is for.
- **Many regions are refinements of RNA domains** (overlap just under 0.5), not entirely new tissue types.
- **Small sample sizes:** 4 control and 3 border-zone slides. The "fold vs control" can be inflated for regions almost absent from controls.
- **Partial circularity:** for regions found with TF as a layer, TF enrichment is partly built in. The ligand-receptor results are not affected when the region came from RNA + TF.

**Statistical confirmation (running at handoff time):**

- **scCODA:** a compositional model of how domain shares per slide differ between zones.
- **Milo:** a neighbourhood test that does not depend on any clustering.

## 8. One-line summary for a slide

> Instead of choosing one set of parameters, we clustered all 27 MI Visium slides 1,000 times, across every combination of layers, granularity, spatial weight and random seed. We kept only regions that are reproducible, spatially coherent, present in several patients, enriched in a disease zone, and not seen with RNA alone or in the authors' niches. The result is 47 regions with zone-specific TF and communication programmes.

## 9. Figures

- **Summary multipanel:** `results/kuppe_2022/grid_search/fig_families_multipanel.png/.pdf`
  - a, b: characteristic TFs and ligand-receptor pairs of 7 representative regions;
  - c: share of slide per zone;
  - d: maps on representative slides;
  - e: how often the stressed-CM domain is recovered across the grid.
- **One figure per region:** `results/kuppe_2022/grid_search/regions/family_<n>/` (all 47).
- **Table of all regions:** `results/kuppe_2022/grid_search/family_tf_lr_summary.csv`.
