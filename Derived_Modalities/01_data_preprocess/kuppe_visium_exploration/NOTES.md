# Kuppe et al. 2022: notes from the Visium exploration

Things that surprised us, or that need handling before the data can be used. Written 23 November 2026, during the descriptive pass (no new analysis).

## Where the usable Visium data is

1. **The Zenodo "Processed Data 1" Visium files (`Visium_*.h5ad`) are not enough on their own.**
   - They are in CELLxGENE schema 2.0. `X` is log-normalized and the raw counts are in `.raw`. They hold 13,490 genes.
   - They do not contain the published cell-type or molecular niches.
   - Spot coordinates (`obsm["X_spatial"]`) are rescaled to 0–1, separately along x and y. This can distort the hexagonal grid.
   - They contain no image.
1. **The full annotation is in the authors' Seurat objects** inside each Space Ranger archive (Zenodo 6580069, `10X_Visium_<id>.tar.gz` → `<id>/<id>.rds`).
   - They hold `celltype_niche`, `molecular_niche`, cell2location abundances (`c2l`) and proportions (`c2l_props`), plus the original spot coordinates.
   - They also contain `SCT` and **`progeny`** assays. PROGENy is shipped but was not used in this pass (out of scope).
   - We extract per-spot tables from them with `01_extract_seurat_tables.R`.
1. **The spot diameter stored in the Seurat objects is wrong.** `scale.factors$spot` equals the hires scale factor (0.203). Use `scalefactors_json.json` from the Space Ranger `spatial.zip` (e.g. 79.5 px for ACH006).
1. **There are 28 Visium runs, and 91,517 spots in total.** This matches the paper's Fig. 2a.
   - The 28 runs cover 27 tissue samples.
   - RZ_P6 has two runs, ACH0021 and ACH0022. Only ACH0022 is in the Zenodo h5ad files and in the paper's sample-level analyses (n = 27).
1. **RZ_P6 also has two H&E images** (`Visium_2_CK280.tif`, `Visium_17_CK295.tif`). No column says which image belongs to which run. Every other image maps one-to-one.

## Labels and names

6. **Zones are per section, not per spot.** No spot-level zone or histology annotation ships with the data.

1. **`major_labl` (sic) collapses mixed regions.**

   - RZ/BZ is counted as BZ, IZ/BZ as IZ, and GT/IZ as IZ.
   - RZ/FZ_P5 is labelled FZ, but its `patient_group` is group_1 (myogenic).
   - The full region name is in `region_novel` / `patient_region_id`. "GT" appears in region names, but the paper's zone list does not define it.

1. **Patient groups:** group_1 = myogenic (control, RZ, BZ), group_2 = ischaemic, group_3 = fibrotic.

1. **Sample identifiers differ between files.**

   | Where | Form | Example |
   |---|---|---|
   | `patient_region_id` | uses `/` | `GT/IZ_P9` |
   | Zenodo h5ad names | use `_` | `Visium_GT_IZ_P9.h5ad` |
   | Space Ranger archives | `hca_sample_id` | `ACH0012` |
   | H&E files (most) | CK numbers | `Visium_9_CK287.tif` |
   | H&E file for control_P1 | AKK number | `157771_2019AKK006.tif` |

   `metadata-Visium.csv` and the H&E `metadata.xlsx` link them through `patient_region_id`.

1. **Cell-type names differ between objects.**

   - The h5ad `obs` uses Cardiomyocyte, Fibroblast, Pericyte, Cycling.cells, …
   - The Seurat `c2l` assays use CM, Fib, PC, prolif, …
   - There are 11 types, not the 8 major types the handoff mentions.

1. **Niche names:** `celltype_niche` values are `ctniche_1`…`ctniche_9` and `molecular_niche` values are `molniche_*` (12 in total). The paper calls them "niche 1–9" and "molecular niche 0–…". Their correspondence to the paper's numbering still needs checking against Fig. 2c before we rely on it.

## Data properties worth knowing

12. **Mitochondrial fraction is very high in myogenic sections.** The median % mitochondrial UMIs per spot is 31–40 % in control, RZ and BZ, but 11 % in IZ. A standard 20 % filter would remove most of the myocardium.
01. **Fibrotic sections are the shallowest, not the ischaemic ones.**
    - Visium medians per zone: FZ 1,797 UMIs and 1.3 estimated cells per spot. IZ: 4,774 UMIs. BZ: 5,032 UMIs.
    - The single-nucleus data is the other way round: according to the authors' sample tables, IZ has the fewest nuclei (median 3,654 snRNA and 434 snATAC nuclei per sample).
01. **The estimated cells per spot (sum of cell2location abundances) rise with depth.** Spearman ρ = 0.84 both per spot and per section. It is therefore not an independent cellularity measure.
    - The paper counted nuclei per spot from the H&E images (VistoSeg), but those counts do not ship with the data.
01. **Sections can contain few niches:** IZ_P3 (10X0017) has only 4 of the 9 cell-type niches.

## What was left out, and why

16. **Out of scope for the Visium-only pass:** snRNA-seq / snATAC-seq QC (handoff 3.3), existing embeddings (3.7) and cross-modality composition concordance (3.8).
    - The snRNA-seq and snATAC-seq objects were still downloading when this note was written.
    - The availability figure takes nuclei counts from the authors' metadata sheets (Supplementary Tables 2 and 3), not from the objects. They should be checked against the objects once these are downloaded.
01. **Time after symptom onset is not in any shipped table**, so the cohort figure has no time axis. It would have to come from the paper's supplementary tables.
