# COSIE mouse brain (6-OHDA), striatum: SpatialLeiden with MSI preprocessing `log`

Step-by-step description of what [`04_all_fmp10_slides.ipynb`](04_all_fmp10_slides.ipynb) does with

```python
DATASET = "striatum"
msi_preprocessing = "log"
```

and of the parameter scan in [`05_param_scan.ipynb`](05_param_scan.ipynb). Both notebooks are committed with their outputs. Results: `results/Cosie_Clustering/striatum/msi_log/` (not in git).

All numbers below come from that run, with a single Leiden seed (`seed = 42`) unless stated otherwise.

---

## 1. Data

Vicari et al., *Nat Biotechnol* 2024 (Mendeley Data `w7nw4km7xd`). We use the FMP-10 mouse slide `V11L12-109` (striatum).

| section | mouse | MALDI | RNA spots | MSI spots |
|---|---|---|---|---|
| A1 | mPD1 | yes | 2,675 | 2,206 |
| B1 | mPD3 | yes | 3,120 | 2,526 |
| C1 | mPD4 | yes | 2,918 | 2,379 |
| D1 | mPD3 | no (control, `Laser = nL`) | 3,017 | – |

- **RNA:** Visium counts from Space Ranger (`filtered_feature_bc_matrix.h5`), in-tissue spots only, 32,285 genes.
- **MSI:** MALDI-MSI already registered to the Visium spots by the original authors (`misc/se.multi.list`). It was exported to Matrix Market with `00_export_msi_spots.R` and has 1,538 m/z features. The values are continuous intensities, not counts. About 18 % of the in-tissue spots of A1–C1 have no MSI value.
- **Annotation (per spot, from the original study):**
  - `RegionLoupe`: 9 regions, namely ACB, CP, CTX, LSX, MSC, OT, VL, aco and cc.
  - `lesion`: intact or lesioned hemisphere.
  - The lesion is also scored within the lesion region only (CP, the dopaminergic terminals; `in_lesion_region`).
- **Cell-type proportions:** the stereoscope output of the original study (`misc/Zeisel_stsc_output.csv`). It is used for plotting only.

## 2. RNA preprocessing

1. Concatenate all four sections.
2. `normalize_total` (library size, target = median), then `log1p`.
3. 3,000 highly variable genes (`batch_key = "section"`).
4. PCA with 50 components on the HVGs.
5. **Harmony** on the PCs (`harmonypy`, batch = section, run within the slide) → `X_pca_harmony`. The RNA sections separate in the UMAP without it, partly because D1 had no MALDI processing.
6. kNN graph (k = 15) on `X_pca_harmony`.

## 3. MSI preprocessing (`log`)

```python
msi = ad.concat(msi_sections, merge="same")      # A1, B1, C1
msi.layers["counts"] = msi.X.copy()
sc.pp.normalize_total(msi)                       # TIC normalisation, target = median TIC
sc.pp.log1p(msi)
sc.pp.pca(msi, n_comps=30)                       # all 1,538 m/z, no scaling
sc.pp.neighbors(msi, n_neighbors=15, use_rep="X_pca")
```

- **No feature selection.** The lesion markers (dopamine m/z 674.28833, its isotopes, m/z 390.17, 435.21 and others) are detected in only 15–40 % of spots, essentially the intact CP only. A detection filter (for example ≥ 50 % of spots) removes exactly these features. With that filter, the lesion ARI in CP dropped from about 0.6 to 0.02.
- **No Harmony for MSI** (`msi_rep = "X_pca"`). The MSI sections separate only moderately. The table below gives the share of each spot's 15 MSI neighbours that come from its own section. Harmony is available as an option (`msi_rep = "X_pca_harmony"`).

| same-section neighbours (MSI) | before Harmony | after Harmony | expected if mixed |
|---|---|---|---|
| A1 | 0.53 | 0.40 | 0.31 |
| B1 | 0.42 | 0.37 | 0.36 |
| C1 | 0.47 | 0.43 | 0.34 |

- **Alternative `msi_preprocessing = "maldipy"`.** This follows MALDIpy `maldi_norm`: TIC, regress out TIC, then scale each m/z, per slide. Its results go to `msi_maldipy/`. It separates the regions somewhat better, and `log` separates the lesion somewhat better.

## 4. Graphs for SpatialLeiden

- **Spatial graph:** the Visium hexagonal grid. Each spot is linked to its 6 direct neighbours (`squidpy.gr.spatial_neighbors`, `coord_type="grid"`). The graph is built per section, so it has no edges between sections.
- **Padded modality graphs:** each modality graph is padded to all spots in `mdata.obsp["<mod>_connectivities"]`. Spots without that modality have no edges in its layer.
- **Smoothed RNA (control without MSI):** the RNA PCs averaged over each spot and 2 rings of spatial neighbours, then standardised and turned into a kNN graph (k = 15).

The RNA and MSI kNN graphs are built jointly over all sections. One `multiplex_leiden` call therefore gives **one shared label set**, so the same cluster means the same domain in every section.

## 5. Sanity checks and marker figure

- **UMAPs:** `rna_umap_harmony.png` and `msi_umap_harmony.png` show the embeddings before and after Harmony, coloured by section, region and lesion, plus one panel per section.
- **Dopamine figure** (`dopamine_msi_pcp4_msn1.pdf/.png`): dopamine (m/z 674.28833, log TIC-normalised), *Pcp4* (RNA) and MSN1 (stereoscope proportion of striatal medium spiny neurons). The figure follows Vicari et al. 2024. Dopamine is high only in the CP of the intact hemisphere.

## 6. SpatialLeiden runs

`multiplex_leiden` is run on A1–C1, using all spots of these sections; spots without MSI get no MSI edges. Each layer has a resolution `r` and a weight `w`. Settings: `n_iterations = 5`, `random_state = 42`.

Scores are computed on the annotated spots **with MSI**:
- ARI against `RegionLoupe`
- ARI against lesion side
- ARI against lesion side within CP

### 6a. Default parameters (all layers r = 1, spatial w = 0.5)

| run | clusters | ARI region | ARI lesion | ARI lesion in CP |
|---|---|---|---|---|
| RNA | 10 | 0.53 | 0.00 | 0.00 |
| MSI | 11 | 0.12 | 0.09 | 0.58 |
| RNA + MSI | 9 | 0.46 | 0.11 | 0.80 |
| RNA + smoothed RNA | 13 | 0.32 | 0.09 | 0.81 |
| RNA + MSI + smoothed RNA | 8 | 0.33 | 0.13 | 0.76 |

### 6b. Same number of clusters (21) for every run

For each run, every combination of layer weights is tested:
- modality layers: 0.5, 1 or 2, relative to the first layer
- spatial layer: 0.25 or 0.5

For each combination, one resolution shared by all layers is bisected until there are exactly 21 clusters. The combination with the highest ARI against "regions, CP split by lesion side" is kept. Because the parameters are chosen with the annotation, these scores are an upper bound. All combinations are saved in `spatialleiden_21_clusters_region+local lesion.csv`.

| run | ARI region | ARI lesion | ARI lesion in CP |
|---|---|---|---|
| RNA | 0.26 | 0.02 | 0.19 |
| MSI | 0.08 | 0.04 | 0.30 |
| RNA + MSI | 0.27 | 0.05 | 0.54 |
| RNA + smoothed RNA | 0.26 | 0.03 | 0.25 |
| RNA + MSI + smoothed RNA | 0.28 | 0.06 | 0.60 |

### 6c. Manual parameters (`MANUAL_RUNS`), used for the figures

| run | parameters | clusters | ARI region | ARI lesion in CP |
|---|---|---|---|---|
| RNA | rna r=1.2 w=1, spatial r=0.8 w=0.8 | 9 | 0.46 | 0.03 |
| MSI | msi r=1.2 w=1, spatial r=0.8 w=0.6 | 8 | 0.14 | 0.61 |
| RNA + MSI | rna r=1.2 w=1, msi r=1.5 w=1, spatial r=0.8 w=0.5 | 10 | 0.37 | 0.80 |
| RNA + MSI + smoothed RNA | rna r=1 w=1, msi r=1 w=0.5, rna_smoothed r=0.7 w=0.5, spatial r=1 w=0.25 | 9 | 0.47 | 0.92 |

Figures:
- `spatialleiden_rna_msi.pdf/.png`: cluster maps per section and run.
- `spatialleiden_rna_msi_with_features.pdf/.png`: the same maps, followed by dopamine, *Pcp4*, MSN1 and the region annotation.

Each cluster gets the colour of the region it overlaps most. If several clusters share a region, the largest keeps the region colour, so the intact and lesioned CP appear in different colours. Colours are comparable within a column (one clustering over all sections), not between columns.

## 7. Masking experiment, one section at a time

Question: what happens when one section has no MSI, and can the MSI of the other sections still help it, without imputing MSI values?

The MSI of **A1, B1 or C1** is hidden in turn. All three MALDI sections are then clustered together again (multi-sample); D1 is not used. Parameters as in 6c.

| method | layers | parameters |
|---|---|---|
| RNA only | RNA + spatial | `RUNS["RNA"]` |
| MSI only | MSI + spatial, real MSI of all sections | `RUNS["MSI"]` |
| Full MSI | RNA + MSI + spatial, real MSI of all sections | `RUNS["RNA + MSI"]` |
| Isolated | as Full MSI, but the masked section has no MSI edges | `RUNS["RNA + MSI"]` |
| Bridge: CCA / PLSCanonical / PLSSVD | as Isolated, plus an RNA bridge layer (see below) | `RUNS["RNA + MSI"]`, bridge layer r = 1.5, w = 1 (same as the MSI layer) |
| Bridge: ridge CCA | as above, with our own ridge-regularised CCA (reference) | as above |

*RNA only*, *MSI only* and *Full MSI* do not depend on the masked section; they are the same clusterings as in 6c. For *Isolated*, the MSI kNN graph is rebuilt from the two other sections only, with the PCA refitted on their spots, so the masked section does not influence the MSI embedding. With `log`, the only part of preprocessing that also sees the masked section is the `normalize_total` target, which is the same constant for every spot.

**RNA bridge.** Without its own MSI, the masked section is linked to the other sections only through RNA edges, and these do not carry the MSI structure. The bridge adds a layer built from the part of the RNA variation that is related to MSI. This is a supervised feature projection: no MSI values or MSI embeddings are predicted. The concept matches the bridge in Seurat bridge integration (Hao et al., *Nat Biotechnol* 2023).

1. Fit a projection between the RNA embedding (`X_pca_harmony`, 50 dimensions) and the MSI PCs (30, refitted without the masked section) on the paired spots of the **other** sections. Inputs are standardised; `bridge_n_components = 10` components. Three established methods from `sklearn.cross_decomposition`:
   - **CCA**: maximises the correlation of the RNA and MSI components.
   - **PLSCanonical**: maximises their covariance (symmetric PLS mode A).
   - **PLSSVD**: SVD of the cross-covariance of the standardised data, the "diagonal CCA" of Seurat v3 (Stuart et al., *Cell* 2019).
2. Project the RNA of **all** spots onto the RNA weights. Each component is weighted by its correlation on the training spots.
3. Build a kNN graph (k = 15) on this projection. It becomes the bridge layer, with edges between all sections.
4. Diagnostic only, not used for clustering: the correlation of the RNA and MSI component scores on the masked section, with its real MSI projected onto the MSI PCs fitted without it.

Held-out correlation on the masked section (ridge CCA; the three sklearn methods give practically the same values on PC space, because PCs are already decorrelated within a modality):

| masked section | train corr (C1–3) | held-out corr on the masked section (C1–3) |
|---|---|---|
| A1 | 0.65, 0.61, 0.34 | 0.65, 0.52, 0.16 |
| B1 | 0.66, 0.57, 0.33 | 0.65, 0.63, 0.25 |
| C1 | 0.62, 0.58, 0.29 | 0.72, 0.64, 0.26 |

The first two components generalise to the held-out section; the later ones carry little signal.

### 7a. Scores (masked section, annotated spots with MSI)

ARI lesion in CP:

| method | A1 | B1 | C1 |
|---|---|---|---|
| RNA only | 0.01 | 0.00 | 0.10 |
| MSI only | 0.80 | 0.57 | 0.70 |
| Full MSI | 0.70 | 0.85 | 0.88 |
| Isolated | 0.12 | 0.00 | 0.16 |
| Bridge: CCA | 0.38 | 0.30 | 0.41 |
| Bridge: PLSCanonical | 0.39 | 0.26 | 0.64 |
| Bridge: PLSSVD | 0.14 | 0.30 | 0.63 |
| Bridge: ridge CCA | 0.36 | 0.32 | 0.66 |

ARI region:

| method | A1 | B1 | C1 |
|---|---|---|---|
| RNA only | 0.49 | 0.45 | 0.44 |
| MSI only | 0.16 | 0.14 | 0.13 |
| Full MSI | 0.38 | 0.38 | 0.36 |
| Isolated | 0.55 | 0.47 | 0.46 |
| Bridge: CCA | 0.45 | 0.45 | 0.36 |
| Bridge: PLSCanonical | 0.46 | 0.45 | 0.35 |
| Bridge: PLSSVD | 0.50 | 0.43 | 0.36 |
| Bridge: ridge CCA | 0.48 | 0.39 | 0.40 |

- **Without the bridge:** a section needs its own MSI to separate intact from lesioned CP (Full MSI 0.70–0.88, Isolated 0.00–0.16).
- **With the bridge:** part of the MSI structure reaches the masked section (A1 and B1 about 0.3–0.4, C1 up to 0.64–0.66). *Full MSI* is not reached.
- **Seed variability:** with 5 Leiden seeds (and 5 instead of 10 components), the bridges on PCs gave 0.24–0.33 in all three sections (mean), *Isolated* 0.03–0.24; single C1 runs scattered between about 0.1 and 0.6, also for *Isolated*. The C1 values above are therefore an example run, not a robust effect size. The notebook keeps this setup switchable (`MASK_SEEDS`, `bridge_n_components`).
- **Bridges on genes × m/z** (optional, `FEATURE_BRIDGES`): unregularised CCA on 3,000 HVGs × 1,538 m/z overfits completely (training correlation 1.0, held-out about 0) and destroys the clustering; PLSCanonical and PLSSVD keep the regions but carry no lesion signal in their first components (dominated by high-variance features).

### 7b. Figures

- `masking_per_section_with_features.pdf/.png`: one row per masked section. Columns: RNA + spatial, MSI + spatial, RNA + MSI + spatial (full), RNA + MSI + spatial with the MSI of this section masked, the same plus the bridge of each of the three methods, then dopamine, MSN1 and the region annotation. The first three columns are the same clustering in every row; the masked columns are a separate clustering per row, so their colours are comparable only within one panel.
- `masking_A1_all_sections.png`, `masking_B1_all_sections.png`, `masking_C1_all_sections.png`: one figure per masked section, with **all three sections of the same clusterings** as rows (masked section labelled in red). Here colours are comparable down each column, so one can see whether the masked section's intact and lesioned CP fall into the same clusters as in the sections with MSI. The dopamine column of the masked section shows its real MSI, for comparison only.
- `rna_bridge_cca*.pdf/.png`: per bridge method, the correlation per component (training sections vs masked section) and the RNA and MSI scores of components 1 and 2 as scatter and as maps (component signs aligned across rows).
- `rna_bridge_methods_panels.pdf/.png`: held-out correlation per component, ARI lesion in CP and ARI region per method (points = seeds).

### 7c. RNA and MSI are complementary (C1)

`rna_msi_synergy_C1.pdf/.png`, PLSCanonical fitted on A1 + B1 and applied to C1 (held out):

| panel | content | AUC intact vs lesioned CP |
|---|---|---|
| a | UMAP of the RNA of C1 | – (mixed) |
| b | UMAP of the MSI of C1 | – (separated) |
| c | RNA PC1 vs PC2 | 0.60 (PC1) |
| d | RNA of C1 projected onto the PLSCanonical axes learned on A1 + B1 | 0.84 (component 1) |
| e | real MSI of C1 projected onto the same axes (reference) | 0.99 (component 1) |
| f | AUC summary | |

The lesion information is in the RNA, but not along its main axes; the MSI-guided projection finds it without any lesion labels and without C1's MSI. Panel e shows where the RNA projection should land: both are coordinates in the same shared space, learned from the paired spots of A1 + B1.

## 8. Parameter scan (`05_param_scan.ipynb`)

Descriptive scan for one section (`SCAN_SECTIONS`, default C1; A1 or B1 by changing one line). RNA and MSI embeddings and kNN graphs are rebuilt within the scanned section (single-sample). For each modality, `multiplex_leiden([modality, spatial])` is run over
- resolution of the modality layer: 0.1–3.0 (14 values; spatial resolution tied to it),
- spatial / modality weight ratio: 0, 0.1, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 2.

MSI + spatial is clustered on the spots with MSI only. No grid point is selected.

- `param_scan_maps_rna_C1.png`, `param_scan_maps_msi_C1.png`: cluster maps, rows = weight ratio, columns = resolution, plus the region annotation.
- `param_scan_heatmaps.png`: number of clusters and ARIs per grid point (for orientation).
- `param_scan_scores.csv`: all grid points.

Within one section, the two CPs are spatially separate, so a strong spatial layer alone can split intact from lesioned CP; the ARI lesion in CP of the scan is therefore not evidence for a modality. Results: `results/Cosie_Clustering/striatum/msi_log/param_scan_C1/`.

## 9. Outputs in `results/Cosie_Clustering/striatum/msi_log/`

| file | content |
|---|---|
| `mouse_brain_6OHDA_striatum.h5mu` | MuData: `rna`, `msi` (`layers["counts"]`, `layers["lognorm"]`), padded graphs, spatial graph |
| `rna_umap_harmony.png`, `msi_umap_harmony.png` | batch effect before and after Harmony |
| `dopamine_msi_pcp4_msn1.pdf/.png` | dopamine, *Pcp4*, MSN1 |
| `spatialleiden_21_clusters_region+local lesion.csv` | tuning grid (6b) |
| `spatialleiden_rna_msi.pdf/.png`, `spatialleiden_rna_msi_with_features.pdf/.png` | cluster maps (6c) |
| `spatialleiden_joint_no_integration.pdf`, `spatialleiden_mosaic_k9_normalized.pdf` | additional joint and mosaic runs including D1 |
| `masking_per_section_scores.csv`, `masking_per_section_with_features.pdf/.png`, `masking_*_all_sections.pdf/.png` | masking experiment (7) |
| `rna_bridge_cca*.pdf/.png`, `rna_bridge_methods_panels.pdf/.png` | bridge diagnostics (7b) |
| `rna_msi_synergy_C1.pdf/.png` | RNA and MSI are complementary (7c) |
| `param_scan_C1/` | parameter scan (8) |

## 10. Limitations

- **Single seed.** Most runs use one Leiden seed (42). Single runs can scatter strongly (see 7a); differences of a few hundredths of ARI are within seed variation.
- **Parameters chosen with the annotation.** The parameters in 6b and 6c were selected with the annotation in view, so they are an upper bound rather than an unbiased estimate.
- **Annotation-based evaluation.** The ARI against lesion side over the whole hemisphere is low for every run. The 6-OHDA lesion acts on the dopaminergic system, so it is visible mainly in CP.
- **Bridge parameters.** The bridge uses fixed parameters (10 components, layer r = 1.5 and w = 1, the same as the MSI layer). They were not tuned. The bridge can only carry MSI structure that is linearly related to RNA.
- **Clusters specific to one section.** A few MSI clusters come mostly from one section. Without Harmony, the MSI layer keeps a moderate section effect.
