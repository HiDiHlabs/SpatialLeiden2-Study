# Brief for a methods slide: SpatialLeiden on RNA, TF activity and CCC

**For Claude in PowerPoint.** Build **one** methods slide (16:9) from this brief. The audience is a hackathon / lab meeting of computational biologists who know scRNA-seq and Leiden clustering. Use only the facts below and do not add numbers. The figures shown next to this slide are the multipanels `fig_rna_tf_combo` (RNA and TF) and `fig_ccc_sp025_combo` (CCC), so this slide explains how they were made.

## Slide title

SpatialLeiden on derived modalities: RNA, TF activity and cell-cell communication

## Data (one line on the slide)

Kuppe et al. 2022 myocardial infarction Visium: P1 control (10X001), P9 ischaemic zone (ACH0012) and P9 ischaemic zone rep2 (ACH0013, a separate piece of tissue, mostly granulation tissue).

## Core idea (one sentence)

SpatialLeiden clusters a multiplex graph. Each modality contributes a kNN graph built from its own features, and a spatial grid graph is added with its own weight. Every layer gives the same spots in a different feature space.

## Layout suggestion

Three columns, one per modality, each a short vertical pipeline from features to graph. Under them sits a shared band, "+ spatial graph → SpatialLeiden → domains". Settings go in a small table at the bottom.

## Column 1: RNA

- Visium counts (spots × genes)
- Genes in ≥ 10 spots; **2,000 highly variable genes** (seurat_v3, per slide)
- Normalise, log1p, scale
- PCA, 30 components
- 15-nearest-neighbour graph

## Column 2: TF activity

- Same RNA counts, normalised
- **decoupler ULM + CollecTRI** regulons (TFs with ≥ 5 targets): **~680 TFs per spot**
- z-score per TF
- PCA, 30 components
- 15-nearest-neighbour graph

## Column 3: Cell-cell communication (CCC)

- **Ligand-receptor local scores per spot** (LIANA+-style, values 0–1), computed by a team member
- 154 / 342 / 620 pairs per slide (P1 / P9 IZ / P9 rep2)
- z-score per pair
- PCA, 30 components
- 15-nearest-neighbour graph
- ⚠ **The scores are already spatially smoothed**, so adding the spatial graph partly counts space twice; **CCC therefore gets a lower spatial weight (0.25)**

## Shared band: spatial graph and clustering

- **Spatial graph:** Visium hexagonal grid, 6 neighbours per spot
- **SpatialLeiden 0.4** (multiplex Leiden): quality = Σ weightₘ × modularityₘ over layers
- **Number of domains K fixed** to the authors' cell-type niche count (8 / 8 / 9) by a resolution search, so all maps are comparable
- **Seed 0** shown; 5 seeds run for stability

## Settings table

| Map in the figure | RNA weight | TF weight | CCC weight | Spatial weight |
|---|---|---|---|---|
| RNA only | 1 | – | – | 0 |
| RNA + space | 1 | – | – | 0.75 |
| TF only / TF + space | – | 1 | – | 0 / 0.75 |
| CCC only / CCC + space | – | – | 1 | 0 / **0.25** |
| RNA + TF + space | 1 | 0.5 | – | 0.75 |
| RNA + CCC + space | 1 | – | 0.5 | **0.25** |

Footer-size notes:

- **Spatial weight 0.75** comes from our layer-ratio sweep (0–2.5) and fits the preprint's ~0.8 for single Visium slides.
- **CCC uses spatial weight 0.25, not 0.75.** About a quarter of each spot's CCC feature-space neighbours are also physical neighbours (vs 4–8 % for RNA and TF), so CCC already carries spatial smoothing. At 0.25, CCC + space is about as smooth as RNA + space at 0.75. RNA + CCC + space at 0.25 agrees with RNA + space as well as at 0.75, without the oversmoothed blobs.
- **The extra-layer weight 0.5 is our choice.** The preprint gives no rule for it.
- **Weights are ratios:** only their proportions matter, so RNA is fixed at 1 as the reference.

## How the maps are read (small text or speaker notes)

- **Domain colours** are matched to the RNA + space map (best overlap), so the same colour means the corresponding domain.
- **Spatial coherence** = the share of neighbouring spots that fall in the same domain.

## Caveats (speaker notes, not the slide body)

- **P1 is shallow** (~3,300 UMIs/spot): a typical TF has ≥ 5 targets detected in only 16 % of P1 spots, against 70 % on P9 rep2. P1's TF layer is therefore noisy.
- **No feature selection was applied to TFs or CCC pairs** (all ~680 TFs and all pairs, z-scored). A filter is planned.
- **The CCC filtering and kernel bandwidth** used upstream are not yet documented.

## Style

- Clean, minimal design; one accent colour per modality (RNA, TF, CCC) used consistently, and grey for the spatial layer.
- No logos, and no figures from the published paper.
