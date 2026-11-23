# Project C: Derived modalities from gene expression (DMGE)

SpaceHack 2026 track C. Branch: `hack/c-derived-modalities`.

## Question

Can functional modalities derived from gene expression act as extra layers in multimodal SpatialLeiden? Examples are pathway activity, transcription-factor activity and cell-cell signalling. We want to know:

- What does each modality find on its own?
- How do the modalities differ?
- Which combination and weighting gives domains that gene expression alone misses?

## Dataset

Kuppe et al. 2022, *Spatial multi-omic map of human myocardial infarction*, Nature 608, 766–777. Visium sections from control, ischaemic and fibrotic hearts, plus snRNA-seq and snATAC-seq. The original study grouped spots with non-spatial clustering.

Readouts that are not used to build any modality, available for evaluation: annotated zones (control, remote, border, ischaemic, fibrotic), patient group and histology.

More datasets: TBD.

## Plan

| Step | Folder | What we do |
|---|---|---|
| 0 | `00_strategy/` | Decide which modalities to derive and with which tools |
| 1 | `01_data_preprocess/` | Download Kuppe et al. 2022; QC, normalization, RNA and spatial graphs |
| 2 | `02_modalities/` | Derive each modality (one subfolder each) |
| 3 | `03_single_modality_clustering/` | SpatialLeiden on each modality separately; compare the results |
| 4 | `04_multimodal_weight_grid/` | `spatialleiden_multimodal` on combinations of modalities across a grid of layer weights |
| 5 | `05_evaluation/` | Domain maps, ARI all-vs-all (redundancy), metrics vs. independent readouts (gain), an example of a functional-only domain |
| 6 | `06_metrics_exploration/` | New metrics for comparing domain maps |

Candidate modalities to start from:

| Modality | Tool |
|---|---|
| Pathway activity | decoupler + PROGENy |
| Transcription-factor activity | decoupler + CollecTRI |
| Cell-cell signalling | LIANA+ |

## Data format

Each modality is saved as its own AnnData file in `data/kuppe_2022/`: `rna.h5ad`, `progeny.h5ad`, `collectri.h5ad`, `liana.h5ad`, and so on.

- All files contain the same spots in the same order.
- `obs["sample"]` gives each spot's sample.
- Spot coordinates are in `obsm["spatial"]`.
- `rna.h5ad` also holds the spatial neighbour graph in `obsp["spatial_connectivities"]`.

For clustering, the files are combined into one MuData object, which is the input `spatialleiden_multimodal` takes.

## Folder structure

```
Derived_Modalities/
├── README.md                        this file
├── 00_strategy/                     which modalities and tools, and why
├── 01_data_preprocess/              download, QC, RNA + spatial graphs
├── 02_modalities/
│   ├── progeny/                     pathway activity
│   ├── collectri/                   TF activity
│   └── liana/                       cell-cell signalling
├── 03_single_modality_clustering/   each modality on its own, compared
├── 04_multimodal_weight_grid/       combinations x weight grid
├── 05_evaluation/
│   ├── domain_maps/
│   ├── redundancy_ari/
│   ├── gain_readouts/
│   └── functional_only_example/
├── 06_metrics_exploration/          new metrics
├── common/                          shared helper code
├── downloads/                       download script per dataset
├── envs/                            conda environment files
├── notes/                           meeting notes and decisions
├── data/                            data (not in git)
└── results/                         domains, metrics, figures (not in git)
```

Everything for project C lives in this folder. We do not change files outside it: the notebooks at the repository top level reproduce the preprint and depend on their current location.

## Working together

| Person | Focus |
|---|---|
| TBD | |
| TBD | |

- We both work directly on `hack/c-derived-modalities`. Pull before you start and before you push (`git pull --rebase`).
- Each person edits their own files, so we don't conflict on the same notebook.
- Name notebooks `<step>_<topic>.ipynb`, e.g. `02_progeny_kuppe.ipynb`.
- Clear notebook outputs before committing and run pre-commit (repository rules).
