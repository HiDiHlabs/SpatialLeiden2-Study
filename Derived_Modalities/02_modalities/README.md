# 02_modalities

One subfolder per modality chosen in `00_strategy`. Each derives its modality from the same normalized expression and builds a k-NN graph in its own feature space.

| Subfolder | Modality | Tool |
|---|---|---|
| `progeny/` | Pathway activity | decoupler + PROGENy |
| `collectri/` | Transcription-factor activity | decoupler + CollecTRI |
| `liana/` | Cell-cell signalling | LIANA+ (spatial bivariate scores) |

To add a modality, create a folder with the same layout.

Output per modality: `../data/kuppe_2022/<modality>.h5ad`, with spots in the same order as `rna.h5ad`. Every parameter (method, resource version, minimum targets, spatial bandwidth) is stored in `uns`.
