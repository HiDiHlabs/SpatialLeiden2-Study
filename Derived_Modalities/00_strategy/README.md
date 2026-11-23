# 00_strategy

Decide which modalities to derive from the Kuppe et al. 2022 expression data, and with which tools.

For each candidate modality, record in `modalities.md`:

| Field | Example |
|---|---|
| Modality | Pathway activity |
| Tool and resource | decoupler + PROGENy |
| What biology it should capture | Signalling states, e.g. hypoxia or TGF-β in the infarct zones |
| Feature count | 14 pathways |
| Spatial or per-spot | Per-spot |
| Known caveats | Few features, so the k-NN graph can be coarse |
| Decision | Include / exclude / maybe, and why |

Starting candidates: pathway activity (PROGENy), transcription-factor activity (CollecTRI), cell-cell signalling (LIANA+). Others can be added, e.g. cell-type proportions from deconvolution using the Kuppe snRNA-seq reference.

Decisions and their reasons also go in `../notes/`.
