# TF activity (decoupler ULM, CollecTRI) on P1 and P9

Panels from `02_modalities/collectri/`. Method and the depth check are described in the README there.

| File | Content |
|---|---|
| `collectri_selected_tfs_zscore.png/.pdf` | 8 TFs (NFKB1, STAT3, SMAD3, HIF1A, SPI1, MEF2C, GATA4, TEAD1), z-scored within each slide |
| `collectri_selected_tfs.png` | The same TFs as raw ULM scores. These are mostly positive on every slide, so they show little within-slide contrast |
| `collectri_qc.csv` | Depth check per slide |

**Key point:** a typical TF has ≥ 5 targets detected in 16 % of P1 spots, against 70 % on P9 rep2. P1's TF layer is noisy mainly because of depth, a high mitochondrial fraction and myocardium-dominated spots. Its salt-and-pepper maps do not by themselves mean the tissue lacks structure.
