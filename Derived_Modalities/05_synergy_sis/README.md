# Synergistic Information Score (SIS)

Richter et al. 2026, *Multimodal cell foundation models require synergistic integration beyond alignment*,
Nat. Mach. Intell. (doi 10.1038/s42256-026-01328-3). SIS = (P(fused) - max P(single)) / max P(single) with linear
probes (eq. 13).

`sis_layers.py` (env spatialleiden): RNA PCA vs TF / CCC / metabolic PCA vs their concatenation; targets
authors' cell-type and molecular niches, cell2location composition, composition 2-3 rings away; 5-fold
spatial-block CV; null = concatenation with the second layer shuffled. Output `results/kuppe_2022/sis/sis_layers.csv`.

Result: on P9 the derived layers are redundant with RNA (SIS within the null range); RNA + CCC gains only for
neighbourhood composition, which reflects CCC's built-in spatial smoothing. SpatialLeiden fuses graphs, not
embeddings, so this tests the layers, not our clusterings.
