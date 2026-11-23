# figures

Finished figures for project C, one subfolder per figure. Each subfolder holds the figure (PDF for print, PNG for slides), its caption and the table behind it.

PDFs are normally ignored by git in this repository. Figures here are added on purpose (`git add -f`).

| Figure | Content |
|---|---|
| `kuppe_cohort/` | Kuppe et al. 2022 cohort: patients, zones, data types, time after symptom onset |
| `kuppe_fig2b/` | Reproduction of Kuppe et al. 2022 Fig. 2b: cell-type niches on P1 (control) and P9 (IZ) |
| `exp0_layer_ratio/` | Experiment 0: SpatialLeiden `layer_ratio` sweep on P9 rep2, P9 and P1, for four choices of the number of domains K |
| `collectri_tf_activity/` | TF activity maps (decoupler ULM + CollecTRI) on P1 and P9, with the depth check |

Not in git: the slide atlas (`results/kuppe_2022/slide_atlas/kuppe_visium_slide_atlas.pdf`, 41 MB). Rebuild it with `01_data_preprocess/kuppe_slide_atlas/slide_atlas.py` (about 90 s).
