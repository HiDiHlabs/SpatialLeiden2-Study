# Multimodal SpatialLeiden prototypes (RNA + TF / CCC / metabolic)

All runs: K = authors' cell-type niche count per slide, 5 seeds, slides 10X001 / ACH0012 / ACH0013. Graphs are
built by the experiment-0 functions (`03_.../exp0_sweep.py`). Run from `Derived_Modalities/`, env spatialleiden.

| Script | What | Output (`results/kuppe_2022/`) |
|---|---|---|
| `mm_proto.py`, `mm_proto_figures.py` | RNA + TF weight grid (spatial 0.5/0.75 x TF 0-2) | `mm_proto_rna_tf/` |
| `mm_configs.py` | named configs rna_tf / rna_ccc / rna_met / rna_all (RNA 1, extra 0.5, spatial 0.75); env vars `MM_SPATIAL_W`, `MM_CONFIGS` | `mm_configs/`, `mm_configs_sp0.25/` |
| `fig_rna_tf_combo.py` | multipanel RNA / TF / RNA + TF (+ SPP1, NPPB) | `mm_proto_rna_tf/figures/` |
| `fig_modality_combo.py {ccc,ccc_sp025,met,all}` | multipanels for CCC, metabolic, all modalities | `mm_configs/figures/` |
| `fig_vessel_band.py`, `vessel_band_subclusters.py`, `fig_vessel_band_subclusters.py` | P9 rep2 artery domain, nested K = 2/3 | `mm_proto_rna_tf/figures/`, `mm_proto_rna_tf/vessel_band/` |
| `vessel_band_across_modalities.py` | does each run recover the artery | `mm_proto_rna_tf/vessel_band/` |
| `fig_kuppe_fig3_compare.py` | Kuppe Fig. 3 redrawn: authors' molecular niches vs our domains | `kuppe_fig3_compare/` |

Interpretation and numbers: `handoff/HANDOFF_2026-11-23_resume_on_new_vm.md` (section 4).
