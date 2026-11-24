# kuppe_fig1h

Reproduction of Kuppe et al. 2022, Fig. 1h–j on four slides. This is the base for the Fig. 1h + SpatialLeiden deliverable; the SpatialLeiden columns come next.

| File | Content |
|---|---|
| `kuppe_fig1h.{pdf,png}` | Panel h only (89 mm) |
| `kuppe_fig1hij.{pdf,png}` | Panels h, i and j (183 mm) |
| `fig1hij_summary.json` | Spot counts and the 1st and 99th percentile used for each map's colour scale |

| Row | Slide | Visium run |
|---|---|---|
| Ctrl tissue (P7) | control_P7 | ACH003 |
| BZ tissue (P3) | RZ/BZ_P3 | 10X0026 |
| IZ tissue (P3) | IZ_P3 | 10X0017 |
| FZ tissue (P14) | FZ_P14 | ACH005 |

## How each panel relates to the paper

The authors' repo (saezlab/visium_heart) has no script for Fig. 1h–j. The panel was styled after plotting, so the choices below were checked against the published panel.

| Panel | Our values | Paper | Same values? |
|---|---|---|---|
| h, cell2location | Authors' `c2l_props` (`prop_<celltype>` in `visium_tables/<run>_spots.csv`) | Same | Yes |
| i, PROGENy (TGFβ, NFκB) | Recomputed with the authors' formula: top-500 human model × log-normalised RNA, z-scored per slide (`02_modalities/progeny/progeny_kuppe.py`) | Their `progeny` assay, in Seurat objects not on this VM | Same method; should be near-identical |
| j, TF | CollecTRI ULM from RNA (`modalities/collectri_all27/`), our SpatialLeiden TF layer | *TF binding* activity: HINT footprinting on snATAC, mapped to spots through cell-type composition | **No: a stand-in.** Panel title changed to "TF activity (CollecTRI)" |

Further details:

- **Panel h uses proportions, not absolute abundance (`c2l`).** Abundance does not reproduce the IZ P3 and FZ P14 vCM patterns; proportions do.
- **FOS–JUNB:** in the paper this is one heterodimer motif. Here it is the mean of the CollecTRI FOS and JUNB scores.
- **Panel j agreement:** SRF in P7 (vessel) and FOS–JUNB in P14 (top right) agree with the paper. MEF2C in P14 and the IZ panels are more diffuse than in the paper, as expected from a method based on RNA rather than footprinting.
- **Colour scale:** per map, clipped at the 1st and 99th percentile (Seurat `min.cutoff = "q1"`, `max.cutoff = "q99"`, as in the authors' `st_snRNAseq/07_structure/plot_spatial_ex_panels.R`).
- **Palettes:** inferno for h. The palettes for i (blue → magenta → yellow) and j (navy → cyan → yellow) are matched by eye.
- **Black background, dashed outline, arrows:** transcribed by eye. Arrow positions are approximate.
- **Font:** Arial is not installed on this VM, so matplotlib falls back to DejaVu Sans.

## How to rebuild

From `Derived_Modalities/`:

```
# PROGENy scores (decoupler env; fetches the PROGENy model once)
python 02_modalities/progeny/progeny_kuppe.py data/kuppe_2022
# figures (spatialleiden env)
python 01_data_preprocess/kuppe_fig1/fig1hij.py data/kuppe_2022 figures/kuppe_fig1h
```

The PROGENy outputs (`data/kuppe_2022/modalities/progeny/`) are not in git, like the rest of `data/`.
