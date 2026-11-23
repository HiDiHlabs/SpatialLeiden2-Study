# Handoff for 24 November 2026: start here

Track C (derived modalities), branch `hack/c-derived-modalities`. Paths are relative to `Derived_Modalities/`. This adds to `HANDOFF_2026-11-23_resume_on_new_vm.md` (setup, data, what was done on 23 Nov). Read that one for the details; this one is the short version plus the last results of the day.

## 1. Pick up on the new VM in four steps

1. **Code:** `git clone` the study repo, then `git checkout hack/c-derived-modalities`.
2. **Objects:** download the OneDrive folder `SpatialLeiden_hackathon/derived_modalities_2026-11-23/` (Charité OneDrive). It contains:
   - `mudata/kuppe_<run>_derived.h5mu`: **work from these.** One per slide, with every inferred matrix and our labels (section 2).
   - `modalities/`: the same layers as separate files, plus the CCC table of the P9 remote zone (ACH0023).
   - `bundles/`: two tar files with `SHA256SUMS`. Unpacking them inside `Derived_Modalities/` restores `data/kuppe_2022/` (the Visium h5ads and spot tables the scripts need) and `results/kuppe_2022/spagrn/`.
3. **Environments:** `spatialleiden` (from `../spatialleiden.yml`) for everything. `envs/sccellfie.yml` and `envs/spagrn.yml` are needed only to recompute those layers, which you don't need to do.
4. **Results:** the cluster labels and summary tables of all runs are already in git under `results/kuppe_2022/`. The figure scripts run on them once `data/` is restored.

## 2. What is in each MuData

```python
import mudata as md
m = md.read_h5mu("kuppe_ACH0013_derived.h5mu")   # 10X001 = P1 control, ACH0012 = P9 IZ, ACH0013 = P9 IZ rep2
m["rna"].layers["counts"]   # raw counts; X = log1p CP10k; obsm X_pca (30); obsp connectivities (15-NN)
m["tf"]                     # CollecTRI ULM TF activity (~680 TFs), layers["padj"]
m["ccc"]                    # team CCC ligand-receptor local scores (pre-smoothed!)
m["met"]                    # scCellFie metabolic task activity (~210 tasks)
m["spagrn"]                 # SpaGRN regulon AUCell (314-355 regulons); uns["regulons"] is JSON {regulon: targets}
m.obs                       # authors' niches, cell2location proportions, QC + our labels sl_* (+ vessel_k2/k3 on ACH0013)
m.obsm["spatial"]           # pixel coordinates (full resolution)
```

**Our label columns** (seed 0, K = the authors' cell-type niche count, so 8 / 8 / 9 domains):

| Column | Layers | Spatial weight |
|---|---|---|
| `sl_rna`, `sl_tf`, `sl_ccc`, `sl_met` | one layer each | 0 |
| `sl_rna_space`, `sl_tf_space`, `sl_met_space` | one layer each | 0.75 |
| `sl_ccc_space025`, `sl_ccc_space075` | CCC | 0.25 / 0.75 |
| `sl_rna_tf_space` | RNA 1 + TF 0.5 | 0.75 (the run behind the vessel findings) |
| `sl_rna_ccc_space025`, `sl_rna_ccc_space075` | RNA 1 + CCC 0.5 | 0.25 / 0.75 |
| `sl_rna_met_space`, `sl_rna_all_space` | RNA 1 + metabolic 0.5 / RNA 1 + all three at 0.5 | 0.75 |

All five seeds of every run are in `results/kuppe_2022/**/labels_*.csv.gz`.

## 3. Last result of 23 Nov: SpaGRN vs CollecTRI

**Code and outputs:** script `02_modalities/spagrn/compare_spagrn_collectri.py`; outputs in `results/kuppe_2022/spagrn/compare_collectri/` (`summary.csv`, `per_tf_agreement.csv`, `vessel_band_top_tfs.csv`, `fig_spagrn_vs_collectri.png`).

**The two methods do different things:**

- **CollecTRI** uses the same curated regulons on every slide. Its scores are signed ULM t-values.
- **SpaGRN** infers regulons per slide. It keeps spatially variable genes, links each TF to the genes co-expressed with it, prunes those links by motif, and keeps activating regulons only. Its scores are AUCell.

| | P1 control | P9 IZ | P9 IZ rep2 |
|---|---|---|---|
| SpaGRN regulons | 314 | 324 | 355 |
| TFs shared with CollecTRI | 196 | 227 | 240 |
| Median regulon size (targets) | 15 | 28 | 33 |
| Median Spearman, same TF (null: mismatched TF) | 0.04 (0.03) | 0.11 (0.03) | 0.14 (0.08) |
| Share of shared TFs with Spearman > 0.3 | 2 % | 23 % | 33 % |
| Median target overlap, Jaccard (null) | 0.000 (0) | 0.002 (0) | 0.002 (0) |

**Findings:**

- **The target sets barely overlap** (median Jaccard ≈ 0). SpaGRN's targets are slide-specific co-expression modules, not curated regulons. The two methods share TF names, not regulons.
- **Spot-wise agreement is weak overall but real on the deep P9 slides,** and absent on shallow P1, which matches the depth confound seen before.
- **Where they agree, it is cell-identity TFs:**
  - TBX5 0.79, GATA4 / GATA5, MEF2A, and MEF2C (0.70 on P9 IZ) for cardiomyocytes;
  - SPI1 0.63–0.68 for myeloid cells;
  - SRF 0.66 on P9 rep2;
  - TWIST1 and FOXC1.
- **P9 rep2 artery band: the two methods are complementary.**

  | Method | Top TFs in the band | What they mark |
  |---|---|---|
  | CollecTRI | MYOCD, TCF21, SRF | Smooth muscle and fibroblasts |
  | SpaGRN | SOX17, FOXC1, FOXL1, then ATF3 / FOS / FOSB / JUNB | Arterial endothelium (SOX17, FOXC1); the second group is immediate-early stress genes |

  - **MYOCD and TCF21 are not SpaGRN regulons on any slide.**
  - **SpaGRN's SRF regulon follows the myocardium, not the artery** (band z = −0.17). It is a co-expression module around SRF in cardiomyocytes, not SRF's smooth-muscle programme.
- **165 regulons recur on all three slides,** and 153 of them are also CollecTRI TFs. These are the candidate set if SpaGRN becomes a SpatialLeiden layer.
- **Caveats:**
  - Some regulons are degenerate (for example STAT3 on P9 rep2 is nearly all zero).
  - AUCell (≥ 0, rank-based) and ULM (signed) are on different scales, so compare ranks.

## 4. Suggested plan for 24 Nov

1. **SpaGRN as a SpatialLeiden layer:** use the AUCell matrix, z-scored, PCA, 15-NN, exactly like TF. Compare `sl_spagrn_space` with TF + space and RNA + space, and run the vessel-recovery check (`04_multimodal_weight_grid/proto_rna_tf/vessel_band_across_modalities.py`). SpaGRN finds the arterial endothelium (SOX17), so it may outline the artery on its own.
2. **Feature-selection grid for TF and CCC** (from the 23 Nov feedback), to fit into an overnight run:
   - TF: filter on targets detected, or unscaled ULM scores;
   - CCC: top N pairs, or collapse pairs by ligand;
   - HVGs: 1,000 / 2,000 / 5,000;
   - PCs per layer;
   - extra-layer weights 0.25 / 0.5 / 1 / 2.
3. **CCC from the colleague:** method, kernel bandwidth, pair filter. They set how much space the CCC layer already contains (we use spatial weight 0.25 because of it).
4. **Slides and figures:**
   - the methods-slide brief is `handoff/HANDOFF_2026-11-23_methods_slide_rna_tf_ccc.md`;
   - key figures: the multipanels in `results/kuppe_2022/mm_configs/figures/`, the Kuppe Fig. 3 comparison, the vessel panels, and the SpaGRN figure above.
5. **Points to keep straight on slides:**
   - **The artery is also in the authors' molecular niches** (niche 11; their Fig. 3b). Our gain is a continuous domain and the wall-vs-margin split, not a new structure.
   - **The margin around the wall is not a fibrotic sheath.** It has the same collagen markers as the surrounding tissue.
   - **SIS:** TF, CCC and metabolic layers add no information beyond RNA for niches and composition on P9. RNA + CCC gains only for neighbourhood composition, which comes from CCC's built-in smoothing.
