# CollecTRI transcription-factor activity

TF activity per Visium spot, from decoupler 2.x (`dc.mt.ulm`, univariate linear model) with the human CollecTRI network (`dc.op.collectri`: 1,185 TFs, 42,990 edges), `tmin = 5`.

- **Expression:** raw counts normalized with `normalize_total(1e4)` + `log1p`, all genes in the file.
- **Environment:** `masld_atlas` (decoupler 2.2.0), used read-only.
- **Method choice:** ULM rather than the paper's wmean/wsum. The paper never computed spot-level TF activity from CollecTRI, so there is nothing to match, and ULM is decoupler's standard method for TF activity.

| Script | Does |
|---|---|
| `collectri_tf_activity.py` | Builds the spot × TF matrix and the depth check, and draws raw-score maps |
| `collectri_maps_zscore.py` | Draws within-slide z-scored maps (raw ULM scores are mostly positive, which hides within-slide structure) |

**Outputs**

- `data/kuppe_2022/modalities/collectri_ulm_<run>.h5ad`: X = ULM score, `layers["padj"]`, spots in h5ad barcode order, `obsm["spatial"]`.
- `results/kuppe_2022/collectri/`: QC tables and figures.

## Depth check (23 Nov 2026)

| Slide | Median genes per spot | TFs scored | Median targets per TF | Targets detected per spot (median TF) | Spots with ≥ 5 targets detected (median TF) | Scores with padj < 0.05 |
|---|---|---|---|---|---|---|
| P1 control (10X001) | 1,649 | 681 | 18 | 3 | **16 %** | 15 % |
| P9 IZ (ACH0012) | 1,985 | 679 | 19 | 4 | 31 % | 23 % |
| P9 IZ rep2 (ACH0013) | 3,458 | 693 | 19 | 6 | **70 %** | 29 % |

673 TFs are scored on all three slides.

**Almost the same TFs pass the filter on every slide.** `tmin` only checks that targets exist in the gene list. What differs is the per-spot signal.

- On P1, a typical TF has fewer than 5 of its targets detected in 84 % of spots. On the deep P9 section, only 30 % of spots fall short.
- P1's TF layer is therefore much noisier. Its salt-and-pepper maps cannot, on their own, separate "no spatial structure" from "no signal". The same slide also has a high mitochondrial fraction (32 %).

**Recommendation:** use the TF layer on P9 for the first multimodal experiment, and treat P1's TF layer as exploratory.
