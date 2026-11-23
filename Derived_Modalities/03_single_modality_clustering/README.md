# 03_single_modality_clustering

Run SpatialLeiden on each modality separately (that modality's graph + the spatial graph), including `rna` as the baseline.

- Use the same resolution search for every modality, so they all produce the same number of domains.
- Compare the domain maps side by side.
- Compute the ARI between every pair of modalities.
- Note where the modalities agree and where they split or merge domains differently.

Output: `../results/kuppe_2022/single/domains_<modality>.parquet`, plus a comparison notebook.

This step answers: what does each modality see on its own?

Start here: `00_intro_spatialleiden_kuppe.ipynb`. Sections 0–6 introduce SpatialLeiden on one Kuppe slide (GT/IZ_P9_rep2): Leiden vs SpatialLeiden at three `layer_ratio` values. Sections 7–14 are experiment 0: a `layer_ratio` sweep (8 values × 5 seeds) on P9, its sister section and the P1 control, with four choices of the number of domains K. Design, scripts and `RESULT.md` are in `exp0_layer_ratio_sweep/`.
