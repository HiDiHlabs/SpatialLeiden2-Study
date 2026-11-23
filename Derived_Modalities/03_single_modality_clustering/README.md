# 03_single_modality_clustering

Run SpatialLeiden on each modality separately (that modality's graph + the spatial graph), including `rna` as the baseline.

- Use the same resolution search for every modality, so they all produce the same number of domains.
- Compare the domain maps side by side.
- Compute the ARI between every pair of modalities.
- Note where the modalities agree and where they split or merge domains differently.

Output: `../results/kuppe_2022/single/domains_<modality>.parquet`, plus a comparison notebook.

This step answers: what does each modality see on its own?
