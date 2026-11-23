# 01_data_preprocess

Load the Kuppe et al. 2022 Visium data, run QC and normalization, and build the graphs that every modality shares.

- Download script: `../downloads/kuppe_2022.py`, which writes to `../data/kuppe_2022/`.
- Spatial neighbour graph: `squidpy.gr.spatial_neighbors(coord_type="grid", n_neighs=6)`, as in the preprint for Visium.
- RNA k-NN graph: same preprocessing as the preprint (2,000 HVGs, PCA with 30 components, 15 neighbours).

Output: `../data/kuppe_2022/rna.h5ad`. This becomes the `rna` modality and defines the spot order used by every other modality.
