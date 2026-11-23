# 04_multimodal_weight_grid

Run `spatialleiden_multimodal` on combinations of modalities across a grid of layer weights.

The grid is defined in `grid.yaml`, which is the only place to change it:

- **Modality combinations**: e.g. rna + one modality, rna + all modalities, all functional modalities without rna.
- **Weight per modality**: e.g. 0, 0.5, 1, 2. A weight of 0 switches a modality off.
- **Spatial weight**: e.g. 0.5, 1, 2, 2.5 (the preprint used 2 for Visium HD).
- **Resolution**: either fixed, or searched so that every run gives the same number of domains.

Write one row per run to `../results/kuppe_2022/grid/runs.parquet`: the parameters, the number of domains, modularity and run time. Domain labels go in `../results/kuppe_2022/grid/domains/<run_id>.parquet`.

Start with a coarse grid on one sample, then refine around the interesting region.
