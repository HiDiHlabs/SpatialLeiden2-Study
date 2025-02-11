# SpatialLeiden2-Study

This repository contains all code to reproduce the results from the
Multi-sample and Multi-omics integration study using SpatialLeiden

TODO: add link to preprint/publication

## Environments

Conda environment configuration files (all .yml files) are available where needed.
The corresponding environment should be clear from the (file)name or
might be directly used within the code.
Environments must be installed manually before executing the code i.e. via

```bash
conda env create -f {env}.yml
```

Furthermore, the `CONDA_PATH` in [utils.py](./utils.py) may need to be
adjusted according to your systems setup before running the notebooks.

TODO: add environment dumps with exact versions

## Data

All data necessary to reproduce the results can be downloaded by executing
the scripts in the [downloads](./downloads) directory.
Note that some of them may need their own environment.

## Analysis

Analysis can be reproduced by running the Jupyter Notebooks in the top-level
directory.
The notebooks analysing multi-modal datasets (SpatialGlue Simulation and
Spatial-Mux-seq) will require to be run in the 'multimodal' conda environment,
all others can be run in the 'multisample' environment.
