# SpatialLeiden2-Study

This repository contains all code to reproduce the results from the
Multi-sample and Multi-modal integration study using SpatialLeiden

TODO: add link to preprint/publication

## Contribution

When contributing to this repository please ensure that

- before you commit a notebook (ipynb) you remove all outputs to ensure that the file sizes and commits don't explode
- install [pre-commit](https://pre-commit.com/) and run it when committing. This can be achieved with a git hook.

Just follow the [documentation](https://pre-commit.com/#3-install-the-git-hook-scripts).

## Environments

Conda environment configuration files (all .yml files) are available where needed.
The corresponding environment should mostly be clear from the (file)name or
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
Note that some of them need their own environment.

## Analysis

Analysis can be reproduced by running the Jupyter Notebooks in the top-level
directory.
The notebooks analysing H&E datasets with `lazyslide` (VisiumHD-CRC) will require to be run in the 'lazyslide' conda environment,
all others can be run in the 'spatialleiden' environment.

### `time`

`time` (v1.9) was installed via

```bash
mkdir -p ~/local/src
cd ~/local/src

wget https://ftp.gnu.org/gnu/time/time-1.9.tar.gz
tar -xzf time-1.9.tar.gz
cd time-1.9

./configure --prefix=$HOME/local
make
make install
```
