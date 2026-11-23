# envs

Conda environment files for project C (e.g. `dmge.yml` for decoupler + LIANA+ + SpatialLeiden).

Write down here any version pins or fixes needed to make an environment work.

## `spagrn` (built 23 November 2026)

Built from `spagrn.yml`: Python 3.8, SpaGRN 1.1.4, pySCENIC 0.12.1, hotspotsc 1.1.1, scanpy 1.9.3, anndata 0.8.0, numpy 1.23.5, pandas 1.5.3, dask/distributed 2023.5.0 (2.3 GB). Jupyter kernel: "Python (spagrn)".

- **Python version:** SpaGRN 1.1.4 requires exactly Python 3.8. dask/distributed 2023.5.0 are the last releases supporting 3.8.

- **pip check:** it only reports that `descartes` (pulled in by pysal) is "not supported on this platform". This is harmless.

- **Databases:** the human cisTarget v10 databases are on the vault and readable from this env:

  - `hg38_*_full_tx_v10_clust.genes_vs_motifs.rankings.feather` (27,090 genes),
  - `motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl`,
  - `allTFs_hg38.txt` (1,892 TFs),

  all in one local folder, given to `run_spagrn.py` via the `CISTARGET_DIR` environment variable (download from resources.aertslab.org/cistarget).

- **Permissions:** the env is on the shared conda install. It was built with `umask 002` and made group-writable afterwards.

**`spagrn` fix (23 Nov 2026).** The first run failed at the Geary's C step with `TypeError: w must be a pysal weights object`.

- **Cause:** conda installed pysal 2.0.0, whose bundled `pysal.lib` weights class is not the `libpysal` 4.7 class that `esda` 2.5.1 checks for.
- **Fix:** one line in `spagrn/autocor.py`, `from pysal.lib import weights` → `from libpysal import weights`. The original is kept as `autocor.py.orig`. The command is in `spagrn.yml`.
- **Tested:** SpaGRN's own weights builder plus `gearys_c` now run.
- **Not used:** upgrading pysal to 2.7 is the alternative, but conda could not solve it in reasonable time on the shared install.
