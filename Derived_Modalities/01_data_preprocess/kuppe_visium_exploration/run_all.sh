#!/bin/bash
# Kuppe 2022 Visium exploration: extract tables, make figures, write inventory.
# Run from Derived_Modalities/. Envs: hbc-scrnaseq (Seurat), masld_atlas (Python with openpyxl).
set -e
D=data/kuppe_2022; O=results/kuppe_2022/visium_exploration; S=01_data_preprocess/kuppe_visium_exploration
T=${TMPDIR:-/tmp}/kuppe_rds_tmp; mkdir -p $T $O
R=~/miniforge3/envs/hbc-scrnaseq/bin/Rscript; PY=~/miniforge3/envs/masld_atlas/bin/python
$R $S/01_extract_seurat_tables.R $D/visium_raw $D/visium_tables $T
$PY -I $S/02_extract_spaceranger_files.py $D/visium_raw $D/visium_tables
$PY -I -W ignore $S/03_figures.py $D $O > $O/03_figures.log
$PY -I -W ignore $S/04_inventory.py $D $O
$PY -I $S/05_report_tables.py $O $S/report
(cd $S/report && ~/miniforge3/envs/latex/bin/tectonic kuppe_visium_exploration.tex)
cp $S/report/kuppe_visium_exploration.pdf $S/NOTES.md $O/
cp $O/data_inventory.csv $S/data_inventory.csv
