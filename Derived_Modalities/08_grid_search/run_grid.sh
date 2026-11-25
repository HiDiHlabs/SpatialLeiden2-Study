#!/bin/bash
# Waits for the RNA + TF + CCC joint run, refreshes its figures, then runs the grid search (seed 0 first).
# Run from Derived_Modalities/ inside screen. Resumable: rerunning skips finished calibrations and runs.
set -u
R=results/kuppe_2022; G=$R/grid_search; PY=${PY:-python}; PD=${PD:-python}  # PY: env spatialleiden_joint; PD: an env with decoupler / liana
log() { echo "$(date -u +%H:%M) UTC $*"; }
log "waiting for exp3b"
while [ ! -f $R/EXP3B_EXIT ]; do sleep 60; done
log "exp3b exit $(cat $R/EXP3B_EXIT)"
$PD -I -W ignore 03_single_modality_clustering/exp2_joint_harmony/domain_enrichment_by_zone.py data/kuppe_2022 $R exp3b_joint_rna_tf_ccc r0.75_s0 > $R/exp3b_enrichment.out 2>&1
for v in "" v2; do $PY -W ignore 07_integration/fig_joint_umaps.py data/kuppe_2022 $R $v > $R/joint_umaps/fig$v.out 2>&1; done
log "figures refreshed"
df -BG --output=avail / | tail -1 | tr -dc 0-9 | { read f; [ "$f" -gt 40 ] || { log "less than 40 GB free, stop"; exit 1; }; }
$PY -W ignore 08_grid_search/grid_search.py data/kuppe_2022 $R $G calibrate 28 > $G/calibrate.out 2>&1; log "calibration done"
GRID_SEEDS=0 $PY -W ignore 08_grid_search/grid_search.py data/kuppe_2022 $R $G grid 28 > $G/grid_s0.out 2>&1; log "grid seed 0 done"
$PY -W ignore 08_grid_search/domain_table.py $G > $G/domain_table_s0.out 2>&1; log "domain table (seed 0) done"
GRID_SEEDS=1,2,3,4 $PY -W ignore 08_grid_search/grid_search.py data/kuppe_2022 $R $G grid 28 > $G/grid_s1-4.out 2>&1; log "grid seeds 1-4 done"
$PY -W ignore 08_grid_search/domain_table.py $G > $G/domain_table.out 2>&1; log "domain table (all seeds) done"
