# TRICARE 3D H&E → Multileiden

3D H&E pathology stacks (TRICARE) tiled with lazyslide → UNI2 features →
SpatialLeiden Multileiden, with OpenST-style 3D figures (slice FacetGrid, depth GIF, 3D scatter).

## Setup

```bash
conda env create -f TRICARE_3D/download_tricare.yml
conda env create -f TRICARE_3D/tricare_tile.yml
conda env create -f TRICARE_3D/tricare_univ2.yml
conda env create -f TRICARE_3D/tricare_multileiden.yml
conda env create -f TRICARE_3D/tricare_plot3d.yml

export TRICARE=./data/TRICARE
export MODEL_DIR=./models
# gated UNI2-h weights (request access on Hugging Face first):
conda activate tricare_univ2
python TRICARE_3D/download_univ2_weights.py --model-dir "$MODEL_DIR"
```

## Download + Multileiden

```bash
conda activate download_tricare
python TRICARE_3D/download_tricare.py --out-dir "$TRICARE" --cohorts prostate_test

# SLURM arrays (one specimen per task)
bash TRICARE_3D/submit_tricare_array.sh tile
bash TRICARE_3D/submit_tricare_array.sh univ2 --skip-existing

conda activate tricare_multileiden
python TRICARE_3D/run_tricare_univ2_multileiden.py \
  --tricare-dir "$TRICARE" --cohorts prostate_test --stage multileiden
```

## Figures

```bash
conda activate tricare_plot3d
python TRICARE_3D/plot_tricare_3d_spatialleiden.py \
  --tricare-dir "$TRICARE" --cohort prostate_test
```

Outputs: `$TRICARE/univ2_multileiden/figures/{cohort}/{specimen}/`
(`*_slices.pdf`, `*_depth.gif`, `*_3d.pdf`).
