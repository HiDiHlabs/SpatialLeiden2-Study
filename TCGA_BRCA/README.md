# TCGA-BRCA (DeepSpot-M → MultiLeiden)

SpatialLeiden2-style multimodal clustering on DeepSpot-M virtual ST for TCGA-BRCA:
`MuData({"rna", "HE"})` + `spatialleiden.spatialleiden_multimodal`.

## Setup

```bash
conda env create -f TCGA_BRCA/download_deepspotm_tcga.yml
conda env create -f TCGA_BRCA/lazyslide.yml
conda env create -f TCGA_BRCA/univ2.yml
conda activate download_deepspotm_tcga
huggingface-cli login   # accept terms for ratschlab/TCGA_virtual_spatial_transcriptomics_atlas
# UNI2 weights: directory containing pytorch_model.bin
export MODEL_DIR=./models
```

## Run (SLURM)

```bash
export OUT_DIR=./data/TCGA_BRCA
export MODEL_DIR=./models
export MAX_PATIENTS=100
bash TCGA_BRCA/submit_deepspotm_tcga_brca_slurm.sh all
# Debug one sample: bash TCGA_BRCA/submit_deepspotm_tcga_brca_slurm.sh test
```

Pipeline per HE array task: fetch SVS → lazyslide tiles (`non_overlap` default) → UNI2 → `he_embeddings/{id}.npy` → MuData.

## Cluster

```bash
python TCGA_BRCA/spatialleiden_mudata_run.py \
  --input "$OUT_DIR/mudata/SLIDE.h5mu" \
  --output "$OUT_DIR/spatialleiden/SLIDE.h5mu" \
  --layer-weights rna:1.5,HE:0.5,spatial:2
```

## Notes

- Clinical: TCGA-CDR OS + PFI; default stages II/III; primary diagnostic DX slides.
- Gene filter off by default; optional `--pcc-cameo` / `--pcc-filter`.
- Envs: `download_deepspotm_tcga` (select/fetch/mudata), `lazyslide` (tile), `univ2` (HE features).
