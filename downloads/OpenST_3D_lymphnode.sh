#!/bin/bash

# download Open-ST 3D reconstructed metastatic lymph node from the original publication
# Schott, Marie, et al. "Open-ST: High-resolution spatial transcriptomics in 3D." Cell 187.15 (2024): 3953-3972.
# https://doi.org/10.1016/j.cell.2024.05.055

# from GEO GSE263498; https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE263498

OUT_DIR=../data/3D_lymphnode

URL="https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE251926&format=file&file=GSE251926%5Fmetastatic%5Flymph%5Fnode%5F3d%2Eh5ad%2Egz"

mkdir -p $OUT_DIR

FILE="$OUT_DIR/GSE251926_metastatic_lymph_node_3d.h5ad.gz"

wget -O $FILE $URL
gunzip $FILE
