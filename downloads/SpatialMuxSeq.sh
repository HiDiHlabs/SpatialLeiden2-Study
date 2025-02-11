#!/bin/bash

# download Spatial-Mux-seq data from the original publication
# Guo, Pengfei, et al. "Multiplexed spatial mapping of chromatin features, transcriptome and proteins in tissues." Nature Methods (2025): 1-10.
# https://doi.org/10.1038/s41592-024-02576-0

# from GEO GSE263333; https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE263333

OUT_DIR=../data/SpatialMuxSeq

URL="https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE263333&format=file"

mkdir -p $OUT_DIR

FILE="$OUT_DIR/GSE263333_RAW.tar"

wget -O $FILE $URL
tar -xf $FILE -C $OUT_DIR
