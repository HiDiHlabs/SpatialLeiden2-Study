#!/bin/bash

# download simulated data from SpatialGLue
# Long, Yahui, et al. "Deciphering spatial domains from spatial multi-omics with SpatialGlue." Nature Methods (2024): 1-10.
# https://doi.org/10.1038/s41592-024-02316-4

# from zenodo; https://doi.org/10.5281/zenodo.10362607

OUT_DIR=../data/SpatialGlue

URL="https://zenodo.org/records/10362607/files/Data_SpatialGlue.zip"

mkdir -p $OUT_DIR

wget -P $OUT_DIR $URL
unzip "$OUT_DIR/Data_SpatialGlue.zip" -d $OUT_DIR
