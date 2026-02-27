#!/bin/bash

OUT_DIR=../data/VisiumHD_CRC

URL_LIST=(
    "https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Colon_Cancer/Visium_HD_Human_Colon_Cancer_tissue_image.btf"
    "https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Colon_Cancer/Visium_HD_Human_Colon_Cancer_feature_slice.h5"
    "https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Colon_Cancer/Visium_HD_Human_Colon_Cancer_spatial.tar.gz"
    "https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Colon_Cancer/Visium_HD_Human_Colon_Cancer_binned_outputs.tar.gz"
)

mkdir -p $OUT_DIR

echo ${URL_LIST[@]} | xargs -n 1 -P 8 wget -q -P $OUT_DIR

tar -xvzf $OUT_DIR/Visium_HD_Human_Colon_Cancer_binned_outputs.tar.gz
tar -xvzf $OUT_DIR/Visium_HD_Human_Colon_Cancer_spatial.tar.gz
