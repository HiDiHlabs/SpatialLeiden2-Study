#!/bin/bash

OUT_DIR=../data/MOSTA_temporal

URL_LIST=(
    "https://ftp.cngb.org/pub/SciRAID/stomics/STDS0000058/stomics/Mouse_embryo_all_stage.h5ad"
)

mkdir -p $OUT_DIR

echo ${URL_LIST[@]} | xargs -n 1 -P 8 wget -q -P $OUT_DIR
