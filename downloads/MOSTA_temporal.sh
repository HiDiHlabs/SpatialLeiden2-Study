#!/bin/bash

OUT_DIR=../data/MOSTA_temporal


URL="https://ftp.cngb.org/pub/SciRAID/stomics/STDS0000058/stomics/Mouse_embryo_all_stage.h5ad"

mkdir -p $OUT_DIR

wget -P $OUT_DIR $URL
