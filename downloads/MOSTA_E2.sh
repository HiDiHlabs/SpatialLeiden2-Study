#!/bin/bash

OUT_DIR=../data/MOSTA

URL_LIST=(
    "https://ftp.cngb.org/pub/SciRAID/stomics/STDS0000058/stomics/E16.5_E2S1.MOSTA.h5ad"
    "https://ftp.cngb.org/pub/SciRAID/stomics/STDS0000058/stomics/E16.5_E2S2.MOSTA.h5ad"
    "https://ftp.cngb.org/pub/SciRAID/stomics/STDS0000058/stomics/E16.5_E2S3.MOSTA.h5ad"
    "https://ftp.cngb.org/pub/SciRAID/stomics/STDS0000058/stomics/E16.5_E2S4.MOSTA.h5ad"
    "https://ftp.cngb.org/pub/SciRAID/stomics/STDS0000058/stomics/E16.5_E2S5.MOSTA.h5ad"
    "https://ftp.cngb.org/pub/SciRAID/stomics/STDS0000058/stomics/E16.5_E2S6.MOSTA.h5ad"
    "https://ftp.cngb.org/pub/SciRAID/stomics/STDS0000058/stomics/E16.5_E2S7.MOSTA.h5ad"
    "https://ftp.cngb.org/pub/SciRAID/stomics/STDS0000058/stomics/E16.5_E2S8.MOSTA.h5ad"
    "https://ftp.cngb.org/pub/SciRAID/stomics/STDS0000058/stomics/E16.5_E2S9.MOSTA.h5ad"
    "https://ftp.cngb.org/pub/SciRAID/stomics/STDS0000058/stomics/E16.5_E2S10.MOSTA.h5ad"
    "https://ftp.cngb.org/pub/SciRAID/stomics/STDS0000058/stomics/E16.5_E2S11.MOSTA.h5ad"
    "https://ftp.cngb.org/pub/SciRAID/stomics/STDS0000058/stomics/E16.5_E2S12.MOSTA.h5ad"
    "https://ftp.cngb.org/pub/SciRAID/stomics/STDS0000058/stomics/E16.5_E2S13.MOSTA.h5ad"
)

mkdir -p $OUT_DIR

echo ${URL_LIST[@]} | xargs -n 1 -P 8 wget -q -P $OUT_DIR
