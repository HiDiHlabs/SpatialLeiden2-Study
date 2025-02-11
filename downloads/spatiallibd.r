#!/usr/bin/env Rscript

# Settings
data_dir <- "../data"


# Script
library(spatialLIBD)
library(zellkonverter)


spe <- fetch_data("spe")

libd_dlpfc_dir <- file.path(data_dir, "LIBD_DLPFC")
dir.create(libd_dlpfc_dir, showWarnings = FALSE, recursive = TRUE)

coords <- as.data.frame(spatialCoords(spe))
colData(spe)$x <- coords$pxl_col_in_fullres
colData(spe)$y <- coords$pxl_row_in_fullres

writeH5AD(spe, file = file.path(libd_dlpfc_dir, "raw.h5ad"))
