# Export the spot-registered MALDI-MSI data from the original study (Vicari et al. 2024)
# so it can be read in Python.
#
# misc/se.multi.list holds one Seurat object per section with RNA and MSI already mapped
# onto the Visium spots. The raw MSI csv files in sma/ are still in MSI pixel coordinates
# (rotated and shifted relative to Visium), so this is the registered version to use.
#
# Run with the r44_seurat env from the SpatialLeiden2-Study folder:
#   /dh-projects/ag-ishaque/analysis/malta/envs/envs/r44_seurat/bin/Rscript Cosie_Clustering/00_export_msi_spots.R

suppressMessages({
  library(Seurat)
  library(Matrix)
})

data_dir <- "data/01_mouse_brain_MALDI_visium"
out_dir <- file.path(data_dir, "processed", "msi_spots")

se_list <- readRDS(file.path(data_dir, "misc", "se.multi.list"))

for (section in names(se_list)) {
  se <- se_list[[section]]
  dest <- file.path(out_dir, section)
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)

  counts <- tryCatch(
    GetAssayData(se, assay = "MSI", layer = "counts"),
    error = function(e) GetAssayData(se, assay = "MSI", slot = "counts")
  )
  # Seurat appends "_1" when merging; strip it to get the 10x barcodes back
  barcodes <- sub("_[0-9]+$", "", colnames(counts))

  writeMM(as(counts, "CsparseMatrix"), file.path(dest, "matrix.mtx"))  # features x spots
  writeLines(barcodes, file.path(dest, "barcodes.tsv"))
  writeLines(rownames(counts), file.path(dest, "features.tsv"))

  meta <- se@meta.data
  meta$barcode <- barcodes
  write.csv(meta, file.path(dest, "metadata.csv"), row.names = FALSE)

  cat(section, ":", ncol(counts), "spots,", nrow(counts), "m/z features\n")
}
