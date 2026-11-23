# Extract per-spot tables from the authors' Seurat objects shipped inside the
# Kuppe 2022 Space Ranger archives (Zenodo 6580069). Descriptive only: copies
# what ships, computes nothing new.
# Usage: Rscript 01_extract_seurat_tables.R <visium_raw dir> <out dir> <tmp dir>
suppressPackageStartupMessages({ library(SeuratObject); library(Seurat) })
args <- commandArgs(trailingOnly = TRUE)
raw_dir <- args[1]; out_dir <- args[2]; tmp_dir <- args[3]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
for (tarf in sort(Sys.glob(file.path(raw_dir, "10X_Visium_*.tar.gz")))) {
  id <- sub("10X_Visium_(.*)\\.tar\\.gz", "\\1", basename(tarf))
  out <- file.path(out_dir, paste0(id, "_spots.csv"))
  if (file.exists(out)) { cat("skip", id, "\n"); next }
  member <- paste0("./", id, "/", id, ".rds")
  untar(tarf, files = member, exdir = tmp_dir)
  rds <- file.path(tmp_dir, id, paste0(id, ".rds"))
  x <- readRDS(rds)
  md <- x@meta.data
  img <- x@images[[1]]
  co <- img@coordinates  # tissue, row, col, imagerow, imagecol (full-resolution pixels)
  co <- co[rownames(md), c("row", "col", "imagerow", "imagecol")]
  colnames(co) <- c("array_row", "array_col", "pxl_row_fullres", "pxl_col_fullres")
  props <- t(as.matrix(GetAssayData(x, assay = "c2l_props", layer = "data")))
  colnames(props) <- paste0("prop_", colnames(props))
  abund <- t(as.matrix(GetAssayData(x, assay = "c2l", layer = "data")))
  tab <- data.frame(barcode = rownames(md), md, co, props[rownames(md), , drop = FALSE],
                    c2l_total_abundance = rowSums(abund[rownames(md), , drop = FALSE]),
                    check.names = FALSE)
  tab$hca_sample_id <- id
  write.csv(tab, out, row.names = FALSE)
  sf <- img@scale.factors
  write.csv(data.frame(hca_sample_id = id, spot_diameter_fullres = sf$spot,
                       tissue_hires_scalef = sf$hires, tissue_lowres_scalef = sf$lowres,
                       assays = paste(Assays(x), collapse = ";")),
            file.path(out_dir, paste0(id, "_meta.csv")), row.names = FALSE)
  unlink(file.path(tmp_dir, id), recursive = TRUE)
  cat("done", id, nrow(tab), "spots\n")
}
