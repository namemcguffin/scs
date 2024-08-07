suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(purrr)
})

cl_args <- commandArgs(trailingOnly = TRUE)

so <- readRDS(cl_args[[1]])
out_path <- cl_args[[2]]
spatial_data <- GetTissueCoordinates(so, if (length(cl_args) >= 3) {
  cl_args[[3]]
} else {
  NULL
})[c("cell", "x", "y")]
mtx <- GetAssayData(
  so,
  assay = if (length(cl_args) >= 4) {
    cl_args[[4]]
  } else {
    NULL
  },
  layer = if (length(cl_args) >= 5) {
    cl_args[[5]]
  } else {
    NULL
  }
) |>
  t() |>
  as.data.table(keep.rownames = "cell")

dir.create(out_path)
dir.create(file.path(out_path, "feat"))

fwrite(
  spatial_data,
  file.path(out_path, "cells.tsv"),
  sep = "\t"
)

walk(
  Features(so),
  function(feat) {
    fwrite(
      mtx[get(feat) > 0, .(cell, expr = get(feat))],
      file.path(out_path, "feat", sprintf("%s.tsv", feat)),
      sep = "\t"
    )
  }
)
