suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(purrr)
})

cl_args <- commandArgs(trailingOnly = TRUE)
so <- readRDS(cl_args[[1]])
out_path <- cl_args[[2]]

selected_fov <- if (length(cl_args) >= 3) {
  cl_args[[3]]
} else {
  NULL
}
selected_assay <- if (length(cl_args) >= 4) {
  cl_args[[4]]
} else {
  NULL
}
selected_layer <- if (length(cl_args) >= 5) {
  cl_args[[5]]
} else {
  NULL
}

dir.create(out_path)
dir.create(sprintf("%s/feat", out_path))

fwrite(
  GetTissueCoordinates(so, selected_fov)[c("cell", "x", "y")],
  sprintf("%s/cells.tsv", out_path),
  sep = "\t"
)

mtx <- GetAssayData(so, assay = selected_assay, layer = selected_layer)
walk(
  Features(so),
  function(feat) {
    mtx[feat, ] |>
      as.data.table(keep.rownames = T) |>
      `colnames<-`(c("cell", "expr")) |>
      _[expr > 0] |>
      fwrite(sprintf("%s/feat/%s.tsv", out_path, feat), sep = "\t")
  }
)
