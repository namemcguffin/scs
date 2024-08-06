# scs: spatial cell selector

## description

`scs` is a program that allows for the selection of spatially annotated points within a dataset using arbitrary bounding boxes.

given a dataset as input, once selection is complete, it will generate as output a TSV file with two columns, `cell` and `within`,
which specify the relevant cell identifier and if that cell is within the defined bounding boxe(s).

## usage

requires two arguments: the input directory path and an output path.

note: it is heavily recommended to compile `scs` on release mode to achieve consistent performance.

example:
```bash
cargo run --release 'input/public' 'output.tsv'
```

## input directory structure

`scs` expects the input directory to have the following structure:
- a top-level `cells.tsv` file with three columns: `cell` (specifying a cell identifier), `x` (specifying the cell's x coordinate), and `y` (specifying the cell's y coordinate)
- a directory called `feat`, which contains feature files
  - each feature file should have the name `{FEATURE_NAME}.tsv` (e.g. `Yap1.tsv` for the `Yap1` feature)
  - each feature file should have two columns: `cell` (specifying a cell identifier) and `expr` (specifying the associated value of that feature for that cell)
  - any cell not present in a feature file is assumed to have a value of `0` for that feature (given feature data is often very sparse this can keep feature file size small)
- all input TSV files should contain a header line

the `scripts` directory contains scripts that can generate input directories for `scs` using common existing data formats (such as seurat)

### generating an input directory from seurat objects: `scripts/seurat.r`

generates an input directory from a seurat object saved as an RDS file.

requires the `Seurat`, `data.table`, `purrr` packages to be installed.

arguments:
1. input RDS file path
1. desired output directory
1. centroid FOV to use for x/y coordinates (optional)
1. assay to use for features (optional)
1. layer to use from selected assay for features (optional)

example:
```bash
Rscript 'scripts/seurat.r' 'input/so.rds' 'input/data' # uses seurat's default values for FOV, assay, and layer
Rscript 'scripts/seurat.r' 'input/so.rds' 'input/data' 'centroids.specific' 'vizgen' 'data' # specific example values
```
