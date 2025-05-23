# scs: spatial cell selector

## description

`scs` is a program that allows for the selection of spatially annotated points within a dataset using arbitrary bounding boxes.

given a dataset as input, once selection is complete, it will generate as output a TSV file with two columns, `cell` and `within_bb`,
which specify the relevant cell identifier and if that cell is within the defined bounding boxe(s).

to compile this project, you will need to install a rust compiler, instructions for which can be found [here](https://www.rust-lang.org/tools/install)

## usage

pressing `h` will bring up a small help menu explaining controls.

takes two arguments: the input directory path (mandatory) and an output path (optional, will default to standard out if no path is provided)

note: it is heavily recommended to compile `scs` on the release profile to achieve consistent performance.

example:

```bash
cargo run --release 'input/public' 'output.tsv'
```

### workflow demonstration

a jupyter notebook (`scripts/demo.ipynb`) is provided to demonstrate how to use `scs` within the context of a scanpy-driven analysis.

the notebook requires the `pandas`, `anndata`, `numpy`, and `scanpy` packages to be installed, which can be done via the following command using the uv tool (installation instructions [here](https://docs.astral.sh/uv/getting-started/installation/)):

```bash
uv pip install pandas anndata numpy scanpy
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

### generating an input directory from seurat objects: `scripts/from_seurat.r`

generates an input directory from a seurat object saved as an RDS file.

requires the `Seurat`, `data.table`, `purrr` packages to be installed, which can be done via the following command:

```bash
R -q -e 'if (system.file(package = "pak") == "") { install.packages("pak") }; pak::pak(c("Seurat", "data.table", "purrr"))'
```

arguments:

1. input RDS file path
1. desired output directory
1. centroid FOV to use for x/y coordinates (optional, defers to seurat's defaults)
1. assay to use for features (optional, defers to seurat's defaults)
1. layer to use from selected assay for features (optional, defers to seurat's defaults)

example:

```bash
Rscript 'scripts/seurat.r' 'input/so.rds' 'input/data' # uses seurat's default values for FOV, assay, and layer
Rscript 'scripts/seurat.r' 'input/so.rds' 'input/data' 'centroids.specific' 'vizgen' 'data' # specific values
```

### generating an input directory from anndata objects: `scripts/from_anndata.py`

generates an input directory from an anndata object saved as an h5ad file.

note: expects spatial coordinates to be in an `obsm` slot, not as a set of `obs` columns.

requires the `anndata` and `pandas` packages to be installed, which can be done via the following command using the uv tool (installation instructions [here](https://docs.astral.sh/uv/getting-started/installation/)):

```bash
uv pip install anndata scanpy
```

arguments:

1. input h5ad file path
1. desired output directory
1. `obsm` index to use for spatial coordinates (optional, defaults to `spatial`)
1. layer index to use (optional, defaults to using `X` matrix)

example (make sure to activate the virtual environment where `anndata` and `scanpy` were installed before running):

```bash
python 'scripts/from_anndata.py' 'input/adata.h5ad' 'input/data' # uses default values
python 'scripts/from_anndata.py' 'input/adata.h5ad' 'input/data' 'spatial_centroids' 'counts' # specific values
```
