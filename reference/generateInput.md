# Load RaCInG input matrices from disk

Internal helper called by
[`prepare_input_files()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/prepare_input_files.md)
to read back the generated CSV files and assemble the normalised
matrices and 3-D tensor. Normalisation happens in the readers it calls:
[`createCellTypeDistr()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/createCellTypeDistr.md)
rescales each patient's cell-type abundances to sum to 1, and
[`createInteractionDistr()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/createInteractionDistr.md)
rescales each patient's ligand-receptor weights to sum to 1. Users
should call
[`prepare_input_files()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/prepare_input_files.md)
instead.

## Usage

``` r
generateInput(file_name, output_folder, read_signs = FALSE)
```

## Arguments

- file_name:

  File stem used for the `Lmatrix_`, `Rmatrix_`, `Cmatrix_`, and
  `LRmatrix_` CSV files.

- output_folder:

  Directory containing the exported input files.

- read_signs:

  Logical; if `TRUE`, attempts to read `Sign_matrix_<file_name>.csv`
  from `output_folder`.

## Value

A named list containing the input matrices and their labels.
