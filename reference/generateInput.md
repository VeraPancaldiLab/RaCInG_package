# Load RaCInG input matrices from disk –\> Is there any normalization done here? I dont see it in the function

Internal helper called by
[`prepare_input_files()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/prepare_input_files.md)
to read back the generated CSV files and assemble the normalised
matrices and 3-D tensor. Users should call
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
