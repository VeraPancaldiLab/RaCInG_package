# Run the full kernel-based RaCInG workflow

Run the full kernel-based RaCInG workflow

## Usage

``` r
compute_racing_kernel(
  counts = NULL,
  output_folder = "~/Documents/racing/vignettes/",
  deconv = NULL,
  cc_network = NULL,
  fun_LR = min,
  cell_expr_profile = NULL,
  source = "source_genesymbol",
  target = "target_genesymbol",
  signed = FALSE,
  deconv_method = "Quantiseq",
  cbsx.name = NULL,
  cbsx.token = NULL,
  file_name = NULL,
  nPatients = "all",
  communication_type = "W",
  norm = TRUE,
  pt_idx = NULL,
  remove_direction = TRUE,
  input_data = NULL
)
```

## Arguments

- counts:

  Gene-by-sample count matrix. Required when `input_data` is not
  supplied; ignored otherwise.

- output_folder:

  Directory used to write and read intermediate input files.

- deconv:

  Optional patient-by-cell-type abundance matrix. If omitted, it is
  computed via
  [`multideconv::compute.deconvolution()`](https://rdrr.io/pkg/multideconv/man/compute.deconvolution.html)
  followed by
  [`multideconv::compute.deconvolution.analysis()`](https://rdrr.io/pkg/multideconv/man/compute.deconvolution.analysis.html)
  (which identifies and collapses correlated cell-type subgroups) and
  [`multideconv::standardize_celltype_colnames()`](https://rdrr.io/pkg/multideconv/man/standardize_celltype_colnames.html).
  See
  [`prepare_input_files()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/prepare_input_files.md).

- cc_network:

  Optional ligand-receptor prior network. If omitted, it is retrieved
  via
  [`liana::get_curated_omni()`](https://saezlab.github.io/liana/reference/get_curated_omni.html).
  See
  [`prepare_input_files()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/prepare_input_files.md).

- fun_LR:

  Function used to combine ligand and receptor expression values.

- cell_expr_profile:

  Optional gene-by-cell-type expression profile matrix. If omitted, it
  is estimated from `counts` and `deconv` via per-gene non-negative
  least squares. See
  [`prepare_input_files()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/prepare_input_files.md).

- source, target:

  Column names to use as ligand and receptor identifiers in
  `cc_network`.

- signed:

  Logical; if `TRUE`, also try to load a sign matrix.

- deconv_method:

  Deconvolution method(s) used when `deconv` is not supplied.

- cbsx.name, cbsx.token:

  Optional credentials for the deconvolution workflow.

- file_name:

  File stem used for intermediate input files.

- nPatients:

  Number of patients to process, or `"all"`.

- communication_type:

  Feature family to compute.

- norm:

  Logical; if `TRUE`, compute a normalized baseline kernel.

- pt_idx:

  Optional single patient index to process.

- remove_direction:

  Logical; if `TRUE`, merge directionally equivalent features.

- input_data:

  Optional named list of pre-computed input matrices as returned by
  [`prepare_input_files()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/prepare_input_files.md).
  Must contain `Lmatrix`, `Rmatrix`, `Cmatrix`, `LRmatrix`, `celltypes`,
  `ligands`, and `receptors`. When supplied, the `counts` argument and
  all preprocessing parameters (`deconv`, `cc_network`, etc.) are
  ignored.

## Value

A list with the kernel arrays and the derived feature matrix.
