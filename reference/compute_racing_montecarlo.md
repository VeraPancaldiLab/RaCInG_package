# Run the full Monte Carlo RaCInG workflow

Run the full Monte Carlo RaCInG workflow

## Usage

``` r
compute_racing_montecarlo(
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
  pt_idx = NULL,
  file_name = NULL,
  nPatients = "all",
  communication_type = "W",
  Ncells = 10000,
  Ngraphs = 100,
  Ndegree = 20,
  remove_direction = TRUE,
  norm = TRUE,
  input_data = NULL,
  ncores = 1
)
```

## Arguments

- counts:

  Gene-by-sample count matrix. Required when `input_data` is not
  supplied; ignored otherwise.

- output_folder:

  Directory used to write intermediate and output files.

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

- pt_idx:

  Optional single patient index to simulate.

- file_name:

  File stem used for intermediate files.

- nPatients:

  Number of patients to process, or `"all"`.

- communication_type:

  Feature family to simulate.

- Ncells:

  Number of cells per simulated graph.

- Ngraphs:

  Number of Monte Carlo iterations.

- Ndegree:

  Target average degree.

- remove_direction:

  Logical; if `TRUE`, merge directionally equivalent features.

- norm:

  Logical; if `TRUE`, also run a uniformized baseline simulation for
  normalization.

- input_data:

  Optional named list of pre-computed input matrices as returned by
  [`prepare_input_files()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/prepare_input_files.md).
  Must contain `Lmatrix`, `Rmatrix`, `Cmatrix`, `LRmatrix`, `celltypes`,
  `ligands`, and `receptors`. When supplied, the `counts` argument and
  all preprocessing parameters (`deconv`, `cc_network`, etc.) are
  ignored.

- ncores:

  Number of cores to compute patients on in parallel. Passed through to
  [`runSim()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/runSim.md);
  see its documentation for details.

## Value

A list with the generated inputs and processed feature matrices.
