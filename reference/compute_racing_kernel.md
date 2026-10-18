# Run the full kernel-based RaCInG workflow

Run the full kernel-based RaCInG workflow

## Usage

``` r
compute_racing_kernel(
  counts = NULL,
  output_folder = "Results/",
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
  communication_type = "D",
  norm = FALSE,
  pt_idx = NULL,
  remove_direction = TRUE,
  input_data = NULL,
  zero_threshold = 1,
  expr_threshold = 10,
  already_normalized = FALSE
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
  [`multideconv::compute.deconvolution()`](https://verapancaldilab.github.io/multideconv/reference/compute.deconvolution.html)
  followed by
  [`multideconv::compute.deconvolution.analysis()`](https://verapancaldilab.github.io/multideconv/reference/compute.deconvolution.analysis.html)
  (which identifies and collapses correlated cell-type subgroups;
  requires multideconv \>= 0.2.0). Columns keep multideconv's names:
  subgroups (`"B.cells_Subgroup.1"`) and unclustered features
  (`"Quantiseq_TIL10_B.cells"`). See
  [`prepare_input_files()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/prepare_input_files.md).

- cc_network:

  Optional ligand-receptor prior network. If omitted, it is retrieved
  via
  [`liana::get_curated_omni()`](https://saezlab.github.io/liana/reference/get_curated_omni.html).
  See
  [`prepare_input_files()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/prepare_input_files.md).

- fun_LR:

  Function combining a ligand's and receptor's expression values (as a
  length-2 vector) into one interaction strength. Default `min`: the
  weaker-expressed partner limits the interaction. Any function
  collapsing 2 values to 1 works, e.g. `prod` for the product of the two
  instead.

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

  Feature family to compute. The kernel method only computes direct
  communication (`"D"`, the default); GSCC is available with
  [`compute_racing_montecarlo()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/compute_racing_montecarlo.md).

- norm:

  Logical; if `TRUE`, also compute a normalized baseline kernel and
  express features as an enrichment ratio over it (isolates specificity
  from abundance/topology). Default `FALSE` returns the raw,
  abundance-weighted communication magnitude.

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

- zero_threshold:

  Drop a feature once its fraction of zero-valued patients reaches this
  threshold. Default `1` only drops features that are zero for every
  patient (the original behavior); lower it (e.g. `0.9`) to also drop
  merely zero-inflated features. See
  [`compute_kernel_features()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/compute_kernel_features.md).

- expr_threshold:

  Minimum expression (in linear TPM, not log-scale) a ligand or receptor
  must reach in its sending/receiving cell type to be kept. See
  [`prepare_input_files()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/prepare_input_files.md).

- already_normalized:

  Logical; if `TRUE`, `counts` is taken as linear-scale TPM and not
  renormalized. See
  [`prepare_input_files()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/prepare_input_files.md).

## Value

A list with three elements:

- `kernel`:

  The raw kernel array `[sender x receiver x patient]` from
  [`compute_kernel()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/compute_kernel.md).

- `kernel_norm`:

  The uniformized baseline kernel when `norm = TRUE`, otherwise `NULL`.

- `features`:

  A patient-by-feature data frame of direct communication scores, with
  columns named `"Dir_A_B"` and rows named by patient.

## Examples

``` r
data(skcm_example)
res <- compute_racing_kernel(input_data = skcm_example)
#> Using pre-computed input matrices; skipping input generation.
#> Computing kernel for 10 patients
#> Calculating kernel...
#> Calculating features...
dim(res$features)
#> [1] 10 45
head(res$features[, 1:4])
#>                Dir_B_B   Dir_B_CAF Dir_B_CD8+ T    Dir_B_DC
#> Patient_1 0.0018678610 0.007140434 0.0004890447 0.001866252
#> Patient_2 0.0012929868 0.003427966 0.0107989029 0.007941141
#> Patient_3 0.0151493949 0.009916129 0.0217169419 0.013045575
#> Patient_4 0.0000000000 0.000000000 0.0000000000 0.000000000
#> Patient_5 0.0012159134 0.005778287 0.0000000000 0.002329358
#> Patient_6 0.0004382255 0.002512236 0.0056910729 0.002374521

# Enrichment ratio over a uniform ligand-receptor baseline
res_norm <- compute_racing_kernel(input_data = skcm_example, norm = TRUE)
#> Using pre-computed input matrices; skipping input generation.
#> Computing kernel for 10 patients
#> Calculating kernel...
#> Calculating features...
head(res_norm$features[, 1:4])
#>             Dir_B_B Dir_B_CAF Dir_B_CD8+ T   Dir_B_DC
#> Patient_1 0.3655633 0.6472004    0.2193690 0.33179690
#> Patient_2 0.7127917 0.5965352    1.0316271 0.84989477
#> Patient_3 1.1609442 0.7959182    0.8781016 0.92070117
#> Patient_4 0.4675364 0.4900159    0.1068774 0.09173062
#> Patient_5 0.2765913 0.5799559    0.1661719 0.37846202
#> Patient_6 0.4209898 0.5559673    0.7787024 0.51126999
```
