# Run the full Monte Carlo RaCInG workflow

Run the full Monte Carlo RaCInG workflow

## Usage

``` r
compute_racing_montecarlo(
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
  pt_idx = NULL,
  file_name = NULL,
  nPatients = "all",
  communication_type = "D",
  Ncells = 10000,
  Ngraphs = 100,
  Ndegree = 20,
  remove_direction = TRUE,
  norm = FALSE,
  input_data = NULL,
  ncores = 1,
  expr_threshold = 10,
  already_normalized = FALSE
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

  Feature family to simulate: `"D"` (direct communication, the default),
  `"GSCC"`, or both (`c("D", "GSCC")`). When both are given, every
  requested type is extracted from the same simulated graphs (via
  [`countAllTypes()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/countAllTypes.md))
  in one pass instead of re-simulating a fresh set of graphs per type –
  graph generation, not feature extraction, is the expensive part of a
  Monte Carlo run, so this amortizes that cost across every type
  requested. `output` in the return value is then a named list (one
  entry per type) instead of a single result.

- Ncells:

  Number of cells per simulated graph.

- Ngraphs:

  Number of Monte Carlo iterations.

- Ndegree:

  Target average degree.

- remove_direction:

  Logical; if `TRUE`, merge the two directions of each cell-type pair in
  the direct communication features.

- norm:

  Logical; if `TRUE`, also run a uniformized baseline simulation and
  express features as an enrichment ratio over it (isolates specificity
  from abundance/topology, at the cost of a second, noisier simulation
  pass – budget a larger `Ngraphs` if enabling this). Default `FALSE`
  returns the raw, abundance-weighted communication magnitude, which is
  simpler, cheaper (one pass), and less noisy at a given `Ngraphs`.

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

- expr_threshold:

  Minimum expression (in linear TPM, not log-scale) a ligand or receptor
  must reach in its sending/receiving cell type to be kept. See
  [`prepare_input_files()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/prepare_input_files.md).

- already_normalized:

  Logical; if `TRUE`, `counts` is taken as linear-scale TPM and not
  renormalized. See
  [`prepare_input_files()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/prepare_input_files.md).

## Value

A list with two elements:

- `input`:

  The input matrices used: `LRmatrix`, `Lmatrix`, `Rmatrix`, `Cmatrix`,
  `celltypes`, `ligands`, `receptors` – the same names as
  [`prepare_input_files()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/prepare_input_files.md),
  so it can be passed back as `input_data`.

- `output`:

  For a single `communication_type`, a list with the patient-by-feature
  data frames `mean` and `sd` (see
  [`compute_results_processing()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/compute_results_processing.md));
  for `c("D", "GSCC")`, a named list with one such entry per type.

As a side effect, the raw simulations (`<file_name>.out`, plus
`<file_name>_norm.out` when `norm = TRUE`) and the processed features
(`<file_name>.csv`, `<file_name>_sd.csv`) are written to
`output_folder`; with several types, each file name gets a `_<type>`
suffix after `file_name`.

## Examples

``` r
data(skcm_example)
set.seed(1)
mc <- compute_racing_montecarlo(
  input_data = skcm_example, output_folder = tempdir(),
  nPatients = 2, communication_type = c("D", "GSCC"),
  Ncells = 100, Ngraphs = 3, Ndegree = 3
)
#> Using pre-computed input matrices; skipping input generation.
#> Running Monte Carlo simulation for 2 patients
#> Processing interaction distribution and generating CSV output...
head(mc$output$D$mean[, 1:4])
#>             Dir_B_B Dir_B_CAF Dir_B_CD8+ T Dir_B_DC
#> Patient_1 0.0000000 2.0000000            0        0
#> Patient_2 0.3333333 0.6666667            4        1
mc$output$GSCC$mean
#>                GSCC_B   GSCC_CAF GSCC_CD8+ T GSCC_DC GSCC_Endo     GSCC_M
#> Patient_1 0.003333333 0.01333333        0.00    0.00      0.01 0.02333333
#> Patient_2 0.013333333 0.01666667        0.08    0.07      0.00 0.07666667
#>           GSCC_NK   GSCC_Treg GSCC_Tumor
#> Patient_1       0 0.006666667  0.5233333
#> Patient_2       0 0.016666667  0.2833333
```
