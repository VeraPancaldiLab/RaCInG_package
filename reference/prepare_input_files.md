# Build RaCInG input files from raw count data

This function combines the preprocessing and input-loading steps into a
single call. It generates the `L`, `R`, `C`, and `LR` CSV files from raw
counts, then reads them back to produce the normalised matrices and 3-D
tensor required by the kernel and Monte Carlo workflows.

## Usage

``` r
prepare_input_files(
  counts,
  output_folder = "Results/",
  deconv = NULL,
  cc_network = NULL,
  fun_LR = min,
  cell_expr_profile = NULL,
  source = "source_genesymbol",
  target = "target_genesymbol",
  deconv_method = "Quantiseq",
  cbsx.name = NULL,
  cbsx.token = NULL,
  file_name = NULL,
  signed = FALSE,
  already_normalized = FALSE,
  expr_threshold = 10
)
```

## Arguments

- counts:

  Gene-by-sample count matrix.

- output_folder:

  Directory where the generated `L`, `R`, `C`, and `LR` files are
  written.

- deconv:

  Optional patient-by-cell-type abundance matrix. If omitted, it is
  computed via
  [`multideconv::compute.deconvolution()`](https://rdrr.io/pkg/multideconv/man/compute.deconvolution.html)
  followed by
  [`multideconv::compute.deconvolution.analysis()`](https://rdrr.io/pkg/multideconv/man/compute.deconvolution.analysis.html)
  (which identifies and collapses correlated cell-type subgroups) and
  [`multideconv::standardize_celltype_colnames()`](https://rdrr.io/pkg/multideconv/man/standardize_celltype_colnames.html).

- cc_network:

  Optional ligand-receptor prior network. If omitted, it is retrieved
  via
  [`liana::get_curated_omni()`](https://saezlab.github.io/liana/reference/get_curated_omni.html)
  (a curated, OmniPath-backed consensus of
  CellPhoneDB/CellChatDB/ICELLNET/connectomeDB2020/CellTalkDB) the first
  time this is called, then cached indefinitely as an RDS file under
  `tools::R_user_dir("RaCInG", "cache")` and loaded from there on every
  subsequent call – avoiding a live network round-trip (and its
  associated transient-failure risk) on every invocation. Delete that
  cache file to force a refresh from OmniPath.

- fun_LR:

  Function combining a ligand's and receptor's expression values (as a
  length-2 vector) into one interaction strength. Default `min`: the
  weaker-expressed partner limits the interaction. Any function
  collapsing 2 values to 1 works, e.g. `prod` for the product of the two
  instead.

- cell_expr_profile:

  Optional gene-by-cell-type expression profile matrix. If omitted, it
  is estimated from `counts` and `deconv` via per-gene non-negative
  least squares (see
  [`.estimate_expression_profiles()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/dot-estimate_expression_profiles.md)).

- source, target:

  Column names to use as ligand and receptor identifiers when
  `cc_network` is supplied.

- deconv_method:

  Deconvolution method(s) passed to
  [`multideconv::compute.deconvolution()`](https://rdrr.io/pkg/multideconv/man/compute.deconvolution.html).

- cbsx.name, cbsx.token:

  Optional credentials forwarded to the deconvolution workflow.

- file_name:

  File stem used when exporting the generated CSV files.

- signed:

  Logical; if `TRUE`, also try to load a sign matrix from
  `output_folder`.

- already_normalized:

  Logical; if `TRUE`, `counts` is treated as already TPM-normalized,
  **linear scale, not logged**, and the internal
  [`ADImpute::NormalizeTPM()`](https://rdrr.io/pkg/ADImpute/man/NormalizeTPM.html)
  gene-length/library-size normalization is skipped. Use this for
  datasets where only pre-normalized expression is available (e.g.
  public cohorts distributed as log2(TPM+1) matrices with no raw counts
  – convert those back to linear scale, e.g. `2^x - 1`, before passing
  them here) – applying `NormalizeTPM()` to already- normalized data
  would re-apply its gene-length correction on top of the original
  normalization, distorting relative gene expression levels. Emits a
  warning, since `expr_threshold` below is only meaningful if this
  linear-scale assumption holds – comparing a log-scale value against it
  would silently apply the wrong cutoff.

- expr_threshold:

  Minimum expression value, in linear-scale TPM (in `cell_expr_profile`
  units), a ligand or receptor must reach in its sending/receiving cell
  type for that ligand-receptor-celltype triple to be kept in
  `CC_table`.

## Value

A named list with the processed input matrices and their labels:
`Lmatrix`, `Rmatrix`, `Cmatrix` (normalised), `LRmatrix` (3-D tensor),
`celltypes`, `ligands`, `receptors`, `Sign_matrix`, and `CC_table`.
