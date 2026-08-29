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
  signed = FALSE
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
  CellPhoneDB/CellChatDB/ICELLNET/connectomeDB2020/CellTalkDB).

- fun_LR:

  Function used to combine ligand and receptor expression values.

- cell_expr_profile:

  Optional gene-by-cell-type expression profile matrix. If omitted, it
  is estimated from `counts` and `deconv` via per-gene non-negative
  least squares (see
  [`.estimate_expression_profiles()`](https://mhurtado13.github.io/racing/reference/dot-estimate_expression_profiles.md)).

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

## Value

A named list with the processed input matrices and their labels:
`Lmatrix`, `Rmatrix`, `Cmatrix` (normalised), `LRmatrix` (3-D tensor),
`celltypes`, `ligands`, `receptors`, `Sign_matrix`, and `CC_table`.
