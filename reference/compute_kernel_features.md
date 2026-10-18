# Derive communication features from a kernel

Derive communication features from a kernel

## Usage

``` r
compute_kernel_features(
  kernel,
  unifKernel = NULL,
  celltypes,
  communication_type = "D",
  bundle = TRUE,
  patient_names = NULL,
  Dcell = NULL,
  patient_idx = NULL,
  zero_threshold = 1
)
```

## Arguments

- kernel:

  Kernel array from
  [`compute_kernel()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/compute_kernel.md).

- unifKernel:

  Optional normalized baseline kernel. When supplied, features are the
  enrichment ratio over it (see
  [`calculateDirect()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/calculateDirect.md)).

- celltypes:

  Character vector of cell-type labels.

- communication_type:

  Feature family to compute. The kernel method only computes direct
  communication (`"D"`); GSCC is available with
  [`compute_racing_montecarlo()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/compute_racing_montecarlo.md).

- bundle:

  Logical; if `TRUE`, merge the two directions of each cell-type pair.

- patient_names:

  Optional patient labels.

- Dcell:

  Patient-by-cell-type abundance matrix. Required when `unifKernel` is
  not supplied (raw, unnormalized features).

- patient_idx:

  Optional patient index subset.

- zero_threshold:

  Drop a feature once its fraction of zero-valued patients reaches this
  threshold. Default `1` only drops features that are zero for every
  patient; lower it (e.g. `0.9`) to also drop merely zero-inflated
  features. Passed through to
  [`calculateDirect()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/calculateDirect.md).

## Value

A patient-by-feature data frame of direct communication scores for the
selected patients.
