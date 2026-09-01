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
  norm = FALSE,
  patient_idx = NULL
)
```

## Arguments

- kernel:

  Kernel array from
  [`compute_kernel()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/compute_kernel.md).

- unifKernel:

  Optional normalized baseline kernel.

- celltypes:

  Character vector of cell-type labels.

- communication_type:

  Feature family to compute (`"D"`, `"W"`, `"TT"`, or `"GSCC"`; `"TT"`
  returns both trust- and cycle-triangle columns together, see
  [`computeTriangles()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/computeTriangles.md)),
  or a character vector of several of these. Since the kernel itself
  (`kernel`/`unifKernel`) is passed in already computed, requesting
  multiple types here costs nothing extra to compute per type beyond
  that one kernel – no re-derivation happens for any of them. With more
  than one type, the return value is a named list (one data frame per
  type) instead of a single data frame.

- bundle:

  Logical; if `TRUE`, merge directionally equivalent features where
  appropriate.

- patient_names:

  Optional patient labels.

- Dcell:

  Patient-by-cell-type abundance matrix. Always required for `"GSCC"`;
  required for `"D"`, `"W"`, `"TT"` only when `unifKernel` is not
  supplied (raw, unnormalized features).

- norm:

  Logical; if `TRUE`, compute normalized features when a baseline is
  supplied.

- patient_idx:

  Optional patient index subset.

## Value

A data frame of feature values for the selected patients, or (when
`communication_type` has more than one entry) a named list of such data
frames, one per requested type.
