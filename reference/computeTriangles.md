# Compute triangle features from kernel matrices

Compute triangle features from kernel matrices

## Usage

``` r
computeTriangles(
  kernel,
  cell_names,
  patient_names,
  Dcell = NULL,
  unifKernel = NULL,
  norm = FALSE,
  bundle = TRUE
)
```

## Arguments

- kernel:

  Kernel array from
  [`compute_kernel()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/compute_kernel.md).

- cell_names:

  Character vector of cell-type names.

- patient_names:

  Character vector of patient names.

- Dcell:

  Patient-by-cell-type abundance matrix, used to weight raw
  (unnormalized) scores. Ignored when `unifKernel`/`norm` is used.

- unifKernel:

  Optional normalized baseline kernel.

- norm:

  Logical; if `TRUE`, divide by the baseline triangle scores.

- bundle:

  Logical; if `TRUE`, aggregate all directions into a single "Tr"
  triangle score per triple. If `FALSE`, two scores are returned per
  triple: "TT" (trust/transitive triangle, `i->j->k` and `i->k`) and,
  for `i<=j<=k` only, "CT" (cycle triangle, `i->j->k->i`).

## Value

A patient-by-feature data frame of triangle scores.
