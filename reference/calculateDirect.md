# Calculate direct communication features from a kernel

Calculate direct communication features from a kernel

## Usage

``` r
calculateDirect(kernel, unifKernel = NULL, cells, Dcell = NULL, bundle = TRUE)
```

## Arguments

- kernel:

  Kernel array returned by
  [`compute_kernel()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/compute_kernel.md).

- unifKernel:

  Optional normalized baseline kernel.

- cells:

  Character vector of cell-type names.

- Dcell:

  Patient-by-cell-type abundance matrix, used to weight raw
  (unnormalized) scores. Ignored when `unifKernel` is supplied, since
  the abundance weights cancel out of the ratio.

- bundle:

  Logical; if `TRUE`, combine reciprocal directions.

## Value

A patient-by-feature data frame of direct communication scores.
