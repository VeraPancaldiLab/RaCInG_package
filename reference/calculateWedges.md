# Calculate wedge features from a kernel

Calculate wedge features from a kernel

## Usage

``` r
calculateWedges(kernel, unifKernel = NULL, cells, Dcell = NULL, bundle = TRUE)
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
  (unnormalized) scores. Ignored when `unifKernel` is supplied.

- bundle:

  Logical; if `TRUE`, combine directionally equivalent wedges.

## Value

A patient-by-feature data frame of wedge scores.
