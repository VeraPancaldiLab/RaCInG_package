# Compute GSCC features from kernel matrices

Compute GSCC features from kernel matrices

## Usage

``` r
computeGSCC(
  kernel,
  Dcell,
  cell_names,
  patient_names,
  unifKernel = NULL,
  norm = FALSE,
  lab = 1
)
```

## Arguments

- kernel:

  Kernel array from
  [`compute_kernel()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/compute_kernel.md).

- Dcell:

  Patient-by-cell-type abundance matrix.

- cell_names:

  Character vector of cell-type labels.

- patient_names:

  Character vector of patient names.

- unifKernel:

  Optional normalized baseline kernel.

- norm:

  Logical; if `TRUE`, divide by the baseline GSCC values.

- lab:

  Scaling factor for interaction strengths.

## Value

A data frame with GSCC feature values per patient.
