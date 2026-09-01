# Calculate direct communication features from a kernel

Calculate direct communication features from a kernel

## Usage

``` r
calculateDirect(
  kernel,
  unifKernel = NULL,
  cells,
  Dcell = NULL,
  bundle = TRUE,
  patient_names = NULL
)
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

- patient_names:

  Optional character vector of patient labels, matching the patient
  order in `kernel`. Defaults to `Patient_1`, `Patient_2`, ... when
  omitted, matching
  [`compute_racing_montecarlo()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/compute_racing_montecarlo.md)'s
  output when it is also given real names.

## Value

A patient-by-feature data frame of direct communication scores.
Cell-type pairs with no possible ligand-receptor pathway in any patient
(zero for every patient) are dropped rather than returned as all-zero
columns.
