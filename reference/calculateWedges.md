# Calculate wedge features from a kernel

Calculate wedge features from a kernel

## Usage

``` r
calculateWedges(
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
  (unnormalized) scores. Ignored when `unifKernel` is supplied.

- bundle:

  Logical; if `TRUE`, combine directionally equivalent wedges.

- patient_names:

  Optional character vector of patient labels, matching the patient
  order in `kernel`. Defaults to `Patient_1`, `Patient_2`, ... when
  omitted.

## Value

A patient-by-feature data frame of wedge scores. Triplets with no
possible ligand-receptor pathway in any patient (zero for every patient)
are dropped rather than returned as all-zero columns.
