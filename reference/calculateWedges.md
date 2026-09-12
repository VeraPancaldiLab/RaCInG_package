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
  patient_names = NULL,
  zero_threshold = 1
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

- zero_threshold:

  Drop a wedge once its fraction of zero-valued patients reaches this
  threshold. Default `1` only drops wedges that are zero for every
  patient (the original behavior); lower it (e.g. `0.9`) to also drop
  merely zero-inflated wedges.

## Value

A patient-by-feature data frame of wedge scores. Triplets with no
possible ligand-receptor pathway in any patient (zero for every patient)
are dropped rather than returned as all-zero columns.
