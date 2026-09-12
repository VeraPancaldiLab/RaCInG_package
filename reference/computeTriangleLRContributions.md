# Decompose a Triangle's edges into per-edge ligand-receptor contributions

A bundled Triangle (see
[`computeTriangles()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/computeTriangles.md))
combines all 6 directed edges among 3 cell types (`i`, `j`, `k`) into 8
direction-combination products. As with Wedges, edges are multiplied
together rather than summed, so there is no single well-defined "share
of the Triangle score" attributable to one edge's ligand-receptor pair
alone – instead, this decomposes each of the 6 directed edges separately
via
[`computeLRContributions()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/computeLRContributions.md)
(`bundle = FALSE`), so `PctOfTotal` in the result is relative to each
edge's own kernel value.

## Usage

``` r
computeTriangleLRContributions(
  Lmatrix,
  Rmatrix,
  Cmatrix,
  LRmatrix,
  celltypes,
  i,
  j,
  k,
  Dcell = Cmatrix,
  patient_names = NULL,
  top_n = NULL
)
```

## Arguments

- Lmatrix, Rmatrix, Cmatrix, LRmatrix:

  Same inputs as
  [`compute_kernel()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/compute_kernel.md).

- celltypes:

  Character vector of cell-type names, in the row order of
  `Lmatrix`/`Rmatrix`/columns of `Cmatrix`.

- i, j, k:

  Cell-type names for the triangle's 3 vertices.

- Dcell:

  Patient-by-cell-type abundance matrix used to scale raw contributions
  to the same units as
  [`calculateDirect()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/calculateDirect.md)'s
  raw score. Defaults to `Cmatrix`.

- patient_names:

  Optional character vector of patient labels, matching the patient
  order in `LRmatrix`. Defaults to `Patient_1`, `Patient_2`, ...

- top_n:

  Optional; if given, keep only the top `top_n` ligand-receptor pairs
  per patient (by `Contribution`).

## Value

A data frame like
[`computeLRContributions()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/computeLRContributions.md)'s,
pooled across the 6 directed edges (`Edge` identifies which one, e.g.
`"i -> j"`).
