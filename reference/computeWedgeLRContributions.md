# Decompose a Wedge's two hub paths into per-edge ligand-receptor contributions

A Wedge (see
[`calculateWedges()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/calculateWedges.md))
sums two directed 2-hop paths through hub `j`: `i -> j -> k` and
`k -> j -> i`. Since Wedges multiply their two edges together rather
than sum them, there is no single well-defined "share of the Wedge
score" attributable to one edge's ligand-receptor pair alone – instead,
this decomposes each of the 4 edges across both paths separately via
[`computeLRContributions()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/computeLRContributions.md)
(`bundle = FALSE`, one direction per edge), so `PctOfTotal` in the
result is relative to each edge's own kernel value.

## Usage

``` r
computeWedgeLRContributions(
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

- i, k:

  Cell-type names for the wedge's two leaves.

- j:

  Cell-type name for the wedge's hub.

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
[`computeLRContributions()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/computeLRContributions.md)'s
(pooled across the 4 edges), with an added `Path` column
(`"i -> j -> k"` or `"k -> j -> i"`, using the actual cell-type names).
