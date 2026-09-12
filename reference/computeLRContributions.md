# Decompose a Direct communication score into per-ligand-receptor-pair contributions

[`compute_kernel()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/compute_kernel.md)
aggregates every ligand-receptor pair into
`kernel[sender, receiver] = lig_weight[sender, ] %*% LRmatrix %*% rec_weight[receiver, ]`
– a sum over all (ligand, receptor) pairs – so once the kernel is built,
individual ligand-receptor identity is gone from the number alone.
Because that formula is linear in `LRmatrix`, the contribution of one
specific ligand-receptor pair `(lig, rec)` to `kernel[sender, receiver]`
is just the single term
`lig_weight[sender, lig] * LRmatrix[lig, rec] * rec_weight[receiver, rec]`,
and these terms sum exactly back to the full kernel value. This function
computes every such term for one sender-receiver pair, so the (weighted)
importance of each ligand-receptor pair to that pair's communication
score can be ranked.

## Usage

``` r
computeLRContributions(
  Lmatrix,
  Rmatrix,
  Cmatrix,
  LRmatrix,
  celltypes,
  sender,
  receiver,
  Dcell = Cmatrix,
  bundle = TRUE,
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

- sender, receiver:

  Cell-type names (must be in `celltypes`) for the pair of interest.

- Dcell:

  Patient-by-cell-type abundance matrix used to scale raw contributions
  to the same units as
  [`calculateDirect()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/calculateDirect.md)'s
  raw score. Defaults to `Cmatrix`.

- bundle:

  Logical; if `TRUE` (default), pools both directions
  (`sender -> receiver` and `receiver -> sender`) under one `Edge` label
  and one shared `PctOfTotal` denominator, matching
  [`calculateDirect()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/calculateDirect.md)'s
  bundled score exactly (a self-pair, `sender == receiver`, is doubled
  the same way
  [`calculateDirect()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/calculateDirect.md)
  doubles it). If `FALSE`, only the `sender -> receiver` direction is
  returned.

- patient_names:

  Optional character vector of patient labels, matching the patient
  order in `LRmatrix`. Defaults to `Patient_1`, `Patient_2`, ...

- top_n:

  Optional; if given, keep only the top `top_n` ligand-receptor pairs
  per patient (by `Contribution`).

## Value

A data frame with one row per (patient, ligand-receptor pair):
`Patient`, `Edge` (a label shared by rows that pool into one
`PctOfTotal` denominator), `Sender`, `Receiver`, `Ligand`, `Receptor`,
`Contribution` (in the same units as
[`calculateDirect()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/calculateDirect.md)'s
raw score), and `PctOfTotal` (share of `Contribution` within that
patient and `Edge`). Contributions that are exactly zero (no possible
pathway) are omitted. `NULL` (with a message) if no ligand-receptor pair
has a nonzero contribution at all.
