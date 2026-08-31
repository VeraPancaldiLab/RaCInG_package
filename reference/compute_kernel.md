# Compute the RaCInG kernel for one or more patients

The kernel is the *unweighted* structural kernel: cell-type abundances
are not baked into it.
[`calculateDirect()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/calculateDirect.md),
[`calculateWedges()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/calculateWedges.md),
[`computeTriangles()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/computeTriangles.md),
and
[`computeGSCC()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/computeGSCC.md)
each multiply in the abundance of every node they use, exactly once per
node, when computing raw (unnormalized) features from this kernel.

## Usage

``` r
compute_kernel(liglist, reclist, Cmatrix, LRmatrix, normalize = FALSE)
```

## Arguments

- liglist:

  Cell-by-ligand compatibility matrix.

- reclist:

  Cell-by-receptor compatibility matrix.

- Cmatrix:

  Patient-by-cell-type abundance matrix.

- LRmatrix:

  Ligand-by-receptor-by-patient interaction tensor.

- normalize:

  Logical; if `TRUE`, also compute a second kernel (`kernel_norm`) under
  a NULL DISTRIBUTION where every ligand-receptor pair that is
  structurally possible (nonzero in `LRmatrix`) is given the same
  uniform strength (`1 / number of active pairs for that patient`),
  instead of its real measured expression-derived strength.
  `kernel_norm` therefore reflects only network topology and cell-type
  abundance – "how much communication would you expect between these two
  cell types if every possible ligand-receptor pair were equally active"
  – with no information about which pairs are actually more or less
  expressed. Downstream functions
  ([`calculateDirect()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/calculateDirect.md),
  [`calculateWedges()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/calculateWedges.md),
  etc.) use this as a baseline to compute a composition-independent
  enrichment score (`kernel / kernel_norm`) instead of the raw,
  abundance-weighted score – the abundance weighting cancels out of that
  ratio algebraically, since both `kernel` and `kernel_norm` are built
  from the same per-patient `lig_weight`/`rec_weight` terms. Use
  `normalize = FALSE` (default) when you want absolute communication
  strength (abundance and real LR expression both matter); use
  `normalize = TRUE` when you want a score that isolates
  specificity/enrichment of a cell-type pair's communication relative to
  what topology alone would predict, independent of how common those
  cell types are or how strong LR expression is overall.

## Value

Either a 3D kernel array `[sender x receiver x patient]`, or (when
`normalize = TRUE`) a list with `kernel` and `kernel_norm`.
