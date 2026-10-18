# Compute the RaCInG kernel for one or more patients

The kernel is the *unweighted* structural kernel: cell-type abundances
are not baked into it.
[`calculateDirect()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/calculateDirect.md)
multiplies in the abundance of both cell types of a pair when computing
raw (unnormalized) features from it.

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
  expressed.
  [`calculateDirect()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/calculateDirect.md)
  uses this as a baseline to compute a composition-independent
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

## Examples

``` r
data(skcm_example)
k <- compute_kernel(skcm_example$Lmatrix, skcm_example$Rmatrix,
                    skcm_example$Cmatrix, skcm_example$LRmatrix)
dim(k)
#> [1]  9  9 10
feats <- compute_kernel_features(k, celltypes = skcm_example$celltypes,
                                 Dcell = skcm_example$Cmatrix,
                                 patient_names = rownames(skcm_example$Cmatrix))
head(feats[, 1:4])
#>                Dir_B_B   Dir_B_CAF Dir_B_CD8+ T    Dir_B_DC
#> Patient_1 0.0018678610 0.007140434 0.0004890447 0.001866252
#> Patient_2 0.0012929868 0.003427966 0.0107989029 0.007941141
#> Patient_3 0.0151493949 0.009916129 0.0217169419 0.013045575
#> Patient_4 0.0000000000 0.000000000 0.0000000000 0.000000000
#> Patient_5 0.0012159134 0.005778287 0.0000000000 0.002329358
#> Patient_6 0.0004382255 0.002512236 0.0056910729 0.002374521
```
