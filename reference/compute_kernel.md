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

  Logical; if `TRUE`, also compute a uniformized baseline kernel.

## Value

Either a 3D kernel array `[sender x receiver x patient]`, or (when
`normalize = TRUE`) a list with `kernel` and `kernel_norm`.
