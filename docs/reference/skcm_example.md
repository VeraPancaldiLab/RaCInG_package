# SKCM melanoma example input for RaCInG

Pre-processed input matrices derived from TCGA skin cutaneous melanoma
(SKCM) bulk RNA-seq data. The object is a named list that can be passed
directly to
[`compute_racing_kernel()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/compute_racing_kernel.md)
or
[`compute_racing_montecarlo()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/compute_racing_montecarlo.md)
via the `input_data` parameter.

## Usage

``` r
skcm_example
```

## Format

A named list with the following elements:

- Lmatrix:

  Numeric matrix (9 cell types x 276 ligands). Binary cell-to-ligand
  compatibility.

- Rmatrix:

  Numeric matrix (9 cell types x 298 receptors). Binary cell-to-receptor
  compatibility.

- Cmatrix:

  Numeric matrix (10 patients x 9 cell types). Row-normalised cell-type
  abundance estimates from deconvolution, with M1 and M2 macrophages
  merged into a single M category.

- LRmatrix:

  3-D numeric array (276 ligands x 298 receptors x 10 patients).
  Normalised ligand–receptor interaction probability tensor.

- celltypes:

  Character vector of 9 cell-type names (alphabetically sorted).

- ligands:

  Character vector of 276 ligand names.

- receptors:

  Character vector of 298 receptor names.

- Sign_matrix:

  Numeric matrix (276 x 298) of zeros (unknown interaction signs).

## Source

Derived from TCGA SKCM data processed with TMEmod deconvolution and
OmniPath ligand–receptor annotations.

## Examples

``` r
data(skcm_example)
str(skcm_example, max.level = 1)
#> List of 8
#>  $ Lmatrix    : num [1:9, 1:276] 1 0 1 1 1 1 1 1 0 1 ...
#>   ..- attr(*, "dimnames")=List of 2
#>  $ Rmatrix    : num [1:9, 1:298] 1 0 1 1 0 1 1 1 0 0 ...
#>   ..- attr(*, "dimnames")=List of 2
#>  $ Cmatrix    : num [1:10, 1:9] 0.00374 0.01407 0.02119 0 0.00313 ...
#>   ..- attr(*, "dimnames")=List of 2
#>  $ LRmatrix   : num [1:276, 1:298, 1:10] 0.000228 0 0 0 0 ...
#>  $ celltypes  : chr [1:9] "B" "CAF" "CD8+ T" "DC" ...
#>  $ ligands    : chr [1:276] "LGALS9" "ADAM10" "TNFSF12" "ICOSLG" ...
#>  $ receptors  : chr [1:298] "PTPRC" "MET" "CD44" "LRP1" ...
#>  $ Sign_matrix: num [1:276, 1:298] 0 0 0 0 0 0 0 0 0 0 ...

# Use directly with the kernel workflow
# result <- compute_racing_kernel(input_data = skcm_example)
```
