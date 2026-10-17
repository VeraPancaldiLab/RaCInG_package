# Warn about cell-type names outside multideconv's nomenclature

A name is accepted if, after dropping a `"_Subgroup.<n>"` suffix, it is
a multideconv cell type or ends in `"_<cell type>"` (e.g.
`"Quantiseq_TIL10_B.cells"`).

## Usage

``` r
.warn_nonstandard_celltypes(celltypes)
```

## Arguments

- celltypes:

  Character vector of cell-type (column) names.

## Value

Invisibly, the non-standard names.
