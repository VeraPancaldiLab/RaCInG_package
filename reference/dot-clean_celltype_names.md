# Strip method/signature prefixes from deconvolution column names

[`multideconv::compute.deconvolution()`](https://rdrr.io/pkg/multideconv/man/compute.deconvolution.html)
names its columns `"{method}_{signature}_{celltype}"` (e.g.
`"Quantiseq_TIL10_B.cells"`). After
[`multideconv::standardize_celltype_colnames()`](https://rdrr.io/pkg/multideconv/man/standardize_celltype_colnames.html)
has normalized the trailing cell-type portion to multideconv's canonical
vocabulary, this matches the longest canonical cell-type name each
column ends with and drops everything before it, leaving a short,
method-agnostic label. Columns that don't end in a recognized cell-type
name are left unchanged.

## Usage

``` r
.clean_celltype_names(x)
```

## Arguments

- x:

  Character vector of (already standardized) column names.

## Value

A character vector of cleaned cell-type names, the same length as `x`.
