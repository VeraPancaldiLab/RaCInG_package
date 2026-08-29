# Estimate per-cell-type average expression profiles

For each gene, solves a non-negative least squares regression of bulk
TPM expression across samples against per-sample cell-type fractions,
giving one expression value per gene per cell type:
`bulk_expr[gene, sample] ~= sum_celltype fractions[sample, celltype] * profile[gene, celltype]`.

## Usage

``` r
.estimate_expression_profiles(counts_tpm, deconv)
```

## Arguments

- counts_tpm:

  Gene-by-sample TPM matrix.

- deconv:

  Sample-by-cell-type fraction matrix.

## Value

A gene-by-cell-type data frame of estimated expression profiles.
