# Estimate per-cell-type average expression profiles from one composition

For each gene, solves a non-negative least squares regression of bulk
TPM expression across samples against per-sample cell-type fractions,
giving one expression value per gene per cell type:
`bulk_expr[gene, sample] ~= sum_celltype fractions[sample, celltype] * profile[gene, celltype]`.
Valid only when every column of `deconv` already belongs to one coherent
composition (fractions summing to ~1 per sample). When `deconv` instead
spans multiple deconvolution methods/signatures concatenated together
(as
[`multideconv::compute.deconvolution()`](https://verapancaldilab.github.io/multideconv/reference/compute.deconvolution.html)
produces before subgroup analysis), use
[`.estimate_expression_profiles_subgroups()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/dot-estimate_expression_profiles_subgroups.md)
instead – fitting NNLS directly across method/signature blocks here
would violate the mixture equation this regression assumes.

## Usage

``` r
.estimate_expression_profiles_samples(counts_tpm, deconv)
```

## Arguments

- counts_tpm:

  Gene-by-sample TPM matrix.

- deconv:

  Sample-by-cell-type fraction matrix.

## Value

A gene-by-cell-type data frame of estimated expression profiles.
