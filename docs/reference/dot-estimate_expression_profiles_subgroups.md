# Estimate per-cell-type expression profiles across multiple deconvolution methods/signatures

[`multideconv::compute.deconvolution()`](https://verapancaldilab.github.io/multideconv/reference/compute.deconvolution.html)'s
raw output concatenates one block of columns per method-signature
combination (named `"{method}_{signature}_{celltype}"`). Each block is,
on its own, a valid composition – multideconv itself verifies every
block's standardized cell-type columns sum to 1 per sample – but
[`multideconv::compute.deconvolution.analysis()`](https://verapancaldilab.github.io/multideconv/reference/compute.deconvolution.analysis.html)'s
merged "Deconvolution matrix" is not: it stitches together
representatives from different blocks per cell type, so fitting
[`.estimate_expression_profiles_samples()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/dot-estimate_expression_profiles_samples.md)'s
NNLS mixture equation directly on it is invalid. This instead fits one
profile per method-signature block (always a valid composition), then
merges the resulting *profiles* – never the raw fractions – according to
whichever subgroups `subgroups` (already computed by the caller via
[`multideconv::compute.deconvolution.analysis()`](https://verapancaldilab.github.io/multideconv/reference/compute.deconvolution.analysis.html))
forms across those blocks.

## Usage

``` r
.estimate_expression_profiles_subgroups(counts_tpm, deconv_wide, subgroups)
```

## Arguments

- counts_tpm:

  Gene-by-sample TPM matrix.

- deconv_wide:

  Sample-by-`{method}_{signature}_{celltype}` raw fraction matrix –
  [`multideconv::compute.deconvolution()`](https://verapancaldilab.github.io/multideconv/reference/compute.deconvolution.html)'s
  direct output, before subgroup analysis.

- subgroups:

  The list returned by
  `multideconv::compute.deconvolution.analysis(deconv_wide)`.

## Value

A gene-by-cell-type data frame of estimated expression profiles, one
column per feature surviving in `subgroups[["Deconvolution matrix"]]`.

## Details

Requires multideconv \>= 0.2.0, whose subgroups come from a single
complete-linkage clustering (members are always raw `deconv_wide`
columns, never other subgroups) and keep feature names unchanged, so
subgroup members and surviving raw features match
`colnames(deconv_wide)` as-is. Each subgroup's profile is the median of
its members' profiles, mirroring how multideconv builds the subgroup
fraction itself.
