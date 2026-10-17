# Correlate network features with an external score

Computes the Spearman correlation between an external per-patient score
(e.g. an immune response score) and every feature column, pooled across
all patients ("All") and, optionally, within each group separately.

## Usage

``` r
correlate_features_with_score(
  feature_matrix,
  score,
  group = NULL,
  p_adjust_method = "fdr"
)
```

## Arguments

- feature_matrix:

  Numeric matrix or data frame with patients in rows (identified by
  [`rownames()`](https://rdrr.io/r/base/colnames.html)) and features in
  columns.

- score:

  Named numeric vector of per-patient scores. Names must match
  `rownames(feature_matrix)`; patients missing a score are dropped.

- group:

  Optional grouping vector (e.g. cancer type or cohort), in the same row
  order as `feature_matrix`. When supplied, correlations are also
  computed separately within each group.

- p_adjust_method:

  Multiple-testing correction method passed to
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html).

## Value

A data frame with `Group`, `Feature`, `Rho`, `P_value`, and
`Adjusted_P_value`.

## Examples

``` r
data(skcm_example)
feats <- compute_racing_kernel(input_data = skcm_example)$features
#> Using pre-computed input matrices; skipping input generation.
#> Computing kernel for 10 patients
#> Calculating kernel...
#> Calculating features...
set.seed(1)
score <- setNames(rnorm(nrow(feats)), rownames(feats))
head(correlate_features_with_score(feats, score))
#> Warning: Cannot compute exact p-value with ties
#> Warning: Cannot compute exact p-value with ties
#> Warning: Cannot compute exact p-value with ties
#> Warning: Cannot compute exact p-value with ties
#> Warning: Cannot compute exact p-value with ties
#> Warning: Cannot compute exact p-value with ties
#> Warning: Cannot compute exact p-value with ties
#> Warning: Cannot compute exact p-value with ties
#> Warning: Cannot compute exact p-value with ties
#> Warning: Cannot compute exact p-value with ties
#> Warning: Cannot compute exact p-value with ties
#> Warning: Cannot compute exact p-value with ties
#> Warning: Cannot compute exact p-value with ties
#> Warning: Cannot compute exact p-value with ties
#> Warning: Cannot compute exact p-value with ties
#> Warning: Cannot compute exact p-value with ties
#> Warning: Cannot compute exact p-value with ties
#> Warning: Cannot compute exact p-value with ties
#> Warning: Cannot compute exact p-value with ties
#> Warning: Cannot compute exact p-value with ties
#> Warning: Cannot compute exact p-value with ties
#> Warning: Cannot compute exact p-value with ties
#> Warning: Cannot compute exact p-value with ties
#> Warning: Cannot compute exact p-value with ties
#>   Group      Feature        Rho     P_value Adjusted_P_value
#> 1   All   Dir_B_Endo -0.8267515 0.003179714       0.07154357
#> 2   All     Dir_B_NK -0.8267515 0.003179714       0.07154357
#> 3   All    Dir_B_CAF -0.7538029 0.011794786       0.13269135
#> 4   All  Dir_B_Tumor -0.7659610 0.009787040       0.13269135
#> 5   All      Dir_B_B -0.6565380 0.039204386       0.32126842
#> 6   All Dir_B_CD8+ T -0.6319494 0.049975087       0.32126842
```
