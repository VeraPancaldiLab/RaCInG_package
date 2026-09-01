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
