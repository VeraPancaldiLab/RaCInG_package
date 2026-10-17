# Run Wilcoxon tests across network features

Runs a Wilcoxon rank-sum test per feature (column) between every pair of
groups. Also reports the fold change (ratio of group means) so results
can be plotted with
[`top_features_plot()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/top_features_plot.md).

## Usage

``` r
wilcox_group_test(data_matrix, groups, p_adjust_method = "fdr")
```

## Arguments

- data_matrix:

  Numeric matrix or data frame with patients in rows and features in
  columns.

- groups:

  Vector of group labels with length matching `nrow(data_matrix)`. When
  more than two groups are present, every pairwise comparison is run.

- p_adjust_method:

  Multiple-testing correction method passed to
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html).

## Value

A data frame with one row per feature per group comparison, including
`Comparison`, `Wilcox_statistic`, `Fold_change`, `Log2FC`, `P_value`,
and `Adjusted_P_value` (adjusted jointly across all comparisons).

## Examples

``` r
data(skcm_example)
feats <- compute_racing_kernel(input_data = skcm_example)$features
#> Using pre-computed input matrices; skipping input generation.
#> Computing kernel for 10 patients
#> Calculating kernel...
#> Calculating features...
groups <- rep(c("A", "B"), length.out = nrow(feats))
res <- wilcox_group_test(feats, groups)
#> Warning: cannot compute exact p-value with ties
#> Warning: cannot compute exact p-value with ties
#> Warning: cannot compute exact p-value with ties
#> Warning: cannot compute exact p-value with ties
#> Warning: cannot compute exact p-value with ties
#> Warning: cannot compute exact p-value with ties
#> Warning: cannot compute exact p-value with ties
#> Warning: cannot compute exact p-value with ties
#> Warning: cannot compute exact p-value with ties
#> Warning: cannot compute exact p-value with ties
#> Warning: cannot compute exact p-value with ties
#> Warning: cannot compute exact p-value with ties
#> Warning: cannot compute exact p-value with ties
#> Warning: cannot compute exact p-value with ties
#> Warning: cannot compute exact p-value with ties
#> Warning: cannot compute exact p-value with ties
#> Warning: cannot compute exact p-value with ties
#> Warning: cannot compute exact p-value with ties
#> Warning: cannot compute exact p-value with ties
#> Warning: cannot compute exact p-value with ties
#> Warning: cannot compute exact p-value with ties
#> Warning: cannot compute exact p-value with ties
#> Warning: cannot compute exact p-value with ties
#> Warning: cannot compute exact p-value with ties
head(res)
#>   Comparison    Feature Wilcox_statistic Fold_change   Log2FC    P_value
#> 1     A_vs_B    Dir_B_B               21    7.395894 2.886724 0.09369262
#> 2     A_vs_B  Dir_B_CAF               19    2.025703 1.018423 0.20869044
#> 3     A_vs_B   Dir_B_DC               18    2.012975 1.009329 0.29480199
#> 4     A_vs_B Dir_B_Endo               18    2.556905 1.354399 0.29480199
#> 5     A_vs_B    Dir_B_M               19    2.542964 1.346511 0.20869044
#> 6     A_vs_B   Dir_B_NK               19    2.264156 1.178974 0.20869044
#>   Adjusted_P_value
#> 1        0.9416225
#> 2        0.9416225
#> 3        0.9416225
#> 4        0.9416225
#> 5        0.9416225
#> 6        0.9416225
```
