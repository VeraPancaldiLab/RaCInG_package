# Run Wilcoxon tests across network features

Runs a Wilcoxon rank-sum test per feature (column) between every pair of
groups. Also reports the fold change (ratio of group means) so results
can be plotted as a standard volcano plot with
[`volcano_plot()`](https://mhurtado13.github.io/racing/reference/volcano_plot.md).

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
