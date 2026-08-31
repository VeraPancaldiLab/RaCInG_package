# Ranked bar plot of top up/down features from Wilcoxon results

Shows the top most upregulated and top most downregulated significant
features (by log2 fold change) from
[`wilcox_group_test()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/wilcox_group_test.md)
as a ranked horizontal bar chart. This scales to features with very long
names: each name is drawn once as an axis label rather than needing to
be positioned/repelled on a scatter plot, so labels never get dropped or
overlap regardless of how many features there are.

## Usage

``` r
top_features_plot(
  wilcox_results,
  top_n = 10,
  p_threshold = 0.05,
  title = "Top features"
)
```

## Arguments

- wilcox_results:

  Output of
  [`wilcox_group_test()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/wilcox_group_test.md).
  If it contains more than one `Comparison`, only the first is plotted;
  subset beforehand to plot a specific comparison.

- top_n:

  Number of top upregulated and top downregulated features to show.

- p_threshold:

  Only features with `Adjusted_P_value < p_threshold` are considered.

- title:

  Plot title.

## Value

A `ggplot2` object.
