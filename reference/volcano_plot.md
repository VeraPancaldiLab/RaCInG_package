# Create a volcano plot from Wilcoxon results

Create a volcano plot from Wilcoxon results

## Usage

``` r
volcano_plot(
  wilcox_results,
  top_labels = 10,
  p_threshold = 0.05,
  title = "Volcano Plot"
)
```

## Arguments

- wilcox_results:

  Output of
  [`wilcox_group_test()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/wilcox_group_test.md).
  If it contains more than one `Comparison`, only the first is plotted;
  subset beforehand to plot a specific comparison.

- top_labels:

  Number of top significant features to label.

- p_threshold:

  Adjusted p-value threshold used to mark significance.

- title:

  Plot title.

## Value

A `ggplot2` object.
