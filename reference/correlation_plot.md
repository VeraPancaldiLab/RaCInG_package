# Rainfall plot of feature-score correlations

Shows the top positive and top negative correlations from
[`correlate_features_with_score()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/correlate_features_with_score.md)
as a diverging bar chart.

## Usage

``` r
correlation_plot(
  correlation_results,
  group = NULL,
  top_n = 10,
  title = "Feature correlation"
)
```

## Arguments

- correlation_results:

  Output of
  [`correlate_features_with_score()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/correlate_features_with_score.md).

- group:

  Which `Group` value to plot. Defaults to the first group present.

- top_n:

  Number of top positive and top negative correlations to show.

- title:

  Plot title.

## Value

A `ggplot2` object.
