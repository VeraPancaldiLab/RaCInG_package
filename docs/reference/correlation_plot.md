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
corr <- correlate_features_with_score(feats, score)
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
if (requireNamespace("ggplot2", quietly = TRUE)) {
  correlation_plot(corr, top_n = 5)
}
```
