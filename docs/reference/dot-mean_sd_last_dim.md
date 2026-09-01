# Mean and SD across the last dimension of an array, vectorized

Equivalent to `apply(arr, seq_len(length(dim(arr)) - 1), mean)` /
`apply(arr, ..., stats::sd)`, but avoids calling
[`mean()`](https://rdrr.io/r/base/mean.html)/[`sd()`](https://rdrr.io/r/stats/sd.html)
once per element of the other dimensions – for a `[264,264,264,itNo]`
tensor that's ~18.4 million individual R function calls per statistic
with [`apply()`](https://rdrr.io/r/base/apply.html), which takes well
over a minute regardless of `itNo`. Reshaping to a matrix and using
[`rowMeans()`](https://rdrr.io/r/base/colSums.html)/a vectorized
variance formula does the same computation in a fraction of a second.

## Usage

``` r
.mean_sd_last_dim(arr)
```

## Arguments

- arr:

  Array whose last dimension is averaged over.

## Value

A list with `mean` and `sd`, each an array with one fewer dimension than
`arr` (the last dimension collapsed).
