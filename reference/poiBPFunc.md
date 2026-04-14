# Poisson branching-process helper

Poisson branching-process helper

## Usage

``` r
poiBPFunc(x, M, sens)
```

## Arguments

- x:

  Current estimate vector.

- M:

  Mean offspring matrix.

- sens:

  Length of `x`.

## Value

A numeric vector used by
[`nleqslv::nleqslv()`](https://bertcarnell.github.io/nleqslv/reference/nleqslv.html).
