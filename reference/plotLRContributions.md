# Dot plot of top ligand-receptor pair contributions (LIANA-style)

Ranks ligand-receptor pairs by mean contribution (within each `Edge`
separately, since different edges of a Wedge/Triangle are typically
driven by different pairs) and plots the top ones as a dot per patient x
pair: dot size is `PctOfTotal`, dot color is the raw `Contribution`.

## Usage

``` r
plotLRContributions(lr_contrib, top_n = 15, title = NULL)
```

## Arguments

- lr_contrib:

  Output of
  [`computeLRContributions()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/computeLRContributions.md),
  [`computeWedgeLRContributions()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/computeWedgeLRContributions.md),
  or
  [`computeTriangleLRContributions()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/computeTriangleLRContributions.md).

- top_n:

  Number of top ligand-receptor pairs to show per `Edge`.

- title:

  Plot title. Defaults to `NULL`, which auto-generates one naming the
  cell type(s) involved (e.g.
  `"Top ligand-receptor pairs: CLL - NKcells"`), derived from
  `lr_contrib`'s `Sender`/`Receiver` columns.

## Value

A `ggplot2` object.
