# Legacy GSCC helper

This helper depended on project-specific files from the original
development workflow. The packaged interface now recommends using
[`computeGSCC()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/computeGSCC.md)
directly.

## Usage

``` r
getGSCCAnalytically(cancer, lab = 1, norm = TRUE, test = FALSE)
```

## Arguments

- cancer:

  Legacy dataset identifier.

- lab:

  Scaling factor.

- norm:

  Logical; kept for backward compatibility.

- test:

  Logical; kept for backward compatibility.

## Value

This function stops with a message directing users to
[`computeGSCC()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/computeGSCC.md).
