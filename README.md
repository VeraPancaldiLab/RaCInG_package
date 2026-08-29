RaCInG <img src="man/figures/logo.png" align="right" height="139" />
================

[![R-CMD-check](https://github.com/VeraPancaldiLab/RaCInG_package/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/VeraPancaldiLab/RaCInG_package/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/VeraPancaldiLab/RaCInG_package/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/VeraPancaldiLab/RaCInG_package/actions/workflows/pkgdown.yaml)

**RaCInG** (**Ra**ndom **C**ell-cell **In**teraction **G**enerator)
reconstructs patient-specific cell-cell communication networks from bulk
RNA-seq data and extracts network-level features using either a
kernel-based or Monte Carlo workflow.

This package is the R implementation of the [original RaCInG Python
framework](https://github.com/SysBioOncology/RaCInG), as described in
[van Santvoort et
al. (2025)](https://pubmed.ncbi.nlm.nih.gov/39954673/), repackaged for
seamless integration with R/Bioconductor analysis pipelines. It is
additionally coupled with
[multideconv](https://github.com/VeraPancaldiLab/multideconv) for
deconvolution and cell-type subgroup identification, and with
[OmnipathR](https://omnipathdb.org/)/[liana](https://saezlab.github.io/liana/)
for ligand-receptor prior knowledge.

## Installation

### Install from GitHub

``` r
# install.packages("remotes")
remotes::install_github("VeraPancaldiLab/RaCInG_package")
library(RaCInG)
```

### Install from a local clone

``` r
# install.packages("devtools")
devtools::install(".")
```

### Optional preprocessing dependencies

If you want to start directly from raw counts with
`prepare_input_files()`, install the optional helper packages used
during deconvolution and prior-network construction:

``` r
install.packages(c("ggplot2", "nnls"))
# ADImpute and OmnipathR are available from Bioconductor:
BiocManager::install(c("ADImpute", "OmnipathR"))
# liana and multideconv are GitHub-only:
remotes::install_github(c("saezlab/liana", "VeraPancaldiLab/multideconv"))
```

## Workflow at a glance

| Goal                                 | Function                      | Output                                                  |
|--------------------------------------|-------------------------------|---------------------------------------------------------|
| Build input matrices from raw counts | `prepare_input_files()`       | Named list with `L`, `R`, `C`, `LR` matrices and labels |
| Run deterministic features           | `compute_racing_kernel()`     | Kernel arrays + feature matrix                          |
| Run simulation-based features        | `compute_racing_montecarlo()` | Monte Carlo summaries                                   |
| Compare patient groups               | `wilcox_group_test()`         | Statistics table for downstream plots                   |

## Quick start

``` r
library(RaCInG)

# Build input matrices from a real bulk RNA-seq counts matrix
# (bundled with the optional `multideconv` dependency)
data(raw_counts, package = "multideconv")
counts_matrix <- as.matrix(raw_counts)

input <- prepare_input_files(
  counts = counts_matrix,
  output_folder = "Results/",
  file_name = "example"
)

# Run kernel method (from raw counts)
kernel_res <- compute_racing_kernel(
  counts = counts_matrix,
  file_name = "example",
  output_folder = tempdir(),
  communication_type = "W"
)

# Or pass pre-computed inputs to skip preprocessing
kernel_res <- compute_racing_kernel(
  input_data = input,
  communication_type = "W"
)

# Monte Carlo method
mc_res <- compute_racing_montecarlo(
  input_data = input,
  file_name = "example",
  output_folder = tempdir(),
  communication_type = "W",
  Ncells = 100,
  Ngraphs = 10,
  Ndegree = 3
)

# Compare clinical groups or correlate features with a continuous score
grouping <- c("Responder", "Responder", "Non-responder", "Non-responder")
wilcox_results <- wilcox_group_test(kernel_res$features, grouping)
volcano_plot(wilcox_results)
```

## Documentation

- 📘 Vignette: [Getting started with
  RaCInG](https://VeraPancaldiLab.github.io/RaCInG_package/articles/RaCiNG.html)
- 🌐 Website: <https://VeraPancaldiLab.github.io/RaCInG_package/>
- 🐍 Original Python implementation:
  <https://github.com/SysBioOncology/RaCInG>

## Citation

If you use this package, please cite the RaCInG publication:

> van Santvoort M, Lapuente-Santana Ó, Zopoglou M, Zackl C, Finotello F,
> van der Hoorn P & Eduati F (2025). *Mathematically mapping the network
> of cells in the tumor microenvironment.* Cell Reports Methods, 5(2),
> 100985.

## R Package Development

This R package implementation was developed by [Marcelo
Hurtado](https://github.com/mhurtado13) from the [Pancaldi
team](https://github.com/VeraPancaldiLab), led by Vera Pancaldi. Marcelo
is currently the primary maintainer of the package.
