# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working
with code in this repository.

## Project Overview

**RaCInG** (Random Cell-cell Interaction Generator) is an R package for
reconstructing patient-specific cell-cell communication networks from
bulk RNA-seq data. It implements two main workflows: a deterministic
**kernel method** and a stochastic **Monte Carlo simulation** method.

Reference: van Santvoort et al. (2025, Cell Reports Methods)

## Development Commands

All standard R package development is done via `devtools` in an R
session:

``` r
devtools::load_all()        # Load package during development (replaces install)
devtools::check()           # Run R CMD check (equivalent to CI)
devtools::document()        # Regenerate NAMESPACE and man/*.Rd from roxygen2 comments
roxygen2::roxygenise()      # Alternative to devtools::document()
pkgdown::build_site()       # Build the documentation website locally
```

CI/CD runs `R CMD check --no-manual --compact-vignettes=gs+qpdf` on
push/PR via GitHub Actions (`.github/workflows/R-CMD-check.yaml`). Docs
are auto-deployed via `.github/workflows/pkgdown.yaml`.

``` r
devtools::test()            # Run the testthat regression suite (tests/testthat/)
```

The `testthat` suite (`tests/testthat/test-*.R`) covers correctness
regressions across every module (network generation, kernel/Monte Carlo
feature extraction, statistical analysis, txt/csv I/O, utilities) — run
it after any change to `R/`, not just `R CMD check`.

## Architecture

### Two Core Workflows

**1. Kernel Method** (`R/Kernel_Method.R`) - Entry point:
`compute_racing_kernel(counts, ...)` or
`compute_racing_kernel(input_data = ..., ...)` - Computes a
deterministic 3D kernel `[sender_cell × receiver_cell × patient]` via
weighted matrix multiplication:
`lig_weight %*% LRmatrix %*% t(rec_weight)`. The kernel is computed
exactly once per call regardless of how many feature types are
requested. - Extracts network motif features from the kernel based on
`communication_type` (`"D"`, `"W"`, `"TT"` — returns both trust- and
cycle-triangle columns — or `"GSCC"`), which also accepts a character
vector to extract several types from the same already-computed kernel in
one call (`features` is then a named list, one data frame per type) -
Feature columns with no possible ligand-receptor pathway in any patient
(structurally zero for every patient) are dropped automatically rather
than returned as dead all-zero columns - Returns
`list(kernel, kernel_norm, features)`

**2. Monte Carlo Method** (`R/Monte_Carlo_Method.R`) - Entry point:
`compute_racing_montecarlo(input_data = ..., communication_type, Ncells, Ngraphs, Ndegree, ncores, ...)` -
For each patient, simulates `Ngraphs` random directed graph realizations
(each with `Ncells` cells and target average degree `Ndegree`) via
[`model1()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/model1.md)
(respects cell-type abundances & LR probabilities), then counts motifs
on each realization - `communication_type` accepts a vector (`"D"`,
`"W"`, `"TT"`, `"CT"`, `"GSCC"`); when more than one is given,
[`countAllTypes()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/countAllTypes.md)
extracts every requested type from the *same* simulated graphs instead
of resimulating per type (`output` is then a named list, one entry per
type) - `ncores` parallelizes across patients (independent of each
other) via `parallel`/`doParallel`/`foreach` — portable across
Windows/macOS/Linux; workers are separate processes, so they use the
*installed* package, not uncommitted
[`source()`](https://rdrr.io/r/base/source.html)-loaded changes -
Outputs mean/std of motif counts across realizations; optionally
normalizes against a uniform random baseline - Results written to
`.out`/CSV files

### Data Flow

    Raw Counts (genes × samples)
        ↓ prepare_input_files()   [R/RaCInG_input_generation.R]
        ├─ Lmatrix  (cells × ligands)
        ├─ Rmatrix  (cells × receptors)
        ├─ Cmatrix  (patients × cells, normalized abundances)
        └─ LRmatrix (ligands × receptors × patients)
        ↓
    Kernel Method                  Monte Carlo Method
    compute_kernel()               runSim() → model1() + motif counters
    compute_kernel_features()      compute_results_processing()
        ↓
    wilcox_group_test() → top_features_plot()

### Module Responsibilities

| File                        | Responsibility                                                                                                                                                                                                                                                                                     |
|-----------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `RaCInG_input_generation.R` | Preprocessing raw counts → standardized input matrices; CSV I/O                                                                                                                                                                                                                                    |
| `Kernel_Method.R`           | Deterministic kernel computation and feature extraction                                                                                                                                                                                                                                            |
| `Monte_Carlo_Method.R`      | Monte Carlo simulation orchestration                                                                                                                                                                                                                                                               |
| `network_generation.R`      | Random directed graph generation ([`model1()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/model1.md))                                                                                                                                                                              |
| `feature_extraction.R`      | Graph motif enumeration (wedges, triangles, GSCC)                                                                                                                                                                                                                                                  |
| `Utilities.R`               | Edge list ↔︎ adjacency matrix conversion; motif helpers                                                                                                                                                                                                                                             |
| `statistical_analysis.R`    | Wilcoxon group tests + ranked bar plots ([`top_features_plot()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/top_features_plot.md)), Spearman correlation + rainfall plots ([`correlation_plot()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/correlation_plot.md)) |
| `txt_to_csv.R`              | Read/write simulation outputs                                                                                                                                                                                                                                                                      |
| `data.R`                    | Documentation for bundled `skcm_example` dataset                                                                                                                                                                                                                                                   |

### Key Design Decisions

- **Sparse matrices**:
  [`Matrix::sparseMatrix()`](https://rdrr.io/pkg/Matrix/man/sparseMatrix.html)
  is used in Monte Carlo for memory-efficient adjacency matrices.
- **Optional heavy dependencies**:
  [`prepare_input_files()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/prepare_input_files.md)
  requires ADImpute, multideconv, OmnipathR, and liana — but all
  workflows accept pre-computed input matrices directly, making these
  optional.
- **Normalization**: Both the kernel and Monte Carlo methods can
  optionally compute a uniformized baseline (`norm`/`normalize`) to
  express features as an abundance-independent enrichment ratio, instead
  of raw communication strength.
- **Communication types**: `communication_type` controls which motif is
  extracted — `"D"` (direct), `"W"` (wedge/2-path), `"TT"` (trust
  triangle), `"CT"` (cycle triangle, Monte Carlo only as its own type —
  the kernel method returns it bundled with `"TT"`), `"GSCC"` (global
  strongly connected component). Both
  [`compute_racing_kernel()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/compute_racing_kernel.md)
  and
  [`compute_racing_montecarlo()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/compute_racing_montecarlo.md)
  accept a vector here to extract multiple types without redoing the
  expensive step (kernel computation / graph simulation) per type.
- **Vectorized, not looped, feature extraction**:
  [`calculateWedges()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/calculateWedges.md)/[`computeTriangles()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/computeTriangles.md)
  (kernel) and the array-reshape logic in
  [`compute_results_processing()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/compute_results_processing.md)
  (Monte Carlo) build the full combination grid up front and compute
  over it with matrix operations, instead of a triple-nested loop
  growing an R list one column at a time — the latter pattern is what
  makes cell-type counts in the hundreds (where wedge/triangle
  combinations reach 10⁵⁻¹⁰6) crash or hang.
- **Column naming is shared across methods**: kernel and Monte Carlo
  feature matrices use identical column-naming conventions (`"Dir_A_B"`,
  `"W_A_B_C"`, etc.) and real patient IDs as rownames, so their outputs
  can be fed into the same downstream statistical/plotting functions
  interchangeably.

### Example Data

`skcm_example` (in `data/skcm_example.rda`) is pre-processed TCGA skin
cutaneous melanoma data: 10 patients × 9 cell types × 276 ligands × 298
receptors. Use it for testing and examples.

## Documentation

- Add/edit roxygen2 comments in `R/` files, then run
  [`devtools::document()`](https://devtools.r-lib.org/reference/document.html)
  to regenerate `man/` and `NAMESPACE`.
- The vignette is in `vignettes/RaCiNG.Rmd` and provides the canonical
  usage walkthrough.
- Website structure is defined in `_pkgdown.yml` (Bootstrap 5,
  auto-deployed to GitHub Pages).
