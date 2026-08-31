# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**RaCInG** (Random Cell-cell Interaction Generator) is an R package (v0.1.0) for reconstructing patient-specific cell-cell communication networks from bulk RNA-seq data. It implements two main workflows: a deterministic **kernel method** and a stochastic **Monte Carlo simulation** method.

Reference: van Santvoort et al. (2025, Cell Reports Methods)

## Development Commands

All standard R package development is done via `devtools` in an R session:

```r
devtools::load_all()        # Load package during development (replaces install)
devtools::check()           # Run R CMD check (equivalent to CI)
devtools::document()        # Regenerate NAMESPACE and man/*.Rd from roxygen2 comments
roxygen2::roxygenise()      # Alternative to devtools::document()
pkgdown::build_site()       # Build the documentation website locally
```

CI/CD runs `R CMD check --no-manual --compact-vignettes=gs+qpdf` on push/PR via GitHub Actions (`.github/workflows/R-CMD-check.yaml`). Docs are auto-deployed via `.github/workflows/pkgdown.yaml`.

There is no test suite beyond `R CMD check` — correctness is validated through the vignette and example data.

## Architecture

### Two Core Workflows

**1. Kernel Method** (`R/Kernel_Method.R`)
- Entry point: `compute_racing_kernel(counts, ...)` or `compute_racing_kernel(input_list, ...)`
- Computes a deterministic 3D kernel `[sender_cell × receiver_cell × patient]` via weighted matrix multiplication: `lig_weight %*% LRmatrix %*% t(rec_weight)`
- Extracts network motif features from the kernel based on `communication_type` ("D", "W", "TT", "GSCC")
- Returns `list(kernel, kernel_norm, features)`

**2. Monte Carlo Method** (`R/Monte_Carlo_Method.R`)
- Entry point: `compute_racing_montecarlo(input_list, Ngraphs, Niterations, ...)`
- For each patient × Niterations, generates `Ngraphs` random directed graphs via `model1()` (respects cell-type abundances & LR probabilities), then counts motifs
- Outputs mean/std of motif counts across iterations; optionally normalizes against a uniform random baseline
- Results written to CSV files

### Data Flow

```
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
```

### Module Responsibilities

| File | Responsibility |
|------|---------------|
| `RaCInG_input_generation.R` | Preprocessing raw counts → standardized input matrices; CSV I/O |
| `Kernel_Method.R` | Deterministic kernel computation and feature extraction |
| `Monte_Carlo_Method.R` | Monte Carlo simulation orchestration |
| `network_generation.R` | Random directed graph generation (`model1()`) |
| `feature_extraction.R` | Graph motif enumeration (wedges, triangles, GSCC) |
| `Utilities.R` | Edge list ↔ adjacency matrix conversion; motif helpers |
| `statistical_analysis.R` | Wilcoxon tests and volcano plots |
| `txt_to_csv.R` | Read/write simulation outputs |
| `data.R` | Documentation for bundled `skcm_example` dataset |

### Key Design Decisions

- **Sparse matrices**: `Matrix::sparseMatrix()` is used in Monte Carlo for memory-efficient adjacency matrices.
- **Optional heavy dependencies**: `prepare_input_files()` requires ADImpute, multideconv, OmnipathR, and liana — but all workflows accept pre-computed input matrices directly, making these optional.
- **Normalization**: The kernel method can optionally compute a uniformized baseline kernel to normalize for varying LR interaction strengths across patients.
- **Communication types**: `communication_type` parameter controls which motif is extracted — "D" (direct), "W" (wedge/2-path), "TT" (trust triangle/3-cycle), "GSCC" (global strong connected component).

### Example Data

`skcm_example` (in `data/skcm_example.rda`) is pre-processed TCGA skin cutaneous melanoma data: 10 patients × 9 cell types × 276 ligands × 298 receptors. Use it for testing and examples.

## Documentation

- Add/edit roxygen2 comments in `R/` files, then run `devtools::document()` to regenerate `man/` and `NAMESPACE`.
- The vignette is in `vignettes/RaCiNG.Rmd` and provides the canonical usage walkthrough.
- Website structure is defined in `_pkgdown.yml` (Bootstrap 5, auto-deployed to GitHub Pages).
