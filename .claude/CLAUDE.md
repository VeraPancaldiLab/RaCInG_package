# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**RaCInG** (Random Cell-cell Interaction Generator) is an R package for reconstructing patient-specific cell-cell communication networks from bulk RNA-seq data. It implements two main workflows: a deterministic **kernel method** and a stochastic **Monte Carlo simulation** method.

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

```r
devtools::test()            # Run the testthat regression suite (tests/testthat/)
```

The `testthat` suite (`tests/testthat/test-*.R`) covers correctness regressions across every module (network generation, kernel/Monte Carlo feature extraction, statistical analysis, txt/csv I/O, utilities) — run it after any change to `R/`, not just `R CMD check`.

## Architecture

### Two Core Workflows

**1. Kernel Method** (`R/Kernel_Method.R`)
- Entry point: `compute_racing_kernel(counts, ...)` or `compute_racing_kernel(input_data = ..., ...)`
- Computes a deterministic 3D kernel `[sender_cell × receiver_cell × patient]` via weighted matrix multiplication: `lig_weight %*% LRmatrix %*% t(rec_weight)`.
- Extracts direct communication features only (`communication_type = "D"`); any other type errors. Wedges and triangles were removed from the package, and GSCC is Monte Carlo only.
- Feature columns with no possible ligand-receptor pathway in any patient (structurally zero for every patient) are dropped automatically rather than returned as dead all-zero columns
- Returns `list(kernel, kernel_norm, features)`

**2. Monte Carlo Method** (`R/Monte_Carlo_Method.R`)
- Entry point: `compute_racing_montecarlo(input_data = ..., communication_type, Ncells, Ngraphs, Ndegree, ncores, ...)`
- For each patient, simulates `Ngraphs` random directed graph realizations (each with `Ncells` cells and target average degree `Ndegree`) via `model1()` (respects cell-type abundances & LR probabilities), then counts direct edges and/or GSCC membership on each realization
- `communication_type` is `"D"`, `"GSCC"`, or `c("D", "GSCC")`; with both, `countAllTypes()` extracts them from the *same* simulated graphs instead of resimulating per type (`output` is then a named list, one entry per type). Wedges/triangles were removed.
- `ncores` parallelizes across patients (independent of each other) via `parallel`/`doParallel`/`foreach` — portable across Windows/macOS/Linux; workers are separate processes, so they use the *installed* package, not uncommitted `source()`-loaded changes
- Outputs mean/std across realizations; optionally normalizes against a uniform random baseline
- Results written to `.out`/CSV files

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
compute_kernel()               runSim() → model1() + D/GSCC counters
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
| `LR_Contributions.R` | Per-ligand-receptor-pair decomposition of direct scores (`computeLRContributions()`) + dot plot (`plotLRContributions()`) |
| `feature_extraction.R` | Strongly connected components / GSCC (Monte Carlo) |
| `Utilities.R` | Edge list ↔ adjacency matrix conversion |
| `statistical_analysis.R` | Wilcoxon group tests + ranked bar plots (`top_features_plot()`), Spearman correlation + rainfall plots (`correlation_plot()`) |
| `txt_to_csv.R` | Read/write simulation outputs |
| `data.R` | Documentation for bundled `skcm_example` dataset |

### Key Design Decisions

- **Sparse matrices**: `Matrix::sparseMatrix()` is used in Monte Carlo for memory-efficient adjacency matrices.
- **Optional heavy dependencies**: `prepare_input_files()` requires ADImpute, multideconv, OmnipathR, and liana — but all workflows accept pre-computed input matrices directly, making these optional. multideconv must be >= 0.2.0: `.estimate_expression_profiles_subgroups()` relies on its single-pass subgroups (members are always raw columns) and unchanged feature names. The cell-type vocabulary lives in `multideconv::get_cell_type_nomenclature()` — don't copy it into RaCInG.
- **Normalization**: Both the kernel and Monte Carlo methods can optionally compute a uniformized baseline (`norm`/`normalize`) to express features as an abundance-independent enrichment ratio, instead of raw communication strength.
- **Communication types**: kernel = `"D"` only; Monte Carlo = `"D"` and/or `"GSCC"`. Wedges (`"W"`) and triangles (`"TT"`/`"CT"`), and the kernel's analytic GSCC, were deliberately removed — don't reintroduce them. Unsupported types error via `.check_kernel_communication_type()` / `.check_montecarlo_communication_type()`.
- **Column naming**: direct features are `"Dir_A_B"` (lower-index cell type first) in both methods, GSCC features `"GSCC_A"`; rownames are real patient IDs, or `Patient_<row>` when `Cmatrix` has no rownames (as in `skcm_example`).

### Example Data

`skcm_example` (in `data/skcm_example.rda`) is pre-processed TCGA skin cutaneous melanoma data: 10 patients × 9 cell types × 276 ligands × 298 receptors. Use it for testing and examples.

### Verified Against the Original Python Implementation

The companion project at `~/Documents/Racing_paper` (`scripts/Melanoma_RaCInG_R_portcheck.qmd`,
`scripts/RaCInG_original_melanoma.ipynb`, `scripts/Compare_RaCInG_R_vs_original.qmd`) verified
that `compute_kernel()`/`calculateDirect()` (and the since-removed Wedges/Triangles functions)
reproduce the original, unmodified Python implementation's math exactly on the full 336-patient
melanoma cohort — Pearson r = 1 and relative differences at machine precision. Note the
original's `calculateDirect()` internally truncates `"CD8+ T"` → `"CD8"` before building column
labels; this R port does not. Rerunning the Wedges/Triangles part of that check needs a
package version from before their removal.

### TODO: Finish or Remove the Sign Matrix

The interaction-sign feature (+1 stimulatory / -1 inhibitory) is plumbed in but not
functional: `signed = TRUE` (in `prepare_input_files()`, `compute_racing_kernel()`,
`compute_racing_montecarlo()`) makes `generateInput(read_signs = TRUE)` require a
`Sign_matrix_<file_name>.csv` that nothing writes, so it always errors unless the file is
placed by hand; and the loaded `Sign_matrix` is never used -- the kernel ignores it, and
`runSim()`/`countAllTypes()` never pass it to `model1(Signmatrix = )`, whose extra edge-sign
column nothing reads. Decide whether to implement it end to end (build the sign file from the
prior network, use signs in the features) or remove `signed`/`read_signs`/
`Read_Lig_Rec_Interaction()`/`Sign_matrix`/`model1()`'s `Signmatrix`.

### Known Design Gap: Feature Column Identity

`calculateDirect()` returns sender/receiver cell-type
identity only baked into a glued column-name string (`"Dir_A_B"`) — never as
structured columns (e.g. `CellA`/`CellB`) alongside the values. Every downstream caller that
needs the identity separately (e.g. "does this pair involve cell type X, in either position")
has to reconstruct it by rebuilding the same name in the same internal iteration order and
matching against `colnames()` — this exact reconstruction has been independently reimplemented
at least three times in the companion `Racing_paper` project. Consider having these functions
return the pair identity as real columns (or an attached attribute) directly, so callers don't
need to reverse-engineer a string this package itself constructed.

## Documentation

- Add/edit roxygen2 comments in `R/` files, then run `devtools::document()` to regenerate `man/` and `NAMESPACE`.
- The vignette is in `vignettes/RaCiNG.Rmd` and provides the canonical usage walkthrough.
- Website structure is defined in `_pkgdown.yml` (Bootstrap 5, auto-deployed to GitHub Pages).
