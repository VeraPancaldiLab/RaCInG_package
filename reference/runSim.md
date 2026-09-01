# Run Monte Carlo simulations for one or more patients

Run Monte Carlo simulations for one or more patients

## Usage

``` r
runSim(
  Lmatrix,
  Rmatrix,
  Cmatrix,
  LRmatrix,
  cells,
  communication_type,
  pats = "all",
  N = 10000,
  itNo = 100,
  av = 20,
  output_folder = NULL,
  file.name = NULL,
  norm = FALSE,
  patient_idx = NULL,
  ncores = 1
)
```

## Arguments

- Lmatrix:

  Cell-by-ligand compatibility matrix.

- Rmatrix:

  Cell-by-receptor compatibility matrix.

- Cmatrix:

  Patient-by-cell-type abundance matrix.

- LRmatrix:

  Ligand-receptor-by-patient tensor.

- cells:

  Character vector of cell-type names.

- communication_type:

  Feature family to simulate (`"D"`, `"W"`, `"TT"`, `"CT"`, or
  `"GSCC"`), or a character vector of several of these. When more than
  one is given, every requested type is extracted from the *same*
  simulated graphs (via
  [`countAllTypes()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/countAllTypes.md))
  instead of re-simulating a fresh set of graphs per type, and one
  `.out` file per type is written (file names suffixed with the type,
  e.g. `<file.name>_D.out`). With a single type, output naming is
  unchanged from previous versions (`<file.name>.out`, no suffix).

- pats:

  Number of patients to process, or `"all"`.

- N:

  Number of cells per graph.

- itNo:

  Number of Monte Carlo iterations.

- av:

  Target average degree.

- output_folder:

  Directory used to write the `.out` files.

- file.name:

  Output filename stem.

- norm:

  Logical; if `TRUE`, use a uniformized LR baseline.

- patient_idx:

  Optional single patient index to simulate.

- ncores:

  Number of cores to compute patients on in parallel, via
  [`parallel::makeCluster()`](https://rdrr.io/r/parallel/makeCluster.html) +
  [`doParallel::registerDoParallel()`](https://rdrr.io/pkg/doParallel/man/registerDoParallel.html) +
  `foreach::foreach(...) %dopar% {...}` (the same backend used
  throughout `pipeML`), so it works identically on Windows/macOS/Unix.
  `ncores = 1` (the default) runs sequentially via
  [`lapply()`](https://rdrr.io/r/base/lapply.html) and skips cluster
  setup entirely. Patients are independent of each other, so this
  parallelizes near-linearly. File writing always happens sequentially
  afterward, in patient order, to keep the on-disk format unchanged.
  Because cluster workers are separate R processes (not forks), each one
  loads the installed `RaCInG` package rather than inheriting the
  calling session's state – if you are iterating on package source via
  [`source()`](https://rdrr.io/r/base/source.html) instead of
  [`library(RaCInG)`](https://github.com/VeraPancaldiLab/RaCInG_package),
  reinstall the package first so workers see your latest changes.

## Value

Invisibly writes the simulation outputs to disk.
