# RaCInG: Detailed Comparison Report

## R Package ([mhurtado13/racing](https://github.com/mhurtado13/racing)) vs. Original Python Implementation ([SysBioOncology/RaCInG](https://github.com/SysBioOncology/RaCInG))

**Package documentation**: <https://mhurtado13.github.io/racing/>

------------------------------------------------------------------------

## 1. Summary

The R package **RaCInG** builds upon the original Python-based RaCInG
pipeline, repackaging it as a distributable R package. The original
implementation established the core mathematical framework—random graph
generation, kernel-based feature computation, Monte Carlo simulation,
and graph motif extraction—which the R package faithfully preserves
while adapting the architecture for R-native workflows and package
distribution. This report documents every change, adaptation, and
addition in detail.

------------------------------------------------------------------------

## 2. Architectural Changes

### 2.1 From Script Collection to Formal R Package

| Aspect | Original (Python) | R Package |
|----|----|----|
| **Structure** | Python scripts organized in a `Python_Code/` folder (research-oriented layout) | Formal R package with `DESCRIPTION`, `NAMESPACE`, `man/`, `R/`, `data/`, `vignettes/` |
| **Installation** | Clone repository and manage scripts locally | One-command install via `remotes::install_github("mhurtado13/racing")` |
| **Documentation** | Python docstrings and Jupyter notebook demos | Full roxygen2 documentation + pkgdown website + vignettes |
| **Dependencies** | Standard Python packages (`numpy`, `scipy`, `matplotlib`, etc.) managed via `pip`/`conda` | Declared in `DESCRIPTION` with `Imports:` and `Suggests:` fields; automatic dependency resolution |
| **Testing** | `if __name__ == "__main__"` blocks at bottom of scripts | R CMD check compatible; formal check infrastructure |
| **Namespace** | Standard Python module imports | Formal `NAMESPACE` with selective exports and imports |

### 2.2 Language Translation: Python → R

The entire codebase has been rewritten from Python to R. Key translation
decisions include:

- **NumPy arrays** → R matrices and
  [`array()`](https://rdrr.io/r/base/array.html) for 3D/4D structures
- **SciPy sparse matrices** →
  [`Matrix::sparseMatrix`](https://rdrr.io/pkg/Matrix/man/sparseMatrix.html)
  (dgCMatrix)
- **Pandas DataFrames** → R `data.frame` with `dplyr`/`tibble`
  manipulation
- **Python classes** (`informationMicroEnv`) → Standalone R functions
  with explicit parameter passing
- **Python `csv` module** →
  [`utils::read.csv()`](https://rdrr.io/r/utils/read.table.html) and
  [`data.table::fread()`](https://rdrr.io/pkg/data.table/man/fread.html)
- **Matplotlib plotting** → `ggplot2`-based plotting

------------------------------------------------------------------------

## 3. File-by-File Comparison

### 3.1 `network_generation.py` → `R/network_generation.R`

#### Original Python

- Defines a class `informationMicroEnv` that encapsulates cell-type
  distributions, ligand-receptor distributions, and connectivity rules.
- Class methods:
  [`genRandomEdgeList()`](https://mhurtado13.github.io/racing/reference/genRandomEdgeList.md),
  [`genRandomCellTypeList()`](https://mhurtado13.github.io/racing/reference/genRandomCellTypeList.md)
- Imports `distribution_generation` module (separate file) for
  `genRandomCellTypeDistr`, `genRandomLigRecDistr`,
  `genRandomCellLigands`, `genRandomCellReceptors`.
- [`model1()`](https://mhurtado13.github.io/racing/reference/model1.md)
  function creates an `informationMicroEnv` instance and generates a
  graph.

#### R Package Changes

- **Adapted from class-based to functional design**: The
  `informationMicroEnv` class is restructured as standalone functions
  ([`genRandomCellTypeDistr()`](https://mhurtado13.github.io/racing/reference/genRandomCellTypeDistr.md),
  [`genRandomLigRecDistr()`](https://mhurtado13.github.io/racing/reference/genRandomLigRecDistr.md),
  [`genRandomCellLigands()`](https://mhurtado13.github.io/racing/reference/genRandomCellLigands.md),
  [`genRandomCellReceptors()`](https://mhurtado13.github.io/racing/reference/genRandomCellReceptors.md),
  [`genRandomCellTypeList()`](https://mhurtado13.github.io/racing/reference/genRandomCellTypeList.md),
  [`genRandomEdgeList()`](https://mhurtado13.github.io/racing/reference/genRandomEdgeList.md)),
  each independently documented and exported.
- **Consolidated `distribution_generation` into the same file**: The
  separate Python module `distribution_generation.py` (referenced as
  `import distribution_generation as dg`) is incorporated directly into
  `network_generation.R`, keeping all network generation logic in one
  place.
- **All functions are individually documented** with roxygen2 `@param`,
  `@return`, and `@export` tags.
- **[`model1()`](https://mhurtado13.github.io/racing/reference/model1.md)**:
  Converted from using a class instance to calling standalone functions.
  Parameters use R-style `NULL` defaults instead of Python empty lists
  `[]`.
- **Added
  [`generateUniformLRGraph()`](https://mhurtado13.github.io/racing/reference/generateUniformLRGraph.md)**:
  A new exported function that generates a graph under a uniformized
  (flat) LR distribution for a specific patient. This function does not
  exist in the original Python code and is useful for normalization
  workflows.
- **1-based indexing**: All vertex and cell-type indices use R’s 1-based
  convention vs. Python’s 0-based.

### 3.2 `feature_extraction.py` → `R/feature_extraction.R`

#### Original Python

Contains all graph motif counting functions: -
[`Find_Number_Trust_Triangles_Unique()`](https://mhurtado13.github.io/racing/reference/Find_Number_Trust_Triangles_Unique.md),
[`Find_Number_Triangles()`](https://mhurtado13.github.io/racing/reference/Find_Number_Triangles.md),
[`Find_Number_Triangles_Unique()`](https://mhurtado13.github.io/racing/reference/Find_Number_Triangles_Unique.md) -
[`Find_Number_2Loops()`](https://mhurtado13.github.io/racing/reference/Find_Number_2Loops.md),
[`Find_Number_2Loops_Unique()`](https://mhurtado13.github.io/racing/reference/Find_Number_2Loops_Unique.md) -
[`Find_Number_Wedges()`](https://mhurtado13.github.io/racing/reference/Find_Number_Wedges.md),
[`Find_Number_Wedges_Unique()`](https://mhurtado13.github.io/racing/reference/Find_Number_Wedges_Unique.md) -
[`Trust_Triangles()`](https://mhurtado13.github.io/racing/reference/Trust_Triangles.md),
[`Cycle_Triangles()`](https://mhurtado13.github.io/racing/reference/Cycle_Triangles.md),
[`Wedges()`](https://mhurtado13.github.io/racing/reference/Wedges.md) -
[`BFS()`](https://mhurtado13.github.io/racing/reference/BFS.md),
`StrongConnect()`, `Tarjan()`,
[`TarjanIterative()`](https://mhurtado13.github.io/racing/reference/TarjanIterative.md),
`sconnect()` -
[`GSCC()`](https://mhurtado13.github.io/racing/reference/GSCC.md),
[`IN()`](https://mhurtado13.github.io/racing/reference/IN.md),
[`OUT()`](https://mhurtado13.github.io/racing/reference/OUT.md)

Uses `scipy.sparse` matrix operations. Contains both recursive and
iterative Tarjan implementations.

#### R Package Changes

- **All functions faithfully ported** with identical algorithmic logic.
- **Dense matrix operations by default**: Where Python used
  `scipy.sparse` `.sign()`, `.getrow()`, `.nonzero()`, the R version
  uses standard matrix indexing (`Adj[v, ]`, `which(Adj[v, ] != 0)`).
  This is simpler but may be less memory-efficient for very large
  graphs. The adjacency matrix construction in `Utilities.R` does use
  [`Matrix::sparseMatrix`](https://rdrr.io/pkg/Matrix/man/sparseMatrix.html),
  but feature extraction functions operate on the result using standard
  R matrix operations.
- **Tarjan’s algorithm**: Only the iterative version (`TarjanIterative`)
  is ported. The recursive `Tarjan()` / `StrongConnect()` from Python is
  dropped, avoiding R’s recursion depth limitations for large graphs.
  The iterative version uses R closures (`<<-` assignment) to manage
  state, replacing Python’s mutable list/dict references.
- **Removed the recursive `Tarjan()`/`StrongConnect()` functions**
  entirely. The iterative-only approach is more robust for large graphs
  in R.
- **Added detailed inline comments** explaining the mathematics behind
  each motif counting operation.

### 3.3 `Kernel_Method.py` → `R/Kernel_Method.R`

#### Original Python

- `Calculate_kernel()`: Reads input files from disk (hardcoded naming:
  `{cancer}_TMEmod_cell_fractions.csv`), computes kernels, saves them as
  `.npz` files.
- `saveKernel()`: Saves kernel to disk as compressed NumPy arrays
  (`.npz`).
- [`calculateDirect()`](https://mhurtado13.github.io/racing/reference/calculateDirect.md),
  [`calculateWedges()`](https://mhurtado13.github.io/racing/reference/calculateWedges.md):
  Extract features from kernels and write CSV files to disk (naming
  convention: `{cancer}_{weight}_weight_Dir_bundle_norm.csv`).
- [`getGSCCAnalytically()`](https://mhurtado13.github.io/racing/reference/getGSCCAnalytically.md):
  Computes GSCC analytically using `nleqslv`-equivalent from SciPy,
  reads patient data from files.
- Heavy reliance on file I/O throughout; tightly coupled to a specific
  directory structure.

#### R Package Changes

- **[`compute_kernel()`](https://mhurtado13.github.io/racing/reference/compute_kernel.md)
  replaces `Calculate_kernel()`**: Accepts matrices as function
  arguments instead of reading files from disk. No file I/O inside the
  function—pure computation. Supports an optional `normalize` parameter
  to also return a uniformized baseline kernel in the same call.
- **[`calculateDirect()`](https://mhurtado13.github.io/racing/reference/calculateDirect.md)
  and
  [`calculateWedges()`](https://mhurtado13.github.io/racing/reference/calculateWedges.md)**:
  Rewritten to take kernel arrays and cell-name vectors as arguments and
  return data frames. No file writing. The `bundle` parameter for
  combining reciprocal directions is preserved.
- **[`computeGSCC()`](https://mhurtado13.github.io/racing/reference/computeGSCC.md)
  replaces
  [`getGSCCAnalytically()`](https://mhurtado13.github.io/racing/reference/getGSCCAnalytically.md)**:
  Same mathematical approach (Poisson branching process via `nleqslv`),
  but receives kernel and abundance matrices as arguments rather than
  reading from files. The
  [`poiBPFunc()`](https://mhurtado13.github.io/racing/reference/poiBPFunc.md)
  helper is preserved as an internal function.
- **[`getGSCCAnalytically()`](https://mhurtado13.github.io/racing/reference/getGSCCAnalytically.md)
  retained as legacy stub**: The original function name is kept but now
  [`stop()`](https://rdrr.io/r/base/stop.html)s with a message directing
  users to
  [`computeGSCC()`](https://mhurtado13.github.io/racing/reference/computeGSCC.md),
  ensuring backward-compatible error messages.
- **NEW:
  [`computeTriangles()`](https://mhurtado13.github.io/racing/reference/computeTriangles.md)**:
  Computes triangle features (trust triangles and bundled cycle
  triangles) directly from kernel arrays. This feature extraction
  function has no direct equivalent in the original Python
  `Kernel_Method.py`—the Python version only computed triangles via
  Monte Carlo simulation.
- **NEW:
  [`compute_kernel_features()`](https://mhurtado13.github.io/racing/reference/compute_kernel_features.md)**:
  A high-level dispatcher that routes to
  [`calculateDirect()`](https://mhurtado13.github.io/racing/reference/calculateDirect.md),
  [`calculateWedges()`](https://mhurtado13.github.io/racing/reference/calculateWedges.md),
  [`computeTriangles()`](https://mhurtado13.github.io/racing/reference/computeTriangles.md),
  or
  [`computeGSCC()`](https://mhurtado13.github.io/racing/reference/computeGSCC.md)
  based on a `communication_type` parameter. Supports patient subsetting
  via `patient_idx`. This function simplifies the user API.
- **NEW:
  [`compute_racing_kernel()`](https://mhurtado13.github.io/racing/reference/compute_racing_kernel.md)**:
  An end-to-end workflow function that:
  1.  Accepts either raw count matrices or pre-computed input data
  2.  Calls
      [`prepare_input_files()`](https://mhurtado13.github.io/racing/reference/prepare_input_files.md)
      (if needed)
  3.  Calls
      [`compute_kernel()`](https://mhurtado13.github.io/racing/reference/compute_kernel.md)
      with optional normalization
  4.  Calls
      [`compute_kernel_features()`](https://mhurtado13.github.io/racing/reference/compute_kernel_features.md)
      for the requested feature type
  5.  Returns kernel arrays and the feature data frame

  This function has no equivalent in the original Python code—users had
  to manually orchestrate each step.
- **In-memory kernel handling**: Kernels are returned as R arrays
  instead of being saved to disk as `.npz`, enabling interactive
  exploration.
- **Flexible file naming**: No assumptions about cancer type names,
  folder structures, or file naming conventions, allowing use with any
  dataset.

### 3.4 `Monte_Carlo_Method.py` → `R/Monte_Carlo_Method.R`

#### Original Python

- [`countWedges()`](https://mhurtado13.github.io/racing/reference/countWedges.md),
  [`countTrustTriangles()`](https://mhurtado13.github.io/racing/reference/countTrustTriangles.md),
  [`countCycleTriangles()`](https://mhurtado13.github.io/racing/reference/countCycleTriangles.md),
  [`countDirect()`](https://mhurtado13.github.io/racing/reference/countDirect.md),
  [`countGSCC()`](https://mhurtado13.github.io/racing/reference/countGSCC.md):
  Functions that run Monte Carlo simulations and count graph motifs.
- `runSimOne()`: Runs simulation for a single patient.
- [`runSim()`](https://mhurtado13.github.io/racing/reference/runSim.md):
  Iterates over patients, writes results to structured `.txt`/`.out`
  files with a well-defined layout.
- Output uses a purpose-built text format, with companion parsers in
  `txt_to_csv.py`.
- Provides console output for progress tracking.

#### R Package Changes

- **Core counting functions preserved**:
  [`countWedges()`](https://mhurtado13.github.io/racing/reference/countWedges.md),
  [`countTrustTriangles()`](https://mhurtado13.github.io/racing/reference/countTrustTriangles.md),
  [`countCycleTriangles()`](https://mhurtado13.github.io/racing/reference/countCycleTriangles.md),
  [`countDirect()`](https://mhurtado13.github.io/racing/reference/countDirect.md),
  [`countGSCC()`](https://mhurtado13.github.io/racing/reference/countGSCC.md)
  all faithfully replicate the Python logic with R syntax.
- **`runSimOne()` removed (commented out)**: The per-patient function is
  replaced by the loop inside
  [`runSim()`](https://mhurtado13.github.io/racing/reference/runSim.md)
  directly. The commented-out code is preserved for reference.
- **[`runSim()`](https://mhurtado13.github.io/racing/reference/runSim.md)
  improved**:
  - Added `output_folder` parameter with automatic directory creation.
  - Added `file.name` parameter for customizable output naming.
  - Added `patient_idx` parameter allowing simulation of a single
    specific patient without re-running all patients.
  - Uses R [`file()`](https://rdrr.io/r/base/connections.html)
    connections for writing (proper resource management via `con`).
  - Output format is structured CSV-like, compatible with
    [`Read_Sim_Output()`](https://mhurtado13.github.io/racing/reference/Read_Sim_Output.md).
- **NEW:
  [`compute_racing_montecarlo()`](https://mhurtado13.github.io/racing/reference/compute_racing_montecarlo.md)**:
  An end-to-end workflow function (parallel to
  [`compute_racing_kernel()`](https://mhurtado13.github.io/racing/reference/compute_racing_kernel.md))
  that:
  1.  Accepts raw counts or pre-computed input data
  2.  Optionally generates input files
  3.  Calls
      [`runSim()`](https://mhurtado13.github.io/racing/reference/runSim.md)
      for the requested feature type
  4.  Reads back results and returns as R data structures

  This end-to-end function has no equivalent in the original Python
  code.

### 3.5 `RaCInG_input_generation.py` → `R/RaCInG_input_generation.R`

#### Original Python

- [`sortPermute()`](https://mhurtado13.github.io/racing/reference/sortPermute.md):
  Sorts a list and returns the permutation.
- [`createCellLigList()`](https://mhurtado13.github.io/racing/reference/createCellLigList.md),
  [`createCellRecList()`](https://mhurtado13.github.io/racing/reference/createCellRecList.md),
  [`createCellTypeDistr()`](https://mhurtado13.github.io/racing/reference/createCellTypeDistr.md),
  [`createInteractionDistr()`](https://mhurtado13.github.io/racing/reference/createInteractionDistr.md):
  Read CSV files following a consistent naming convention (e.g.,
  `{cancer}_celltype_ligand_{weight}.csv`).
- [`Read_Lig_Rec_Interaction()`](https://mhurtado13.github.io/racing/reference/Read_Lig_Rec_Interaction.md):
  Reads a ligand-receptor interaction file.
- [`generateInput()`](https://mhurtado13.github.io/racing/reference/generateInput.md):
  Master function that calls all readers with hardcoded file naming.
- `get_patient_names()`: Extracts patient names from file names.

#### R Package Changes

- **All reader functions ported and generalized**:
  [`createCellLigList()`](https://mhurtado13.github.io/racing/reference/createCellLigList.md),
  [`createCellRecList()`](https://mhurtado13.github.io/racing/reference/createCellRecList.md),
  [`createCellTypeDistr()`](https://mhurtado13.github.io/racing/reference/createCellTypeDistr.md),
  [`createInteractionDistr()`](https://mhurtado13.github.io/racing/reference/createInteractionDistr.md)
  all accept explicit `filename` parameters rather than constructing
  paths from cancer type and weight type. This makes them usable with
  any naming convention.
- **[`createCellTypeDistr()`](https://mhurtado13.github.io/racing/reference/createCellTypeDistr.md)
  enhanced**: Includes M1/M2 macrophage merging logic that was done
  externally in the original Python workflow. The function now accepts a
  `cells` parameter and automatically merges M1+M2 if the number of cell
  types doesn’t match.
- **[`createInteractionDistr()`](https://mhurtado13.github.io/racing/reference/createInteractionDistr.md)
  enhanced**: Builds a proper 3D tensor (ligand × receptor × patient)
  instead of a flat matrix. Includes automatic normalization so each
  patient’s interaction distribution sums to 1. Handles NA values by
  replacing them with 0.
- **[`sortPermute()`](https://mhurtado13.github.io/racing/reference/sortPermute.md)
  preserved** as an internal helper.
- **[`generateInput()`](https://mhurtado13.github.io/racing/reference/generateInput.md)
  generalized**: Uses `file_name` and `output_folder` parameters instead
  of cancer type/weight type. Supports an optional `read_signs`
  parameter.
- **`get_patient_names()` removed**: Not needed because patient names
  are extracted directly from data frame row names in R.
- **NEW:
  [`prepare_input_files()`](https://mhurtado13.github.io/racing/reference/prepare_input_files.md)**:
  A comprehensive input preparation function that:
  1.  Normalizes counts to TPM (via
      [`ADImpute::NormalizeTPM`](https://rdrr.io/pkg/ADImpute/man/NormalizeTPM.html))
  2.  Runs deconvolution (via `multideconv`) or accepts pre-computed
      estimates
  3.  Estimates cell-type expression profiles
  4.  Retrieves ligand-receptor prior knowledge from OmniPath/LIANA (or
      accepts custom network)
  5.  Filters interactions by expression threshold (≥10 TPM)
  6.  Builds the cell-cell communication table
  7.  Computes all four input matrices (Lmatrix, Rmatrix, Cmatrix,
      LRmatrix)
  8.  Exports them as CSV files
  9.  Returns processed matrices ready for analysis

  This function replaces the three separate R Markdown notebooks
  (`RaCInG_ccc_prior_knowledge.Rmd`, `RaCInG_input_tcga.Rmd`,
  `RaCInG_input_published_cohorts.Rmd`) and the Python
  [`generateInput()`](https://mhurtado13.github.io/racing/reference/generateInput.md)
  function, consolidating over 500 lines of notebook code into a single
  callable function.

### 3.6 `Utilities.py` → `R/Utilities.R`

#### Original Python

- [`EdgetoAdj()`](https://mhurtado13.github.io/racing/reference/EdgetoAdj.md):
  Builds a sparse matrix from an edge list.
- [`EdgetoAdj_No_loop()`](https://mhurtado13.github.io/racing/reference/EdgetoAdj_No_loop.md):
  Same but removes self-loops.
- [`Count_Types()`](https://mhurtado13.github.io/racing/reference/Count_Types.md):
  Counts motifs by cell-type combinations.
- `createSlurm()`: Generates SLURM job submission scripts for HPC.
- [`poiBPFunc()`](https://mhurtado13.github.io/racing/reference/poiBPFunc.md):
  Poisson branching process helper for GSCC computation.

#### R Package Changes

- **[`EdgetoAdj()`](https://mhurtado13.github.io/racing/reference/EdgetoAdj.md)
  and
  [`EdgetoAdj_No_loop()`](https://mhurtado13.github.io/racing/reference/EdgetoAdj_No_loop.md)
  ported**: Use
  [`Matrix::sparseMatrix()`](https://rdrr.io/pkg/Matrix/man/sparseMatrix.html)
  instead of `scipy.sparse`. The functions are exported and documented.
- **[`Count_Types()`](https://mhurtado13.github.io/racing/reference/Count_Types.md)
  ported and generalized**: Now works with R’s 1-based indexing. Accepts
  an optional `maxTypes` parameter.
- **[`poiBPFunc()`](https://mhurtado13.github.io/racing/reference/poiBPFunc.md)
  preserved** as an internal helper.
- **`createSlurm()` not done yet**: It will be implemented in R using
  parallelization frameworks (`parallel`, `future`, `batchtools`) for
  users who need HPC capabilities.

### 3.7 `statistical_analysis.py` → `R/statistical_analysis.R`

#### Original Python

Comprehensive file (~800+ lines) providing a full downstream analysis
suite: - `remove_celltype()`, `readAllDataExact()`, `dataReadExact()`,
`data_read()`: Data reading/preprocessing functions tailored to the
original pipeline’s file structure. - `Find_Patient_Subtype_Bagaev()`:
Reads metadata from an Excel file (`Bagaev_mmc6.xlsx`). -
`Metadata_csv_read_in()`: Reads metadata CSV files. - `readAllData()`:
Reads Monte Carlo output CSVs. - `addMetadata()`, `removeOutliers()`:
Data cleaning functions. - `computeCorrelation()`: Spearman correlation
between features and immune response. - `groupPatients()`,
`createPanCancerData()`: Patient grouping for comparisons. -
`wilcoxon()`: Full Wilcoxon rank-sum test with Bonferroni correction,
pan-cancer analysis. - `volcanoPan()`, `volcanoInd()`, `volcanoCross()`:
Matplotlib-based volcano plots (~300 lines). - `findLargeFold()`,
`plot_heatmap()`, [`heatmap()`](https://rdrr.io/r/stats/heatmap.html):
Fold-change analysis and heatmaps. - `contributionAnalysisGSCC()`: GSCC
contribution analysis with Mann-Whitney tests.

#### R Package Changes

- **Generalized approach**: The original’s extensive analysis suite is
  distilled into two flexible, reusable functions (TO BE COMPLETED WITH
  ADDITIONAL FUNCTIONS OF demo.py):
  1.  **[`wilcox_group_test()`](https://mhurtado13.github.io/racing/reference/wilcox_group_test.md)**:
      A generic Wilcoxon rank-sum test function that accepts any data
      matrix and group labels. Uses
      [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html) for
      multiple testing correction (supports FDR, Bonferroni, and others
      via a parameter). Returns a tidy data frame sorted by adjusted
      p-value. This generalizes the original’s cancer-specific testing
      pipeline into a dataset-agnostic tool.
  2.  **[`volcano_plot()`](https://mhurtado13.github.io/racing/reference/volcano_plot.md)**:
      A `ggplot2`-based volcano plot that accepts the output of
      [`wilcox_group_test()`](https://mhurtado13.github.io/racing/reference/wilcox_group_test.md).
- **Functions not carried over** (with rationale):
  - Data reading functions (`readAllDataExact`, `dataReadExact`,
    `data_read`, `readAllData`): Superseded by the generic input
    pipeline.
  - `Find_Patient_Subtype_Bagaev()`: Designed for a specific metadata
    source (`Bagaev_mmc6.xlsx`); users can implement study-specific
    metadata loaders as needed.
  - `Metadata_csv_read_in()`: Standard CSV reading is natively supported
    in R via [`read.csv()`](https://rdrr.io/r/utils/read.table.html).
  - `addMetadata()`, `removeOutliers()`: Data cleaning operations
    well-served by R’s existing ecosystem (`dplyr`, `tidyr`).
  - `groupPatients()`, `createPanCancerData()`: Study-specific grouping
    logic that users can implement based on their particular study
    design.
  - `computeCorrelation()`: Correlation analysis is well-supported
    natively in R via
    [`cor.test()`](https://rdrr.io/r/stats/cor.test.html) and related
    functions.
  - `findLargeFold()`, `plot_heatmap()`,
    [`heatmap()`](https://rdrr.io/r/stats/heatmap.html): Mature R
    packages like `ComplexHeatmap` and `pheatmap` offer extensive
    heatmap capabilities.
  - `contributionAnalysisGSCC()`: Study-specific analysis that users can
    implement using
    [`computeGSCC()`](https://mhurtado13.github.io/racing/reference/computeGSCC.md)
    output.
  - `volcanoCross()`, `volcanoPan()`, `volcanoInd()`: Consolidated into
    the single generic
    [`volcano_plot()`](https://mhurtado13.github.io/racing/reference/volcano_plot.md).

### 3.8 `txt_to_csv.py` → `R/txt_to_csv.R`

#### Original Python

- `Triangle_Prop_Read()`: Parses custom `.txt` output files for triangle
  data.
- `Direct_Comm_Read()`: Parses direct communication `.txt` files.
- `GSCC_Read()`: Parses GSCC `.txt` files.
- `Generate_normalised_count_csv()`: Reads raw simulations and produces
  normalized CSV files. Includes direction bundling logic.

#### R Package Changes

- **[`Read_Sim_Output()`](https://mhurtado13.github.io/racing/reference/Read_Sim_Output.md)
  unifies the three readers**: A single function that auto-detects the
  communication type (D, W, TT, CT, GSCC) from the file header and
  adapts its parsing logic accordingly. The original Python version used
  three dedicated parsers, each optimized for its respective data
  dimensionality (1D for GSCC, 2D for Direct, 3D for triangles/wedges).
- **`Generate_normalised_count_csv()` included as
  [`compute_results_processing()`](https://mhurtado13.github.io/racing/reference/compute_results_processing.md)**:
  Handles normalization of Monte Carlo output, direction bundling, and
  CSV export.
- **Output format adapted**: The R package writes a labeled CSV-like
  format from
  [`runSim()`](https://mhurtado13.github.io/racing/reference/runSim.md),
  building on the original’s structured text format with added column
  headers for self-describing output.

### 3.9 Files in the Original with No Direct R Package Equivalent

| Original File | Purpose | R Package Status |
|----|----|----|
| `Circos.py` | Circos plot generation for visualizing cell-cell interactions | **Not yet ported**. Planned for a future release. |
| `HPC_CLI.py` | Command-line interface for HPC cluster submission | **Not yet ported**. The R package currently focuses on interactive use. Planned for a future release. |
| `R_Code/RaCInG_ccc_prior_knowledge.Rmd` | Notebook for building prior knowledge from OmniPath/Ramilowski/CCLE | **Integrated into [`prepare_input_files()`](https://mhurtado13.github.io/racing/reference/prepare_input_files.md)** as automated pipeline steps. |
| `R_Code/RaCInG_input_tcga.Rmd` | Notebook for processing TCGA data | **Integrated into [`prepare_input_files()`](https://mhurtado13.github.io/racing/reference/prepare_input_files.md)**. |
| `R_Code/RaCInG_input_published_cohorts.Rmd` | Notebook for processing published cohort data | **Integrated into [`prepare_input_files()`](https://mhurtado13.github.io/racing/reference/prepare_input_files.md)**. |
| `R_Code/utils/run_TMEmod_deconvolution.R` | Script to run 6 deconvolution methods | **Replaced** by integration with `multideconv` package. |
| `Python_Code/Demo.ipynb` | Jupyter notebook demo | **Replaced** by the R vignette (`vignettes/RaCiNG.Rmd`) and pkgdown website. |

------------------------------------------------------------------------

## 4. New Features and Additions

### 4.1 Bundled Example Dataset (`skcm_example`)

The R package includes a pre-processed SKCM (skin cutaneous melanoma)
dataset as built-in package data (`data/skcm_example.rda`). This allows
users to immediately run RaCInG without downloading or preparing
external data:

``` r

data(skcm_example)
result <- compute_racing_kernel(input_data = skcm_example)
```

The original Python code provided example input files in
`Example input/` but they were not part of a formal data distribution
system.

### 4.2 End-to-End Workflow Functions

Two new high-level entry points that did not exist in the original:

- **[`compute_racing_kernel()`](https://mhurtado13.github.io/racing/reference/compute_racing_kernel.md)**:
  From raw counts → input preparation → kernel computation → feature
  extraction in one call.
- **[`compute_racing_montecarlo()`](https://mhurtado13.github.io/racing/reference/compute_racing_montecarlo.md)**:
  From raw counts → input preparation → Monte Carlo simulation → result
  processing in one call.

These functions support both “from scratch” workflows (providing a count
matrix) and “bring your own data” workflows (providing the `input_data`
list).

### 4.3 `compute_kernel_features()` Dispatcher

A new function that acts as a single entry point for computing any type
of kernel-derived feature (Direct, Wedges, Trust Triangles, GSCC), with
patient subsetting support.

### 4.4 Kernel-Based Triangle Computation

The function
[`computeTriangles()`](https://mhurtado13.github.io/racing/reference/computeTriangles.md)
computes triangle features analytically from the kernel, which was only
available via Monte Carlo simulation in the original Python code. This
eliminates the need for computationally expensive graph simulations when
only triangle features are needed.

### 4.5 `generateUniformLRGraph()`

A new helper function for generating graphs under a uniform
ligand-receptor distribution. This is useful for normalization studies
and understanding the null model.

### 4.6 pkgdown Documentation Website

A full documentation website built with pkgdown (Bootstrap 5), available
at <https://mhurtado13.github.io/racing/>, including: - Getting Started
vignette with step-by-step walkthrough - Input files explanation section

### 4.7 Dependency Validation

The
[`prepare_input_files()`](https://mhurtado13.github.io/racing/reference/prepare_input_files.md)
function uses
[`.check_installed_packages()`](https://mhurtado13.github.io/racing/reference/dot-check_installed_packages.md)
to verify that optional dependencies (`ADImpute`, `multideconv`,
`OmnipathR`, `liana`) are installed before attempting to use them,
providing clear error messages instead of cryptic import failures.

------------------------------------------------------------------------

## 5. Features Not Included in the R Package

| Feature | Original Implementation | Rationale |
|----|----|----|
| SLURM job generation | `createSlurm()` in `Utilities.py` | R provides its own parallelization ecosystem (`parallel`, `future`, `batchtools`) |
| Circos plots | `Circos.py` | Planned for a future release; meanwhile the `circlize` R package is available |
| HPC CLI | `HPC_CLI.py` | Planned for a future release |
| Bagaev metadata loader | `Find_Patient_Subtype_Bagaev()` | Designed for a specific dataset; users can implement study-specific loaders |
| Pan-cancer pipeline | `groupPatients()`, `createPanCancerData()`, cross-cancer volcano plots | Study-specific analysis logic; users implement based on their study design |
| Heatmap pipeline | `findLargeFold()`, `plot_heatmap()`, [`heatmap()`](https://rdrr.io/r/stats/heatmap.html) | Mature R packages (`ComplexHeatmap`, `pheatmap`) provide extensive heatmap support |
| Correlation analysis | `computeCorrelation()` | Well-supported natively in R via [`cor.test()`](https://rdrr.io/r/stats/cor.test.html) and related functions |
| Recursive Tarjan | `Tarjan()`, `StrongConnect()` | Iterative version preferred to accommodate R’s default recursion limits for large graphs |

------------------------------------------------------------------------

## 6. API Design Improvements

### 6.1 Pure Functional Design

The original Python code integrated computation with file I/O, which was
well-suited for its batch processing workflow: - `Calculate_kernel()`
reads CSV files, computes kernels, and writes `.npz` files. -
[`calculateDirect()`](https://mhurtado13.github.io/racing/reference/calculateDirect.md)
computes features and writes CSV files. -
[`runSim()`](https://mhurtado13.github.io/racing/reference/runSim.md)
computes features and writes `.txt` files.

The R package separates concerns: - **Computation functions**
([`compute_kernel()`](https://mhurtado13.github.io/racing/reference/compute_kernel.md),
[`calculateDirect()`](https://mhurtado13.github.io/racing/reference/calculateDirect.md),
etc.) are pure: they accept data as arguments and return results. -
**I/O functions**
([`prepare_input_files()`](https://mhurtado13.github.io/racing/reference/prepare_input_files.md),
[`Read_Sim_Output()`](https://mhurtado13.github.io/racing/reference/Read_Sim_Output.md))
handle file operations explicitly. - **Workflow functions**
([`compute_racing_kernel()`](https://mhurtado13.github.io/racing/reference/compute_racing_kernel.md),
[`compute_racing_montecarlo()`](https://mhurtado13.github.io/racing/reference/compute_racing_montecarlo.md))
compose computation and I/O.

### 6.2 Consistent Parameter Naming

The original Python code uses concise parameter names suited to its
script-based workflow: - `cellLigList` vs `structurelig`, `Dcelltype` vs
`Dligrec` - `N` (cells per graph), `itNo` (iterations), `av` (average
degree)

The R package maintains backward compatibility with the original
parameter names where possible while adding descriptive roxygen2
documentation for each parameter.

### 6.3 Flexible Input Modes

The
[`compute_racing_kernel()`](https://mhurtado13.github.io/racing/reference/compute_racing_kernel.md)
and
[`compute_racing_montecarlo()`](https://mhurtado13.github.io/racing/reference/compute_racing_montecarlo.md)
functions accept either: 1. Raw count matrices (triggering the full
preprocessing pipeline) 2. Pre-computed input data lists (bypassing
preprocessing)

This dual-mode design accommodates both new users and those with
existing preprocessed data.

------------------------------------------------------------------------

## 7. Data Flow Comparison

### Original Python Workflow

    CSV files on disk
      → generateInput() [reads files with hardcoded names]
      → model1() [generates graph]
      → feature extraction functions
      → runSim() [writes .txt files]
      → txt_to_csv functions [reads .txt, writes .csv]
      → statistical_analysis functions [reads .csv, writes plots]

Each step uses disk I/O, providing a structured and reproducible
file-based pipeline.

### R Package Workflow (Kernel)

    Gene count matrix (or pre-computed input_data list)
      → compute_racing_kernel()
        → prepare_input_files() [optional; one-time I/O]
        → compute_kernel() [in-memory]
        → compute_kernel_features() [in-memory]
      → returns list(kernel, kernel_norm, features)

Minimal disk I/O. Data stays in memory throughout the analysis.

### R Package Workflow (Monte Carlo)

    Gene count matrix (or pre-computed input_data list)
      → compute_racing_montecarlo()
        → prepare_input_files() [optional]
        → runSim() [writes .out file]
        → Read_Sim_Output() + compute_results_processing()
      → returns processed results

File I/O only for Monte Carlo output persistence (necessary due to
simulation size).

------------------------------------------------------------------------

## 8. Summary of All Exported Functions

### Functions Ported from Python (with same or similar name):

| R Function | Python Origin | Changes |
|----|----|----|
| [`sortPermute()`](https://mhurtado13.github.io/racing/reference/sortPermute.md) | [`sortPermute()`](https://mhurtado13.github.io/racing/reference/sortPermute.md) | Minimal; R syntax |
| [`createCellLigList()`](https://mhurtado13.github.io/racing/reference/createCellLigList.md) | [`createCellLigList()`](https://mhurtado13.github.io/racing/reference/createCellLigList.md) | Generic filename; adds sorting |
| [`createCellRecList()`](https://mhurtado13.github.io/racing/reference/createCellRecList.md) | [`createCellRecList()`](https://mhurtado13.github.io/racing/reference/createCellRecList.md) | Generic filename; adds sorting |
| [`createCellTypeDistr()`](https://mhurtado13.github.io/racing/reference/createCellTypeDistr.md) | [`createCellTypeDistr()`](https://mhurtado13.github.io/racing/reference/createCellTypeDistr.md) | Adds M1/M2 merge logic |
| [`createInteractionDistr()`](https://mhurtado13.github.io/racing/reference/createInteractionDistr.md) | [`createInteractionDistr()`](https://mhurtado13.github.io/racing/reference/createInteractionDistr.md) | Returns 3D tensor; auto-normalizes |
| [`Read_Lig_Rec_Interaction()`](https://mhurtado13.github.io/racing/reference/Read_Lig_Rec_Interaction.md) | [`Read_Lig_Rec_Interaction()`](https://mhurtado13.github.io/racing/reference/Read_Lig_Rec_Interaction.md) | Uses [`data.table::fread`](https://rdrr.io/pkg/data.table/man/fread.html) |
| [`generateInput()`](https://mhurtado13.github.io/racing/reference/generateInput.md) | [`generateInput()`](https://mhurtado13.github.io/racing/reference/generateInput.md) | Generalized file naming |
| [`model1()`](https://mhurtado13.github.io/racing/reference/model1.md) | [`model1()`](https://mhurtado13.github.io/racing/reference/model1.md) | Standalone functions instead of class |
| [`genRandomCellTypeDistr()`](https://mhurtado13.github.io/racing/reference/genRandomCellTypeDistr.md) | `dg.genRandomCellTypeDistr()` | Direct port |
| [`genRandomLigRecDistr()`](https://mhurtado13.github.io/racing/reference/genRandomLigRecDistr.md) | `dg.genRandomLigRecDistr()` | Direct port |
| [`genRandomCellLigands()`](https://mhurtado13.github.io/racing/reference/genRandomCellLigands.md) | `dg.genRandomCellLigands()` | Direct port |
| [`genRandomCellReceptors()`](https://mhurtado13.github.io/racing/reference/genRandomCellReceptors.md) | `dg.genRandomCellReceptors()` | Direct port |
| [`genRandomCellTypeList()`](https://mhurtado13.github.io/racing/reference/genRandomCellTypeList.md) | `informationMicroEnv.genRandomCellTypeList()` | Extracted from class |
| [`genRandomEdgeList()`](https://mhurtado13.github.io/racing/reference/genRandomEdgeList.md) | `informationMicroEnv.genRandomEdgeList()` | Extracted from class |
| [`EdgetoAdj()`](https://mhurtado13.github.io/racing/reference/EdgetoAdj.md) | [`EdgetoAdj()`](https://mhurtado13.github.io/racing/reference/EdgetoAdj.md) | Uses [`Matrix::sparseMatrix`](https://rdrr.io/pkg/Matrix/man/sparseMatrix.html) |
| [`EdgetoAdj_No_loop()`](https://mhurtado13.github.io/racing/reference/EdgetoAdj_No_loop.md) | [`EdgetoAdj_No_loop()`](https://mhurtado13.github.io/racing/reference/EdgetoAdj_No_loop.md) | Uses [`Matrix::sparseMatrix`](https://rdrr.io/pkg/Matrix/man/sparseMatrix.html) |
| [`Count_Types()`](https://mhurtado13.github.io/racing/reference/Count_Types.md) | [`Count_Types()`](https://mhurtado13.github.io/racing/reference/Count_Types.md) | 1-based indexing |
| [`poiBPFunc()`](https://mhurtado13.github.io/racing/reference/poiBPFunc.md) | [`poiBPFunc()`](https://mhurtado13.github.io/racing/reference/poiBPFunc.md) | Internal; same math |
| [`Find_Number_Trust_Triangles_Unique()`](https://mhurtado13.github.io/racing/reference/Find_Number_Trust_Triangles_Unique.md) | Same | Dense matrix operations |
| [`Find_Number_Triangles()`](https://mhurtado13.github.io/racing/reference/Find_Number_Triangles.md) | Same | Direct port |
| [`Find_Number_Triangles_Unique()`](https://mhurtado13.github.io/racing/reference/Find_Number_Triangles_Unique.md) | Same | Direct port |
| [`Find_Number_2Loops()`](https://mhurtado13.github.io/racing/reference/Find_Number_2Loops.md) | Same | Direct port |
| [`Find_Number_2Loops_Unique()`](https://mhurtado13.github.io/racing/reference/Find_Number_2Loops_Unique.md) | Same | Direct port |
| [`Find_Number_Wedges()`](https://mhurtado13.github.io/racing/reference/Find_Number_Wedges.md) | Same | Direct port |
| [`Find_Number_Wedges_Unique()`](https://mhurtado13.github.io/racing/reference/Find_Number_Wedges_Unique.md) | Same | Direct port |
| [`Trust_Triangles()`](https://mhurtado13.github.io/racing/reference/Trust_Triangles.md) | Same | Returns list instead of tuple |
| [`Cycle_Triangles()`](https://mhurtado13.github.io/racing/reference/Cycle_Triangles.md) | Same | Returns list instead of tuple |
| [`Wedges()`](https://mhurtado13.github.io/racing/reference/Wedges.md) | Same | Returns list instead of tuple |
| [`BFS()`](https://mhurtado13.github.io/racing/reference/BFS.md) | Same | Internal; 1-based indexing |
| [`TarjanIterative()`](https://mhurtado13.github.io/racing/reference/TarjanIterative.md) | Same | Uses R closures |
| [`GSCC()`](https://mhurtado13.github.io/racing/reference/GSCC.md) | Same | Direct port |
| [`IN()`](https://mhurtado13.github.io/racing/reference/IN.md) | Same | Direct port |
| [`OUT()`](https://mhurtado13.github.io/racing/reference/OUT.md) | Same | Direct port |
| [`countWedges()`](https://mhurtado13.github.io/racing/reference/countWedges.md) | Same | Direct port |
| [`countTrustTriangles()`](https://mhurtado13.github.io/racing/reference/countTrustTriangles.md) | Same | Direct port |
| [`countCycleTriangles()`](https://mhurtado13.github.io/racing/reference/countCycleTriangles.md) | Same | Direct port |
| [`countDirect()`](https://mhurtado13.github.io/racing/reference/countDirect.md) | Same | Direct port |
| [`countGSCC()`](https://mhurtado13.github.io/racing/reference/countGSCC.md) | Same | Direct port |
| [`runSim()`](https://mhurtado13.github.io/racing/reference/runSim.md) | Same | Enhanced with output_folder, file.name, patient_idx |
| [`calculateDirect()`](https://mhurtado13.github.io/racing/reference/calculateDirect.md) | Same | No file I/O; returns data.frame |
| [`calculateWedges()`](https://mhurtado13.github.io/racing/reference/calculateWedges.md) | Same | No file I/O; returns data.frame |
| [`getGSCCAnalytically()`](https://mhurtado13.github.io/racing/reference/getGSCCAnalytically.md) | Same | Legacy stub; redirects to computeGSCC() |
| [`Read_Sim_Output()`](https://mhurtado13.github.io/racing/reference/Read_Sim_Output.md) | `Triangle_Prop_Read()` + `Direct_Comm_Read()` + `GSCC_Read()` | Unified auto-detecting parser |
| [`wilcox_group_test()`](https://mhurtado13.github.io/racing/reference/wilcox_group_test.md) | `wilcoxon()` | Simplified, generic, FDR-corrected |
| [`volcano_plot()`](https://mhurtado13.github.io/racing/reference/volcano_plot.md) | `volcanoPan()` + `volcanoInd()` + `volcanoCross()` | Unified ggplot2 version |

### New Functions (no Python equivalent):

| R Function | Purpose |
|----|----|
| [`compute_kernel()`](https://mhurtado13.github.io/racing/reference/compute_kernel.md) | Pure kernel computation from matrices |
| [`computeGSCC()`](https://mhurtado13.github.io/racing/reference/computeGSCC.md) | Analytical GSCC from kernel (function-argument based) |
| [`computeTriangles()`](https://mhurtado13.github.io/racing/reference/computeTriangles.md) | Kernel-based triangle feature computation |
| [`compute_kernel_features()`](https://mhurtado13.github.io/racing/reference/compute_kernel_features.md) | Feature dispatcher for all kernel feature types |
| [`compute_racing_kernel()`](https://mhurtado13.github.io/racing/reference/compute_racing_kernel.md) | End-to-end kernel workflow |
| [`compute_racing_montecarlo()`](https://mhurtado13.github.io/racing/reference/compute_racing_montecarlo.md) | End-to-end Monte Carlo workflow |
| [`compute_results_processing()`](https://mhurtado13.github.io/racing/reference/compute_results_processing.md) | Normalize and bundle Monte Carlo results |
| [`prepare_input_files()`](https://mhurtado13.github.io/racing/reference/prepare_input_files.md) | Automated input preparation from raw counts |
| [`generateUniformLRGraph()`](https://mhurtado13.github.io/racing/reference/generateUniformLRGraph.md) | Generate graph under uniform LR null model |
| [`.check_installed_packages()`](https://mhurtado13.github.io/racing/reference/dot-check_installed_packages.md) | Internal dependency checker |

------------------------------------------------------------------------

## 9. Conclusion

The R package builds on the solid mathematical foundation established by
the original Python implementation, adapting it for the R ecosystem:

1.  **Packaging**: Organized as a formal, installable R package with
    automated documentation and a pkgdown website.
2.  **Usability**: Introduced end-to-end workflow functions that
    streamline multi-step analyses into single calls.
3.  **Generality**: Parameterized file naming and input handling, making
    the package applicable to any dataset without modification.
4.  **Separation of concerns**: Separated computation from I/O, enabling
    both interactive exploration and batch processing.
5.  **New capabilities**: Kernel-based triangle computation, unified
    feature dispatcher, bundled example data, and end-to-end workflow
    functions.
6.  **Focused scope**: Downstream statistical analysis distilled into
    two generic, reusable functions, leveraging R’s rich existing
    ecosystem for specialized tasks (to be completed with additional
    functions).

The core mathematics (random graph generation, kernel computation, motif
counting, branching-process GSCC) is preserved faithfully from the
original, while the software engineering has been adapted for R-native
workflows and package distribution.
