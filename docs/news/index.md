# Changelog

## RaCInG 0.3.0

### Removed

- Removed 11 functions confirmed unused anywhere in the package, its
  tests, vignette, or README, and unreferenced by any downstream
  project: the count-only `Find_Number_Trust_Triangles_Unique()`,
  `Find_Number_Triangles()`, `Find_Number_Triangles_Unique()`,
  `Find_Number_2Loops()`, `Find_Number_2Loops_Unique()`,
  `Find_Number_Wedges()`, `Find_Number_Wedges_Unique()` (superseded by
  [`Trust_Triangles()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/Trust_Triangles.md)/
  [`Cycle_Triangles()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/Cycle_Triangles.md)/[`Wedges()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/Wedges.md),
  which already return the count plus the full motif list); the bow-tie
  decomposition `IN()`/`OUT()` and their sole dependency `BFS()`;
  `generateUniformLRGraph()` (superseded by the `norm`/ `normalize`
  logic built directly into
  [`runSim()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/runSim.md)/[`compute_kernel()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/compute_kernel.md));
  and
  [`model1()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/model1.md)’s
  `genRandom = TRUE` demo path along with the 4 fabricated-data
  generators it called (`genRandomCellTypeDistr()`,
  `genRandomLigRecDistr()`, `genRandomCellLigands()`,
  `genRandomCellReceptors()`) –
  [`model1()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/model1.md)
  now requires real `cellLigList`/`cellRecList`/`Dcelltype`/`Dligrec`
  arguments. Several of these were direct ports of code the original
  Python implementation’s own authors had already set aside as test-only
  scaffolding (`distribution_generation.py`, explicitly documented there
  as “not used in paper”).
- Removed `volcano_plot()` (labels were dropped for busy plots and it
  couldn’t scale to features with long names); replaced by
  [`top_features_plot()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/top_features_plot.md),
  a ranked horizontal bar chart of the top up/down features by
  `log2(fold change)`.

### Multiple feature types without recomputing

- [`compute_racing_kernel()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/compute_racing_kernel.md)/[`compute_kernel_features()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/compute_kernel_features.md)
  now accept `communication_type` as a character vector (e.g.
  `c("D", "W", "TT", "GSCC")`). The kernel is still computed exactly
  once per call regardless of how many types are requested; `features`
  becomes a named list (one data frame per type) when more than one is
  given.
- [`compute_racing_montecarlo()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/compute_racing_montecarlo.md)/[`runSim()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/runSim.md)
  now accept `communication_type` as a character vector too (any
  combination of `"D"`, `"W"`, `"TT"`, `"CT"`, `"GSCC"`), via a new
  [`countAllTypes()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/countAllTypes.md)
  that extracts every requested type from the *same* simulated graphs
  instead of resimulating a fresh set of graphs per type – graph
  generation, not feature extraction, is the expensive part of a Monte
  Carlo run. `output` becomes a named list (one entry per type, one
  `.out`/CSV file pair per type on disk) when more than one is given.
  Single-type usage is unchanged (verified bit-identical against the
  previous single-type implementation with a fixed seed).

### Fixed a catastrophic O(n^3) crash/hang in kernel Wedge/Triangle features

- [`calculateWedges()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/calculateWedges.md)
  and
  [`computeTriangles()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/computeTriangles.md)
  built one named-list entry per cell-type combination inside a
  triple-nested loop, then converted the whole list to a data frame –
  for realistic cell-type counts (e.g. 100+) this reaches 10⁵⁻¹⁰6+
  individual list insertions and made `"W"`/`"TT"` effectively unusable
  via the kernel method, unlike `"D"`. Rewritten as vectorized array
  operations (index grid + flattened-matrix lookups, mirroring the
  equivalent Monte Carlo fix below); verified numerically identical to
  the original loop-based output across bundled/unbundled and
  raw/normalized cases. At 106 real cell types, `"TT"` (previously
  unusable) now completes in a few seconds.
- [`calculateDirect()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/calculateDirect.md),
  [`calculateWedges()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/calculateWedges.md),
  and
  [`computeTriangles()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/computeTriangles.md)
  now drop cell-type combinations with no possible ligand-receptor
  pathway in any patient (structurally zero for every patient) instead
  of returning them as dead all-zero columns – for
  [`calculateWedges()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/calculateWedges.md)/[`computeTriangles()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/computeTriangles.md)
  this filtering happens *before* the abundance-weighting/ratio
  computation, so it also reduces the actual memory/compute cost, not
  just the output size.

### Monte Carlo memory and parallelization fixes

- [`runSim()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/runSim.md)’s
  `ncores > 1` path
  ([`parallel::makeCluster()`](https://rdrr.io/r/parallel/makeCluster.html) +
  [`doParallel::registerDoParallel()`](https://rdrr.io/pkg/doParallel/man/registerDoParallel.html) +
  `foreach`, replacing the previous fork-based
  [`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html) for
  Windows portability) was shipping a full copy of the entire
  patient-by-ligand-by-receptor `LRmatrix` tensor to *every* worker
  process, regardless of patient count – for a 106 cell-type,
  116-patient cohort this meant ~460MB duplicated per worker even for a
  trivial single-patient run. Fixed by pre-slicing `LRmatrix`/`Cmatrix`
  per patient and making those small slices the `foreach` iteration
  variable (dispatched per task) instead of a name `foreach` broadcasts
  whole to every worker; confirmed via monitored testing this dropped
  per-worker overhead from ~460MB to a few MB.
- [`compute_results_processing()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/compute_results_processing.md)’s
  `remove_direction` column-merging grouped and relabeled columns by
  splitting the pasted column-name string on `"_"`, which silently
  scrambles cell-type identity whenever a cell-type name itself contains
  `"_"` (e.g. multideconv’s method-tagged names like
  `"CBSX_CBSX.HNSCC.scRNAseq_B.cells"`) – there is no way to tell, from
  the string alone, which underscores separate the two cell types in a
  pair from underscores inside one cell type’s own name. Fixed by
  grouping on the underlying integer index grid instead of re-parsing
  labels; added a regression test.
- [`calculateDirect()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/calculateDirect.md)/[`calculateWedges()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/calculateWedges.md)
  silently dropped real patient names (always fell back to
  `"Patient_1"`, `"Patient_2"`, …) even when
  [`compute_kernel_features()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/compute_kernel_features.md)/[`compute_racing_kernel()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/compute_racing_kernel.md)
  had real ones available – `"TT"`/`"GSCC"` already used them correctly.
  Fixed, and kernel and Monte Carlo feature matrices now use identical
  column-naming conventions (`"Dir_A_B"`, `"W_A_B_C"`, etc.) too, so
  their outputs can be fed into the same downstream statistical/plotting
  functions interchangeably.

### Other changes

- Added `ncores` to
  [`runSim()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/runSim.md)/[`compute_racing_montecarlo()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/compute_racing_montecarlo.md)
  for patient-level parallelization.
- Added `already_normalized` and `expr_threshold` parameters to
  [`prepare_input_files()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/prepare_input_files.md).
- Vectorized several remaining O(n^(2)-O(n)3) hotspots:
  [`prepare_input_files()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/prepare_input_files.md)’s
  `ccc_table` construction,
  [`Read_Sim_Output()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/Read_Sim_Output.md)’s
  `.out`-file parsing, and
  [`runSim()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/runSim.md)’s
  composition-block writing.
- [`createCellTypeDistr()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/createCellTypeDistr.md)
  now drops `Cmatrix` columns not present in `Lmatrix`/`Rmatrix` and
  renormalizes, instead of leaving mismatched columns in place.
- [`calculateDirect()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/calculateDirect.md)/[`calculateWedges()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/calculateWedges.md)
  fill structurally-impossible (0/0) normalized ratios with `0`.
- Expanded the `testthat` regression suite to cover all of the above.

## RaCInG 0.2.0

### Critical fixes

- **[`model1()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/model1.md)
  /
  [`genRandomEdgeList()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/genRandomEdgeList.md)
  ignored the requested average degree.**
  [`genRandomEdgeList()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/genRandomEdgeList.md)
  never received the caller’s edge count and instead always sampled
  `nrow(Dligrec) * ncol(Dligrec)` edges, regardless of `N` or
  `av`/`avdeg`. Every Monte Carlo simulation run against a real
  ligand-receptor network generated graphs with the wrong density (often
  drastically so) in either direction, unrelated to the requested
  average degree. Simulation results computed with earlier versions
  should be treated as invalid and re-run.
- **[`TarjanIterative()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/TarjanIterative.md)
  ([`GSCC()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/GSCC.md)/`IN()`/`OUT()`/[`countGSCC()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/countGSCC.md))
  was broken.** A scoping bug (`stack <- ...` instead of `<<-` inside a
  closure) desynced the SCC stack from tracked state, and
  `seq(i + 1, length(neigh))` produced a descending sequence for any
  zero-out-degree vertex - both caused crashes or silently wrong
  strongly-connected-component results on ordinary graphs. Rewritten as
  a standard iterative Tarjan and cross-validated against a brute-force
  reachability oracle on 200 random graphs.
- **Kernel-based Wedge/Triangle/GSCC features double-counted shared
  nodes.**
  [`compute_kernel()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/compute_kernel.md)
  baked cell-type abundance into both kernel indices, which is correct
  for Direct communication but double-weights any node shared by more
  than one edge in a motif (the middle node of a wedge, all three nodes
  of a triangle, and - since the underlying branching-process equations
  are nonlinear - every GSCC computation regardless of normalization).
  [`compute_kernel()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/compute_kernel.md)
  now returns the unweighted structural kernel, and each feature
  function multiplies in abundance explicitly, once per node, matching
  the original implementation.
- **Triangle features were mislabeled.** The unbundled “trust triangle”
  (TT) computation actually computed the cycle-triangle (CT) formula,
  and the true trust-triangle metric was missing entirely; the bundled
  8-term formula dropped 6 of its 8 terms.
  [`computeTriangles()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/computeTriangles.md)
  now returns both `TT_*` and `CT_*` columns when unbundled, and the
  correct 8-term sum when bundled.
- Self-loops were not excluded before Monte Carlo motif counting
  ([`EdgetoAdj()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/EdgetoAdj.md)
  was used where
  [`EdgetoAdj_No_loop()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/EdgetoAdj_No_loop.md)
  was required), inflating direct/wedge/triangle counts.
- `Find_Number_Trust_Triangles_Unique()` compared against the raw wedge
  count instead of its sign, undercounting whenever two
  directly-connected nodes had two or more common neighbors.
- Cell-type identifiers (patient row names) were silently dropped from
  `Cmatrix` in
  [`createCellTypeDistr()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/createCellTypeDistr.md),
  breaking `patient_names` in downstream outputs.
- [`Read_Sim_Output()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/Read_Sim_Output.md)
  inferred the number of cell types via
  [`ncol()`](https://rdrr.io/r/base/nrow.html) on a ragged, NA-padded
  data frame, which broke for datasets with fewer than 4 cell types.
- Empty motif lists
  ([`Wedges()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/Wedges.md)/[`Trust_Triangles()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/Trust_Triangles.md)/[`Cycle_Triangles()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/Cycle_Triangles.md)
  returning zero rows) crashed
  [`Count_Types()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/Count_Types.md).
- [`Cycle_Triangles()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/Cycle_Triangles.md)/`IN()`
  called bare `t(Adj)`, which only dispatches to `Matrix::t()` correctly
  if the `Matrix` package happens to be attached; fixed to call
  `Matrix::t()` explicitly.

### multideconv / OmnipathR / liana integration

- [`prepare_input_files()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/prepare_input_files.md)
  now calls
  [`multideconv::compute.deconvolution.analysis()`](https://rdrr.io/pkg/multideconv/man/compute.deconvolution.analysis.html)
  to identify and collapse correlated cell-type subgroups (the actual
  “subgroups” feature of multideconv), rather than using raw,
  unclustered per-method deconvolution columns directly.
- Replaced the dependency on
  [`multideconv::estimate_expression_profiles()`](https://rdrr.io/pkg/multideconv/man/estimate_expression_profiles.html)
  (removed from multideconv’s current source, and mathematically
  unsound) with a simple, correct per-gene non-negative least squares
  estimator.
- Cell-type names are now cleaned of method/signature prefixes (e.g.
  `"Quantiseq_TIL10_B.cells"` -\> `"B.cells"`) after subgroup analysis.
- Switched the ligand-receptor prior network from a generic OmnipathR
  confidence filter to
  [`liana::get_curated_omni()`](https://saezlab.github.io/liana/reference/get_curated_omni.html),
  matching the original R notebook’s curated resource.

### New functions

- [`correlate_features_with_score()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/correlate_features_with_score.md)
  and
  [`correlation_plot()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/correlation_plot.md):
  Spearman correlation between network features and a continuous
  per-patient score (e.g. an immune response score), with pooled and
  per-group results and a rainfall-style visualization.

### Other changes

- [`wilcox_group_test()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/wilcox_group_test.md)
  now reports fold change and `Log2FC` per feature and runs every
  pairwise comparison automatically for more than two groups;
  `volcano_plot()` plots `log2(fold change)` instead of the raw Wilcoxon
  statistic.
- Rewrote
  [`Wedges()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/Wedges.md),
  [`Trust_Triangles()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/Trust_Triangles.md),
  and
  [`Cycle_Triangles()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/Cycle_Triangles.md)
  to cache neighbor lookups and vectorize row construction instead of
  repeated sparse-matrix row slicing and O(n^2) list growth;
  cross-validated as exactly equivalent to the previous implementation
  on 900 random-graph trials, roughly 100-1000x faster on realistic
  graph sizes.
- Added a `testthat` regression suite covering the fixes above.
- Added `stats`, `utils`, and `nnls` to package dependencies; removed
  stray committed run artifacts (`Results/*.csv`, `Rplots.pdf`).

## RaCInG 0.1.0

- Initial R port of the original Python RaCInG implementation.
