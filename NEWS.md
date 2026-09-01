# RaCInG 0.3.0

## Breaking change: `norm` now defaults to `FALSE`

* `compute_racing_kernel()` and `compute_racing_montecarlo()` previously
  defaulted to `norm = TRUE` (a normalized enrichment ratio against a
  uniformized baseline), inconsistent with the lower-level functions they
  wrap (`compute_kernel()`'s `normalize`, `runSim()`'s `norm`), which already
  defaulted to `FALSE`. Since the normalized value is a ratio of two Monte
  Carlo estimates, it is inherently noisier than the raw value at a given
  `Ngraphs` -- in practice this meant the previous default silently required
  a much larger `Ngraphs` than most callers were using to detect real
  differences (confirmed on real data: a comparison with clear signal in the
  raw features showed *zero* significant features after normalization at
  `Ngraphs = 10`). Both functions now default to `norm = FALSE` (the raw,
  abundance-weighted communication magnitude), matching the rest of the
  package; pass `norm = TRUE` explicitly when you specifically want the
  enrichment-ratio behavior, and budget a larger `Ngraphs` for the Monte
  Carlo workflow if you do.

## Fixed a real bug behind blank/incomplete `top_features_plot()`/`correlation_plot()` output

* `wilcox.test()` returns an `NA` p-value for a feature that is (near-)
  constant within one or both groups -- common after normalization, where
  many features collapse to the same `0/0` fallback value. `top_features_plot()`'s
  filter (`wilcox_results[is.finite(Log2FC) & Adjusted_P_value < p_threshold, ]`)
  did not account for this: in R, `TRUE & NA` is `NA`, and indexing a data
  frame with a logical vector containing `NA` **inserts an all-`NA` row**
  at that position rather than dropping it (unlike `FALSE`, which correctly
  excludes the row) -- so those NA-p-value features silently became blank,
  labelless bars in the output instead of being excluded. Fixed by coercing
  `NA` to `FALSE` in the filter before subsetting. `correlation_plot()` had
  the same latent risk via `cor.test()`'s NA `Rho` for constant features
  (surfaced via `order()`'s NA-last default rather than row-insertion, but
  fixed the same way: explicit `!is.na()` filtering before ranking).

## Fixed a severe O(n^2) inefficiency in `countDirect()`

* `countDirect()`'s cell-type-pair edge-counting step re-sliced the full
  `N x N` adjacency matrix from scratch for every one of `nCells^2` cell-type
  pairs per simulated graph -- redundant even within its own inner loop,
  since the row-slice for a given source type doesn't depend on the target
  type, yet was recomputed once per target type anyway. At 106 cell types
  this is 11,236 redundant matrix re-slices per graph. Replaced with a single
  sparse matrix product (`t(type_indicator) %*% Adj %*% type_indicator`,
  where `type_indicator` is a one-hot cell-type indicator matrix), which
  computes the full cell-type x cell-type edge-count matrix in one pass.
  Verified numerically identical to the previous implementation; **~182x**
  faster on real data (106 cell types, 500 cells, average degree 20: 1.22s ->
  0.007s per graph).

## Removed

* Removed 11 functions confirmed unused anywhere in the package, its tests,
  vignette, or README, and unreferenced by any downstream project: the
  count-only `Find_Number_Trust_Triangles_Unique()`, `Find_Number_Triangles()`,
  `Find_Number_Triangles_Unique()`, `Find_Number_2Loops()`,
  `Find_Number_2Loops_Unique()`, `Find_Number_Wedges()`,
  `Find_Number_Wedges_Unique()` (superseded by `Trust_Triangles()`/
  `Cycle_Triangles()`/`Wedges()`, which already return the count plus the
  full motif list); the bow-tie decomposition `IN()`/`OUT()` and their sole
  dependency `BFS()`; `generateUniformLRGraph()` (superseded by the `norm`/
  `normalize` logic built directly into `runSim()`/`compute_kernel()`); and
  `model1()`'s `genRandom = TRUE` demo path along with the 4 fabricated-data
  generators it called (`genRandomCellTypeDistr()`, `genRandomLigRecDistr()`,
  `genRandomCellLigands()`, `genRandomCellReceptors()`) -- `model1()` now
  requires real `cellLigList`/`cellRecList`/`Dcelltype`/`Dligrec` arguments.
  Several of these were direct ports of code the original Python
  implementation's own authors had already set aside as test-only scaffolding
  (`distribution_generation.py`, explicitly documented there as "not used in
  paper").
* Removed `volcano_plot()` (labels were dropped for busy plots and it
  couldn't scale to features with long names); replaced by
  `top_features_plot()`, a ranked horizontal bar chart of the top up/down
  features by `log2(fold change)`.

## Multiple feature types without recomputing

* `compute_racing_kernel()`/`compute_kernel_features()` now accept
  `communication_type` as a character vector (e.g.
  `c("D", "W", "TT", "GSCC")`). The kernel is still computed exactly once per
  call regardless of how many types are requested; `features` becomes a
  named list (one data frame per type) when more than one is given.
* `compute_racing_montecarlo()`/`runSim()` now accept `communication_type` as
  a character vector too (any combination of `"D"`, `"W"`, `"TT"`, `"CT"`,
  `"GSCC"`), via a new `countAllTypes()` that extracts every requested type
  from the *same* simulated graphs instead of resimulating a fresh set of
  graphs per type -- graph generation, not feature extraction, is the
  expensive part of a Monte Carlo run. `output` becomes a named list (one
  entry per type, one `.out`/CSV file pair per type on disk) when more than
  one is given. Single-type usage is unchanged (verified bit-identical
  against the previous single-type implementation with a fixed seed).

## Fixed a catastrophic O(n^3) crash/hang in kernel Wedge/Triangle features

* `calculateWedges()` and `computeTriangles()` built one named-list entry per
  cell-type combination inside a triple-nested loop, then converted the
  whole list to a data frame -- for realistic cell-type counts (e.g. 100+)
  this reaches 10^5-10^6+ individual list insertions and made `"W"`/`"TT"`
  effectively unusable via the kernel method, unlike `"D"`. Rewritten as
  vectorized array operations (index grid + flattened-matrix lookups,
  mirroring the equivalent Monte Carlo fix below); verified numerically
  identical to the original loop-based output across bundled/unbundled and
  raw/normalized cases. At 106 real cell types, `"TT"` (previously
  unusable) now completes in a few seconds.
* `calculateDirect()`, `calculateWedges()`, and `computeTriangles()` now drop
  cell-type combinations with no possible ligand-receptor pathway in any
  patient (structurally zero for every patient) instead of returning them as
  dead all-zero columns -- for `calculateWedges()`/`computeTriangles()` this
  filtering happens *before* the abundance-weighting/ratio computation, so it
  also reduces the actual memory/compute cost, not just the output size.

## Monte Carlo memory and parallelization fixes

* `runSim()`'s `ncores > 1` path (`parallel::makeCluster()` +
  `doParallel::registerDoParallel()` + `foreach`, replacing the previous
  fork-based `parallel::mclapply()` for Windows portability) was shipping a
  full copy of the entire patient-by-ligand-by-receptor `LRmatrix` tensor
  to *every* worker process, regardless of patient count -- for a 106
  cell-type, 116-patient cohort this meant ~460MB duplicated per worker even
  for a trivial single-patient run. Fixed by pre-slicing `LRmatrix`/`Cmatrix`
  per patient and making those small slices the `foreach` iteration variable
  (dispatched per task) instead of a name `foreach` broadcasts whole to every
  worker; confirmed via monitored testing this dropped per-worker overhead
  from ~460MB to a few MB.
* `compute_results_processing()`'s `remove_direction` column-merging grouped
  and relabeled columns by splitting the pasted column-name string on `"_"`,
  which silently scrambles cell-type identity whenever a cell-type name
  itself contains `"_"` (e.g. multideconv's method-tagged names like
  `"CBSX_CBSX.HNSCC.scRNAseq_B.cells"`) -- there is no way to tell, from the
  string alone, which underscores separate the two cell types in a pair from
  underscores inside one cell type's own name. Fixed by grouping on the
  underlying integer index grid instead of re-parsing labels; added a
  regression test.
* `calculateDirect()`/`calculateWedges()` silently dropped real patient
  names (always fell back to `"Patient_1"`, `"Patient_2"`, ...) even when
  `compute_kernel_features()`/`compute_racing_kernel()` had real ones
  available -- `"TT"`/`"GSCC"` already used them correctly. Fixed, and kernel
  and Monte Carlo feature matrices now use identical column-naming
  conventions (`"Dir_A_B"`, `"W_A_B_C"`, etc.) too, so their outputs can be
  fed into the same downstream statistical/plotting functions
  interchangeably.

## Other changes

* Added `ncores` to `runSim()`/`compute_racing_montecarlo()` for patient-level
  parallelization.
* Added `already_normalized` and `expr_threshold` parameters to
  `prepare_input_files()`.
* Vectorized several remaining O(n^2)-O(n^3) hotspots: `prepare_input_files()`'s
  `ccc_table` construction, `Read_Sim_Output()`'s `.out`-file parsing, and
  `runSim()`'s composition-block writing.
* `createCellTypeDistr()` now drops `Cmatrix` columns not present in
  `Lmatrix`/`Rmatrix` and renormalizes, instead of leaving mismatched columns
  in place.
* `calculateDirect()`/`calculateWedges()` fill structurally-impossible (0/0)
  normalized ratios with `0`.
* Expanded the `testthat` regression suite to cover all of the above.

# RaCInG 0.2.0

## Critical fixes

* **`model1()` / `genRandomEdgeList()` ignored the requested average degree.**
  `genRandomEdgeList()` never received the caller's edge count and instead
  always sampled `nrow(Dligrec) * ncol(Dligrec)` edges, regardless of `N` or
  `av`/`avdeg`. Every Monte Carlo simulation run against a real ligand-receptor
  network generated graphs with the wrong density (often drastically so) in
  either direction, unrelated to the requested average degree. Simulation
  results computed with earlier versions should be treated as invalid and
  re-run.
* **`TarjanIterative()` (`GSCC()`/`IN()`/`OUT()`/`countGSCC()`) was broken.**
  A scoping bug (`stack <- ...` instead of `<<-` inside a closure) desynced
  the SCC stack from tracked state, and `seq(i + 1, length(neigh))` produced
  a descending sequence for any zero-out-degree vertex - both caused crashes
  or silently wrong strongly-connected-component results on ordinary graphs.
  Rewritten as a standard iterative Tarjan and cross-validated against a
  brute-force reachability oracle on 200 random graphs.
* **Kernel-based Wedge/Triangle/GSCC features double-counted shared nodes.**
  `compute_kernel()` baked cell-type abundance into both kernel indices,
  which is correct for Direct communication but double-weights any node
  shared by more than one edge in a motif (the middle node of a wedge, all
  three nodes of a triangle, and - since the underlying branching-process
  equations are nonlinear - every GSCC computation regardless of
  normalization). `compute_kernel()` now returns the unweighted structural
  kernel, and each feature function multiplies in abundance explicitly, once
  per node, matching the original implementation.
* **Triangle features were mislabeled.** The unbundled "trust triangle" (TT)
  computation actually computed the cycle-triangle (CT) formula, and the true
  trust-triangle metric was missing entirely; the bundled 8-term formula
  dropped 6 of its 8 terms. `computeTriangles()` now returns both `TT_*` and
  `CT_*` columns when unbundled, and the correct 8-term sum when bundled.
* Self-loops were not excluded before Monte Carlo motif counting
  (`EdgetoAdj()` was used where `EdgetoAdj_No_loop()` was required),
  inflating direct/wedge/triangle counts.
* `Find_Number_Trust_Triangles_Unique()` compared against the raw wedge count
  instead of its sign, undercounting whenever two directly-connected nodes
  had two or more common neighbors.
* Cell-type identifiers (patient row names) were silently dropped from
  `Cmatrix` in `createCellTypeDistr()`, breaking `patient_names` in
  downstream outputs.
* `Read_Sim_Output()` inferred the number of cell types via `ncol()` on a
  ragged, NA-padded data frame, which broke for datasets with fewer than 4
  cell types.
* Empty motif lists (`Wedges()`/`Trust_Triangles()`/`Cycle_Triangles()`
  returning zero rows) crashed `Count_Types()`.
* `Cycle_Triangles()`/`IN()` called bare `t(Adj)`, which only dispatches to
  `Matrix::t()` correctly if the `Matrix` package happens to be attached;
  fixed to call `Matrix::t()` explicitly.

## multideconv / OmnipathR / liana integration

* `prepare_input_files()` now calls `multideconv::compute.deconvolution.analysis()`
  to identify and collapse correlated cell-type subgroups (the actual
  "subgroups" feature of multideconv), rather than using raw, unclustered
  per-method deconvolution columns directly.
* Replaced the dependency on `multideconv::estimate_expression_profiles()`
  (removed from multideconv's current source, and mathematically unsound) with
  a simple, correct per-gene non-negative least squares estimator.
* Cell-type names are now cleaned of method/signature prefixes (e.g.
  `"Quantiseq_TIL10_B.cells"` -> `"B.cells"`) after subgroup analysis.
* Switched the ligand-receptor prior network from a generic OmnipathR
  confidence filter to `liana::get_curated_omni()`, matching the original
  R notebook's curated resource.

## New functions

* `correlate_features_with_score()` and `correlation_plot()`: Spearman
  correlation between network features and a continuous per-patient score
  (e.g. an immune response score), with pooled and per-group results and a
  rainfall-style visualization.

## Other changes

* `wilcox_group_test()` now reports fold change and `Log2FC` per feature and
  runs every pairwise comparison automatically for more than two groups;
  `volcano_plot()` plots `log2(fold change)` instead of the raw Wilcoxon
  statistic.
* Rewrote `Wedges()`, `Trust_Triangles()`, and `Cycle_Triangles()` to cache
  neighbor lookups and vectorize row construction instead of repeated
  sparse-matrix row slicing and O(n^2) list growth; cross-validated as
  exactly equivalent to the previous implementation on 900 random-graph
  trials, roughly 100-1000x faster on realistic graph sizes.
* Added a `testthat` regression suite covering the fixes above.
* Added `stats`, `utils`, and `nnls` to package dependencies; removed stray
  committed run artifacts (`Results/*.csv`, `Rplots.pdf`).

# RaCInG 0.1.0

* Initial R port of the original Python RaCInG implementation.
