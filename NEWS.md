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
