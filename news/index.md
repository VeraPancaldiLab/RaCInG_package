# Changelog

## RaCInG 0.2.0

### Critical fixes

- **[`model1()`](https://mhurtado13.github.io/racing/reference/model1.md)
  /
  [`genRandomEdgeList()`](https://mhurtado13.github.io/racing/reference/genRandomEdgeList.md)
  ignored the requested average degree.**
  [`genRandomEdgeList()`](https://mhurtado13.github.io/racing/reference/genRandomEdgeList.md)
  never received the caller’s edge count and instead always sampled
  `nrow(Dligrec) * ncol(Dligrec)` edges, regardless of `N` or
  `av`/`avdeg`. Every Monte Carlo simulation run against a real
  ligand-receptor network generated graphs with the wrong density (often
  drastically so) in either direction, unrelated to the requested
  average degree. Simulation results computed with earlier versions
  should be treated as invalid and re-run.
- **[`TarjanIterative()`](https://mhurtado13.github.io/racing/reference/TarjanIterative.md)
  ([`GSCC()`](https://mhurtado13.github.io/racing/reference/GSCC.md)/[`IN()`](https://mhurtado13.github.io/racing/reference/IN.md)/[`OUT()`](https://mhurtado13.github.io/racing/reference/OUT.md)/[`countGSCC()`](https://mhurtado13.github.io/racing/reference/countGSCC.md))
  was broken.** A scoping bug (`stack <- ...` instead of `<<-` inside a
  closure) desynced the SCC stack from tracked state, and
  `seq(i + 1, length(neigh))` produced a descending sequence for any
  zero-out-degree vertex - both caused crashes or silently wrong
  strongly-connected-component results on ordinary graphs. Rewritten as
  a standard iterative Tarjan and cross-validated against a brute-force
  reachability oracle on 200 random graphs.
- **Kernel-based Wedge/Triangle/GSCC features double-counted shared
  nodes.**
  [`compute_kernel()`](https://mhurtado13.github.io/racing/reference/compute_kernel.md)
  baked cell-type abundance into both kernel indices, which is correct
  for Direct communication but double-weights any node shared by more
  than one edge in a motif (the middle node of a wedge, all three nodes
  of a triangle, and - since the underlying branching-process equations
  are nonlinear - every GSCC computation regardless of normalization).
  [`compute_kernel()`](https://mhurtado13.github.io/racing/reference/compute_kernel.md)
  now returns the unweighted structural kernel, and each feature
  function multiplies in abundance explicitly, once per node, matching
  the original implementation.
- **Triangle features were mislabeled.** The unbundled “trust triangle”
  (TT) computation actually computed the cycle-triangle (CT) formula,
  and the true trust-triangle metric was missing entirely; the bundled
  8-term formula dropped 6 of its 8 terms.
  [`computeTriangles()`](https://mhurtado13.github.io/racing/reference/computeTriangles.md)
  now returns both `TT_*` and `CT_*` columns when unbundled, and the
  correct 8-term sum when bundled.
- Self-loops were not excluded before Monte Carlo motif counting
  ([`EdgetoAdj()`](https://mhurtado13.github.io/racing/reference/EdgetoAdj.md)
  was used where
  [`EdgetoAdj_No_loop()`](https://mhurtado13.github.io/racing/reference/EdgetoAdj_No_loop.md)
  was required), inflating direct/wedge/triangle counts.
- [`Find_Number_Trust_Triangles_Unique()`](https://mhurtado13.github.io/racing/reference/Find_Number_Trust_Triangles_Unique.md)
  compared against the raw wedge count instead of its sign,
  undercounting whenever two directly-connected nodes had two or more
  common neighbors.
- Cell-type identifiers (patient row names) were silently dropped from
  `Cmatrix` in
  [`createCellTypeDistr()`](https://mhurtado13.github.io/racing/reference/createCellTypeDistr.md),
  breaking `patient_names` in downstream outputs.
- [`Read_Sim_Output()`](https://mhurtado13.github.io/racing/reference/Read_Sim_Output.md)
  inferred the number of cell types via
  [`ncol()`](https://rdrr.io/r/base/nrow.html) on a ragged, NA-padded
  data frame, which broke for datasets with fewer than 4 cell types.
- Empty motif lists
  ([`Wedges()`](https://mhurtado13.github.io/racing/reference/Wedges.md)/[`Trust_Triangles()`](https://mhurtado13.github.io/racing/reference/Trust_Triangles.md)/[`Cycle_Triangles()`](https://mhurtado13.github.io/racing/reference/Cycle_Triangles.md)
  returning zero rows) crashed
  [`Count_Types()`](https://mhurtado13.github.io/racing/reference/Count_Types.md).
- [`Cycle_Triangles()`](https://mhurtado13.github.io/racing/reference/Cycle_Triangles.md)/[`IN()`](https://mhurtado13.github.io/racing/reference/IN.md)
  called bare `t(Adj)`, which only dispatches to `Matrix::t()` correctly
  if the `Matrix` package happens to be attached; fixed to call
  `Matrix::t()` explicitly.

### multideconv / OmnipathR / liana integration

- [`prepare_input_files()`](https://mhurtado13.github.io/racing/reference/prepare_input_files.md)
  now calls
  [`multideconv::compute.deconvolution.analysis()`](https://rdrr.io/pkg/multideconv/man/compute.deconvolution.analysis.html)
  to identify and collapse correlated cell-type subgroups (the actual
  “subgroups” feature of multideconv), rather than using raw,
  unclustered per-method deconvolution columns directly.
- Replaced the dependency on
  `multideconv::estimate_expression_profiles()` (removed from
  multideconv’s current source, and mathematically unsound) with a
  simple, correct per-gene non-negative least squares estimator.
- Cell-type names are now cleaned of method/signature prefixes (e.g.
  `"Quantiseq_TIL10_B.cells"` -\> `"B.cells"`) after subgroup analysis.
- Switched the ligand-receptor prior network from a generic OmnipathR
  confidence filter to
  [`liana::get_curated_omni()`](https://saezlab.github.io/liana/reference/get_curated_omni.html),
  matching the original R notebook’s curated resource.

### New functions

- [`correlate_features_with_score()`](https://mhurtado13.github.io/racing/reference/correlate_features_with_score.md)
  and
  [`correlation_plot()`](https://mhurtado13.github.io/racing/reference/correlation_plot.md):
  Spearman correlation between network features and a continuous
  per-patient score (e.g. an immune response score), with pooled and
  per-group results and a rainfall-style visualization.

### Other changes

- [`wilcox_group_test()`](https://mhurtado13.github.io/racing/reference/wilcox_group_test.md)
  now reports fold change and `Log2FC` per feature and runs every
  pairwise comparison automatically for more than two groups;
  [`volcano_plot()`](https://mhurtado13.github.io/racing/reference/volcano_plot.md)
  plots `log2(fold change)` instead of the raw Wilcoxon statistic.
- Rewrote
  [`Wedges()`](https://mhurtado13.github.io/racing/reference/Wedges.md),
  [`Trust_Triangles()`](https://mhurtado13.github.io/racing/reference/Trust_Triangles.md),
  and
  [`Cycle_Triangles()`](https://mhurtado13.github.io/racing/reference/Cycle_Triangles.md)
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
