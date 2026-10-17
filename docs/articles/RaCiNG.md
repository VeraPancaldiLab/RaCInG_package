# Getting started with RaCInG

## Overview

`RaCInG` reconstructs patient-specific cell-cell communication networks
from bulk RNA-seq data. The package supports two complementary
workflows:

- a **kernel-based** approach for fast, deterministic
  direct-communication features (`"D"`), and
- a **Monte Carlo** approach that simulates networks and reports direct
  communication (`"D"`) and/or the giant strongly connected component
  (`"GSCC"`).

This vignette shows the recommended workflow and the most important
entry points for new users, using both the bundled, ready-made example
data and a real bulk RNA-seq counts matrix run through the full
preprocessing pipeline.

## Installation

### Install from GitHub

``` r
# install.packages("remotes")
remotes::install_github("VeraPancaldiLab/RaCInG_package")
library(RaCInG)
```

### Install from a local checkout

``` r
# install.packages("devtools")
devtools::install(".")
library(RaCInG)
```

If you want to build the RaCInG input matrices directly from raw counts,
install the optional preprocessing dependencies used by
[`prepare_input_files()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/prepare_input_files.md):

``` r
install.packages(c("ggplot2", "nnls"))
# ADImpute and OmnipathR are available from Bioconductor:
BiocManager::install(c("ADImpute", "OmnipathR"))
# liana and multideconv are GitHub-only:
remotes::install_github(c("saezlab/liana", "VeraPancaldiLab/multideconv"))
```

## Workflow at a glance

| Goal                                                   | Function                                                                                                                                                                                                                                  | Output                                                  |
|--------------------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------|
| Build input matrices from raw counts                   | [`prepare_input_files()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/prepare_input_files.md)                                                                                                                              | Named list with `L`, `R`, `C`, `LR` matrices and labels |
| Compute direct communication (deterministic)           | [`compute_racing_kernel()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/compute_racing_kernel.md)                                                                                                                          | Kernel arrays + feature matrix                          |
| Compute direct communication and/or GSCC by simulation | [`compute_racing_montecarlo()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/compute_racing_montecarlo.md)                                                                                                                  | Processed Monte Carlo results                           |
| Decompose a direct score into ligand-receptor pairs    | [`computeLRContributions()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/computeLRContributions.md) + [`plotLRContributions()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/plotLRContributions.md)         | Ranked contributions and a LIANA-style dot plot         |
| Compare clinical groups                                | [`wilcox_group_test()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/wilcox_group_test.md) + [`top_features_plot()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/top_features_plot.md)                       | Statistics table and ranked bar plot                    |
| Relate features to a continuous score                  | [`correlate_features_with_score()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/correlate_features_with_score.md) + [`correlation_plot()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/correlation_plot.md) | Correlation table and rainfall plot                     |

## How `prepare_input_files()` builds its inputs

Starting from a raw gene-by-sample counts matrix,
[`prepare_input_files()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/prepare_input_files.md)
chains together several packages so you never have to run them by hand:

1.  **Normalization**: counts are converted to TPM with
    [`ADImpute::NormalizeTPM()`](https://rdrr.io/pkg/ADImpute/man/NormalizeTPM.html).
2.  **Deconvolution**: cell-type fractions are estimated with
    [`multideconv::compute.deconvolution()`](https://verapancaldilab.github.io/multideconv/reference/compute.deconvolution.html)
    (defaults to the Quantiseq method, but any method(s) supported by
    multideconv can be passed via `deconv_method`), which also
    standardizes cell-type names to multideconv’s vocabulary
    ([`multideconv::get_cell_type_nomenclature()`](https://verapancaldilab.github.io/multideconv/reference/get_cell_type_nomenclature.html)).
3.  **Subgroup analysis**:
    [`multideconv::compute.deconvolution.analysis()`](https://verapancaldilab.github.io/multideconv/reference/compute.deconvolution.analysis.html)
    (multideconv \>= 0.2.0) clusters highly correlated deconvolution
    features of the same cell type (e.g. several methods agreeing on it)
    into subgroups such as `"B.cells_Subgroup.1"`, discarding
    mostly-zero and low-variance features along the way. Features that
    join no subgroup keep their full name
    (e.g. `"Quantiseq_TIL10_B.cells"`); these names are used as the
    RaCInG cell types.
4.  **Cell-type expression profiles**: for every gene, a non-negative
    least squares fit of bulk expression against the cell-type fractions
    of each method-signature combination gives one expression value per
    gene per cell type; each subgroup’s profile is the median of its
    members’ profiles.
5.  **Ligand-receptor prior knowledge**: a curated, OmniPath-backed
    consensus network is retrieved with
    [`liana::get_curated_omni()`](https://saezlab.github.io/liana/reference/get_curated_omni.html)
    and restricted to ligand/receptor roles, then decomplexified into
    single-gene subunits.
6.  **Filtering and assembly**: ligand-receptor pairs are kept only
    where both genes reach `expr_threshold` (default 10 TPM) in the
    sending and receiving cell types, and the `Lmatrix`, `Rmatrix`,
    `Cmatrix`, and `LRmatrix` inputs (see below) are assembled and
    written to `output_folder`. Each ligand-receptor pair’s per-patient
    strength combines the two genes’ bulk `log2(TPM + 1)` expression
    with `fun_LR` (default `min`, i.e. the weaker partner limits the
    interaction).

You can skip any of steps 2-5 by supplying your own `deconv`,
`cell_expr_profile`, or `cc_network` directly. A user-supplied `deconv`
must have exactly `colnames(counts)` as row names, in the same order,
and should use multideconv’s cell-type nomenclature (a warning lists any
names outside it). If your expression data is already linear-scale TPM
rather than raw counts, pass `already_normalized = TRUE`; log-scale data
(e.g. `log2(TPM + 1)`) must be converted back first (`2^x - 1`), since
`expr_threshold` is compared on the linear scale.

## Recommended workflow with a real bulk RNA-seq dataset

The examples below use the raw bulk RNA-seq counts matrix bundled with
the `multideconv` package (`raw_counts`, from Mariathasan et al. 2018, a
urothelial/bladder cancer cohort profiled for anti-PD-L1 response) to
exercise the full pipeline end to end, from raw counts to features to
statistics. These chunks call out to the network (OmniPath/LIANA) and
run real deconvolution, so they are shown but not evaluated when the
vignette is built; run them interactively to try the full pipeline
yourself.

### 1. Build input matrices from raw counts

``` r
data(raw_counts, package = "multideconv")
counts_matrix <- as.matrix(raw_counts)

input <- prepare_input_files(
  counts = counts_matrix,
  output_folder = "Results/",
  file_name = "mariathasan",
  deconv_method = "Quantiseq"
)

str(input, max.level = 1)
input$celltypes
```

### 2. Run the kernel method

The kernel method is the fastest way to derive direct communication
features (`communication_type = "D"`, the only type it computes) across
patients. You can pass `counts` to let the function compute inputs
automatically, or supply previously computed matrices via `input_data`.

``` r
# Option A: from raw counts (runs prepare_input_files internally)
kernel_res <- compute_racing_kernel(
  counts = counts_matrix,
  output_folder = tempdir(),
  file_name = "mariathasan",
  communication_type = "D",
  norm = TRUE
)

# Option B: from pre-computed input matrices (skips preprocessing)
kernel_res <- compute_racing_kernel(
  input_data = input,
  communication_type = "D",
  norm = TRUE
)

head(kernel_res$features[, 1:5])
```

Features are named `"Dir_A_B"`, with the lower-index (alphabetically
first) cell type first; with `remove_direction = TRUE` (the default) the
two directions of each pair are summed into one feature. Feature columns
with no possible ligand-receptor pathway in any patient (structurally
zero for every patient) are dropped automatically, rather than returned
as dead all-zero columns – so `ncol(kernel_res$features)` is usually
well below the full combinatorial count of cell-type pairs. Lower
`zero_threshold` (e.g. `0.9`) to also drop features that are zero in
most patients.

#### Raw or normalized features?

With `norm = FALSE` (the default) features are the raw communication
strength, which grows with both the abundance of the two cell types and
the expression of their ligand-receptor pairs. With `norm = TRUE`,
RaCInG also builds a baseline in which every possible ligand-receptor
pair is equally strong and reports the ratio of the two: an enrichment
score that does not depend on cell-type abundance and isolates which
pairs communicate more (or less) than the network structure alone
predicts. In the Monte Carlo method the ratio of two simulated estimates
is noisier than the raw value, so budget a larger `Ngraphs` when you use
`norm = TRUE` there.

### 3. Run the Monte Carlo method

Use the Monte Carlo workflow when you want simulation-based summaries or
uncertainty estimates from repeated graph realizations. The same
`input_data` shortcut is available here. Simulation cost scales with
`Ncells` and `Ngraphs`; the values below are kept small for a quick
illustration.

``` r
set.seed(1)
mc_res <- compute_racing_montecarlo(
  counts = counts_matrix,
  output_folder = tempdir(),
  deconv_method = "Quantiseq",
  file_name = "mariathasan",
  nPatients = 3,
  communication_type = "D",
  Ncells = 100,
  Ngraphs = 10,
  Ndegree = 3,
  norm = TRUE
)

# Or from pre-computed inputs:
mc_res <- compute_racing_montecarlo(
  input_data = input,
  output_folder = tempdir(),
  file_name = "mariathasan",
  communication_type = "D",
  Ncells = 100,
  Ngraphs = 10,
  Ndegree = 3,
  norm = TRUE
)

head(mc_res$output$mean[, 1:5])
```

`mc_res$output` holds the per-patient mean (`mean`) and standard
deviation (`sd`) across the simulated graphs. The raw simulations
(`<file_name>.out`, and `<file_name>_norm.out` with `norm = TRUE`) and
the processed features (`<file_name>.csv`, `<file_name>_sd.csv`) are
also written to `output_folder`; they can be re-read later with
[`Read_Sim_Output()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/Read_Sim_Output.md)
and
[`compute_results_processing()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/compute_results_processing.md).

`communication_type` can also be `"GSCC"` (the share of cells, per cell
type, in the giant strongly connected component) or both,
`c("D", "GSCC")`. Both are then extracted from the *same* simulated
graphs, instead of re-simulating a fresh set of graphs per type – graph
generation, not feature extraction, is the expensive part of a Monte
Carlo run:

``` r
mc_res <- compute_racing_montecarlo(
  input_data = input,
  communication_type = c("D", "GSCC"),
  Ncells = 100,
  Ngraphs = 10,
  Ndegree = 3,
  norm = TRUE,
  output_folder = tempdir(),
  file_name = "mariathasan"
)

# output is now a named list, one entry per requested type
names(mc_res$output)
head(mc_res$output$D$mean[, 1:5])
```

Patients are independent of each other, so larger runs can be
parallelized across cores with `ncores` (via
`parallel`/`doParallel`/`foreach`, so it works the same way on Windows,
macOS, and Linux):

``` r
mc_res <- compute_racing_montecarlo(
  input_data = input,
  communication_type = "D",
  Ncells = 500,
  Ngraphs = 100,
  Ndegree = 20,
  ncores = 4,
  output_folder = tempdir(),
  file_name = "mariathasan"
)
```

### 4. Compare clinical groups

Once features are available in a patient-by-feature matrix, use the
built-in Wilcoxon workflow to compare groups.
[`wilcox_group_test()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/wilcox_group_test.md)
reports, per feature and per pair of groups, the Wilcoxon statistic,
fold change, and (jointly) FDR-adjusted p-value;
[`top_features_plot()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/top_features_plot.md)
shows the top up/down features by `log2(fold change)` as a ranked bar
chart.

``` r
# One label per patient (row of the feature matrix); replace this placeholder
# with your own clinical annotation, e.g. response to therapy.
grouping <- rep(c("Responder", "Non-responder"), length.out = nrow(kernel_res$features))
wilcox_results <- wilcox_group_test(kernel_res$features, grouping)
head(wilcox_results)
top_features_plot(wilcox_results, top_n = 15)
```

With more than two groups, every pairwise comparison is run
automatically; subset `wilcox_results` by `Comparison` (or pass a
two-group subset to
[`top_features_plot()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/top_features_plot.md))
to inspect one pair at a time.

### 5. Correlate features with a continuous score

If you have a continuous per-patient score (e.g. an immune response or
survival score),
[`correlate_features_with_score()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/correlate_features_with_score.md)
reports the Spearman correlation between every feature and the score,
pooled across all patients and, optionally, within groups.
[`correlation_plot()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/correlation_plot.md)
shows the top positive and negative correlations.

``` r
response_score <- setNames(rnorm(nrow(kernel_res$features)), rownames(kernel_res$features))
corr_results <- correlate_features_with_score(kernel_res$features, response_score)
head(corr_results)
correlation_plot(corr_results, top_n = 10)
```

## Understanding the input files

RaCInG requires four matrices and associated label vectors that describe
the cell-cell communication landscape for a cohort of patients. The
table below summarises each component:

| Component       | Dimensions                     | Description                                                                                                                                                     |
|-----------------|--------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Lmatrix**     | cell types × ligands           | Binary compatibility matrix: `1` if that ligand is expressed above `expr_threshold` in that cell type, `0` otherwise. Rows are cell types; columns are ligands. |
| **Rmatrix**     | cell types × receptors         | Binary compatibility matrix: `1` if that receptor is expressed above `expr_threshold` in that cell type, `0` otherwise. Same row order as `Lmatrix`.            |
| **Cmatrix**     | patients × cell types          | Cell-type fraction for each patient. Each row sums to 1.                                                                                                        |
| **LRmatrix**    | ligands × receptors × patients | 3-D tensor of ligand–receptor interaction weights. Each patient slice is normalised to sum to 1.                                                                |
| **celltypes**   | character vector               | Alphabetically sorted cell-type names (shared across L, R, and C).                                                                                              |
| **ligands**     | character vector               | Ligand names matching the columns of `Lmatrix` and the first axis of `LRmatrix`.                                                                                |
| **receptors**   | character vector               | Receptor names matching the columns of `Rmatrix` and the second axis of `LRmatrix`.                                                                             |
| **Sign_matrix** | ligands × receptors            | Placeholder for stimulatory (+1) / inhibitory (−1) interaction signs; all zeros (unknown) by default. Not used by the current kernel or Monte Carlo features.   |

These inputs are typically generated by
[`prepare_input_files()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/prepare_input_files.md)
from a raw counts matrix, or they can be assembled manually from
pre-existing deconvolution and prior-network data. The bundled
`skcm_example` dataset provides a ready-made example of this structure.

### Inspecting the example inputs

``` r
library(RaCInG)
data(skcm_example)

# Overall structure
str(skcm_example, max.level = 1)
#> List of 8
#>  $ Lmatrix    : num [1:9, 1:276] 1 0 1 1 1 1 1 1 0 1 ...
#>   ..- attr(*, "dimnames")=List of 2
#>  $ Rmatrix    : num [1:9, 1:298] 1 0 1 1 0 1 1 1 0 0 ...
#>   ..- attr(*, "dimnames")=List of 2
#>  $ Cmatrix    : num [1:10, 1:9] 0.00374 0.01407 0.02119 0 0.00313 ...
#>   ..- attr(*, "dimnames")=List of 2
#>  $ LRmatrix   : num [1:276, 1:298, 1:10] 0.000228 0 0 0 0 ...
#>  $ celltypes  : chr [1:9] "B" "CAF" "CD8+ T" "DC" ...
#>  $ ligands    : chr [1:276] "LGALS9" "ADAM10" "TNFSF12" "ICOSLG" ...
#>  $ receptors  : chr [1:298] "PTPRC" "MET" "CD44" "LRP1" ...
#>  $ Sign_matrix: num [1:276, 1:298] 0 0 0 0 0 0 0 0 0 0 ...
```

``` r
# Lmatrix: 9 cell types × 276 ligands
dim(skcm_example$Lmatrix)
#> [1]   9 276
skcm_example$Lmatrix[1:4, 1:5]
#>      LGALS9 ADAM10 TNFSF12 ICOSLG TNF
#> [1,]      1      1       1      1   1
#> [2,]      0      1       1      0   0
#> [3,]      1      1       1      0   1
#> [4,]      1      1       1      1   1
```

``` r
# Rmatrix: 9 cell types × 298 receptors
dim(skcm_example$Rmatrix)
#> [1]   9 298
skcm_example$Rmatrix[1:4, 1:5]
#>      PTPRC MET CD44 LRP1 CD47
#> [1,]     1   0    1    0    1
#> [2,]     0   1    1    1    1
#> [3,]     1   0    1    0    1
#> [4,]     1   0    1    1    1
```

``` r
# Cmatrix: 10 patients × 9 cell types (rows sum to 1)
dim(skcm_example$Cmatrix)
#> [1] 10  9
skcm_example$Cmatrix[1:4, ]
#>                B         CAF       CD8+ T          DC        Endo          M
#> [1,] 0.003742534 0.017784413 0.0004763235 0.002384983 0.007165426 0.01515965
#> [2,] 0.014074525 0.019407444 0.0807858943 0.062715975 0.004718812 0.07897314
#> [3,] 0.021190628 0.007153891 0.0166497814 0.012495442 0.012443088 0.03645401
#> [4,] 0.000000000 0.038687885 0.0000000000 0.000942249 0.025184327 0.02306847
#>                NK       Treg     Tumor
#> [1,] 6.110609e-10 0.01089191 0.9423948
#> [2,] 6.783486e-04 0.02023155 0.7184143
#> [3,] 6.351805e-04 0.00000000 0.8929780
#> [4,] 4.680075e-09 0.00000000 0.9121171
```

``` r
# LRmatrix: 276 ligands × 298 receptors × 10 patients (3-D tensor)
dim(skcm_example$LRmatrix)
#> [1] 276 298  10
# First patient slice, top-left corner:
skcm_example$LRmatrix[1:4, 1:4, 1]
#>              [,1]        [,2]        [,3]        [,4]
#> [1,] 0.0002281731 0.001269455 0.001846635 0.001846635
#> [2,] 0.0000000000 0.001269455 0.003174809 0.000000000
#> [3,] 0.0000000000 0.000000000 0.000000000 0.000000000
#> [4,] 0.0000000000 0.000000000 0.000000000 0.000000000
```

``` r
# Label vectors
skcm_example$celltypes
#> [1] "B"      "CAF"    "CD8+ T" "DC"     "Endo"   "M"      "NK"     "Treg"  
#> [9] "Tumor"
head(skcm_example$ligands, 10)
#>  [1] "LGALS9"   "ADAM10"   "TNFSF12"  "ICOSLG"   "TNF"      "HLA.B"   
#>  [7] "HLA.DRA"  "HLA.DQA2" "HLA.DQA1" "HLA.DQB1"
head(skcm_example$receptors, 10)
#>  [1] "PTPRC"   "MET"     "CD44"    "LRP1"    "CD47"    "PTPRK"   "COLEC12"
#>  [8] "HAVCR2"  "MRC2"    "TSPAN15"
```

## Running with the bundled example data

The `skcm_example` list shown above can be passed directly to the kernel
or Monte Carlo workflows via the `input_data` parameter, and to the
statistical analysis functions once features have been computed.

### Kernel method on the example data

``` r
kernel_res <- compute_racing_kernel(
  input_data   = skcm_example,
  output_folder = tempdir(),
  communication_type = "D",
  norm = TRUE
)
#> Using pre-computed input matrices; skipping input generation.
#> Computing kernel for 10 patients
#> Calculating kernel...
#> Calculating features...

head(kernel_res$features[, 1:5])
#>             Dir_B_B Dir_B_CAF Dir_B_CD8+ T   Dir_B_DC Dir_B_Endo
#> Patient_1 0.3655633 0.6472004    0.2193690 0.33179690  0.6839505
#> Patient_2 0.7127917 0.5965352    1.0316271 0.84989477  0.6860757
#> Patient_3 1.1609442 0.7959182    0.8781016 0.92070117  0.8362820
#> Patient_4 0.4675364 0.4900159    0.1068774 0.09173062  0.3913845
#> Patient_5 0.2765913 0.5799559    0.1661719 0.37846202  0.4973818
#> Patient_6 0.4209898 0.5559673    0.7787024 0.51126999  0.6621874
```

### Decomposing a score into ligand-receptor pairs

[`compute_kernel()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/compute_kernel.md)
sums over every ligand-receptor pair to produce one score per cell-type
pair, so that score alone can’t tell you which ligand-receptor pairs
drove it.
[`computeLRContributions()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/computeLRContributions.md)
recovers that breakdown for one sender-receiver pair (exploiting that
the kernel formula is linear in `LRmatrix`, so each pair’s contribution
can be computed directly and sums back exactly to the original score),
and
[`plotLRContributions()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/plotLRContributions.md)
shows the top contributors as a LIANA-style dot plot (dot size = share
of that pair’s total score, color = raw contribution).

``` r
lr_contrib <- computeLRContributions(
  Lmatrix = skcm_example$Lmatrix, Rmatrix = skcm_example$Rmatrix,
  Cmatrix = skcm_example$Cmatrix, LRmatrix = skcm_example$LRmatrix,
  celltypes = skcm_example$celltypes, sender = "Tumor", receiver = "CD8+ T"
)
head(lr_contrib)
#>     Patient             Edge Sender Receiver Ligand Receptor Contribution
#> 1 Patient_1 Tumor <-> CD8+ T  Tumor   CD8+ T  TGFB1   TGFBR3 0.0010444341
#> 2 Patient_1 Tumor <-> CD8+ T  Tumor   CD8+ T  HLA.E    KLRK1 0.0005353052
#> 3 Patient_1 Tumor <-> CD8+ T  Tumor   CD8+ T    ADM    VIPR1 0.0002670687
#> 4 Patient_1 Tumor <-> CD8+ T  Tumor   CD8+ T  HLA.F     CD8A 0.0001695168
#> 5 Patient_1 Tumor <-> CD8+ T  Tumor   CD8+ T  HLA.A     CD8A 0.0001652873
#> 6 Patient_1 Tumor <-> CD8+ T  Tumor   CD8+ T  HLA.C     CD8A 0.0001652873
#>   PctOfTotal
#> 1  21.519765
#> 2  11.029554
#> 3   5.502746
#> 4   3.492763
#> 5   3.405619
#> 6   3.405619
plotLRContributions(lr_contrib, top_n = 10)
```

![Dot plot of the top ligand-receptor pairs contributing to Tumor and
CD8+ T communication, per patient; dot size is the share of the total
score and color the raw
contribution.](RaCiNG_files/figure-html/unnamed-chunk-19-1.png)

### Statistical analysis on the example data

The example data has no clinical annotation, so for illustration
patients are split by their CD8+ T-cell infiltration (above vs below the
median). In your own data, use a real clinical variable such as response
to therapy. (With only 10 patients many values are tied, so
[`wilcox.test()`](https://rdrr.io/r/stats/wilcox.test.html) warns that
it cannot compute exact p-values; those warnings are hidden here.)

``` r
cd8 <- skcm_example$Cmatrix[, "CD8+ T"]
grouping <- ifelse(cd8 > median(cd8), "CD8_high", "CD8_low")
wilcox_results <- wilcox_group_test(kernel_res$features, grouping)
head(wilcox_results)
#>            Comparison           Feature Wilcox_statistic Fold_change    Log2FC
#> 1 CD8_high_vs_CD8_low      Dir_B_CD8+ T               25    3.644724 1.8658095
#> 2 CD8_high_vs_CD8_low           Dir_B_M               25    1.624455 0.6999561
#> 3 CD8_high_vs_CD8_low        Dir_B_Treg               25    1.853754 0.8904495
#> 4 CD8_high_vs_CD8_low    Dir_CAF_CD8+ T               25    1.998395 0.9988417
#> 5 CD8_high_vs_CD8_low Dir_CD8+ T_CD8+ T               25    4.820017 2.2690382
#> 6 CD8_high_vs_CD8_low     Dir_CD8+ T_DC               25    3.317927 1.7302821
#>       P_value Adjusted_P_value
#> 1 0.007936508       0.03968254
#> 2 0.007936508       0.03968254
#> 3 0.007936508       0.03968254
#> 4 0.007936508       0.03968254
#> 5 0.007936508       0.03968254
#> 6 0.007936508       0.03968254
top_features_plot(wilcox_results, top_n = 10)
```

![Ranked bar chart of the direct-communication features most up- and
down-regulated in CD8-high versus CD8-low
patients.](RaCiNG_files/figure-html/unnamed-chunk-20-1.png)

For a continuous score, a random placeholder is used below; replace it
with a per-patient score of your own (e.g. an immune-response
signature).

``` r
set.seed(1)
response_score <- setNames(rnorm(nrow(kernel_res$features)), rownames(kernel_res$features))
corr_results <- correlate_features_with_score(kernel_res$features, response_score)
head(corr_results)
#>   Group      Feature        Rho   P_value Adjusted_P_value
#> 1   All    Dir_B_CAF -0.3212121 0.3676840        0.8720362
#> 2   All Dir_B_CD8+ T -0.4424242 0.2042006        0.8720362
#> 3   All   Dir_B_Endo -0.2969697 0.4069502        0.8720362
#> 4   All      Dir_B_M -0.4060606 0.2473708        0.8720362
#> 5   All     Dir_B_NK -0.5151515 0.1328231        0.8720362
#> 6   All  Dir_B_Tumor -0.3333333 0.3488462        0.8720362
correlation_plot(corr_results, top_n = 10)
```

![Diverging bar chart of the features with the strongest positive and
negative Spearman correlations with the
score.](RaCiNG_files/figure-html/unnamed-chunk-21-1.png)

### Monte Carlo method on the example data

The settings below are deliberately tiny so the example runs in seconds;
real analyses use thousands of cells (`Ncells`) and around 100 graphs
(`Ngraphs`) per patient.

``` r
set.seed(42)
mc_res <- compute_racing_montecarlo(
  input_data   = skcm_example,
  output_folder = tempdir(),
  file_name     = "skcm_example",
  nPatients     = 2,
  communication_type = c("D", "GSCC"),
  Ncells  = 100,
  Ngraphs = 5,
  Ndegree = 3
)
#> Using pre-computed input matrices; skipping input generation.
#> Running Monte Carlo simulation for 2 patients
#> Processing interaction distribution and generating CSV output...

# Direct communication: mean and standard deviation across the simulated graphs
head(mc_res$output$D$mean[, 1:5])
#>           Dir_B_B Dir_B_CAF Dir_B_CD8+ T Dir_B_DC Dir_B_Endo
#> Patient_1       0       0.6          0.0      0.0        0.2
#> Patient_2       0       0.2          1.4      1.2        0.0
head(mc_res$output$D$sd[, 1:5])
#>           Dir_B_B Dir_B_CAF Dir_B_CD8+ T Dir_B_DC Dir_B_Endo
#> Patient_1       0 1.0000000     0.000000 0.000000  0.4472136
#> Patient_2       0 0.4472136     1.224745 1.224745  0.0000000

# Fraction of all cells that belong to the GSCC, per cell type
mc_res$output$GSCC$mean
#>           GSCC_B GSCC_CAF GSCC_CD8+ T GSCC_DC GSCC_Endo GSCC_M GSCC_NK
#> Patient_1  0.002    0.014       0.000   0.004     0.008  0.010   0.000
#> Patient_2  0.012    0.014       0.068   0.070     0.002  0.088   0.002
#>           GSCC_Treg GSCC_Tumor
#> Patient_1     0.012      0.578
#> Patient_2     0.014      0.374
```

## Notes

- [`compute_racing_kernel()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/compute_racing_kernel.md)
  and
  [`compute_racing_montecarlo()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/compute_racing_montecarlo.md)
  are the main entry points.
- The kernel method computes direct communication only; GSCC is
  available with the Monte Carlo method. Wedges and triangles are no
  longer part of RaCInG.
- Both accept an `input_data` argument with pre-computed matrices (as
  returned by
  [`prepare_input_files()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/prepare_input_files.md)),
  which skips all preprocessing.
- [`prepare_input_files()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/prepare_input_files.md)
  requires the optional packages listed under Installation for
  deconvolution and prior-network assembly.
- The original Python implementation is available at
  <https://github.com/SysBioOncology/RaCInG>.
