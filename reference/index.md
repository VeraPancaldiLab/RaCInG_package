# Package index

## Main workflows

- [`prepare_input_files()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/prepare_input_files.md)
  : Build RaCInG input files from raw count data
- [`compute_racing_kernel()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/compute_racing_kernel.md)
  : Run the full kernel-based RaCInG workflow
- [`compute_racing_montecarlo()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/compute_racing_montecarlo.md)
  : Run the full Monte Carlo RaCInG workflow

## Kernel and Monte Carlo methods

- [`compute_kernel()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/compute_kernel.md)
  : Compute the RaCInG kernel for one or more patients
- [`compute_kernel_features()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/compute_kernel_features.md)
  : Derive communication features from a kernel
- [`calculateDirect()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/calculateDirect.md)
  : Calculate direct communication features from a kernel
- [`countDirect()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/countDirect.md)
  : Count direct edges across Monte Carlo graph simulations
- [`countGSCC()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/countGSCC.md)
  : Count GSCC contributions across Monte Carlo graph simulations
- [`countAllTypes()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/countAllTypes.md)
  : Extract multiple communication feature types from the same simulated
  graphs
- [`model1()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/model1.md)
  : Generate a single RaCInG graph realization
- [`runSim()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/runSim.md)
  : Run Monte Carlo simulations for one or more patients

## Ligand-receptor contribution analysis

- [`computeLRContributions()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/computeLRContributions.md)
  : Decompose a Direct communication score into per-ligand-receptor-pair
  contributions
- [`plotLRContributions()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/plotLRContributions.md)
  : Dot plot of top ligand-receptor pair contributions (LIANA-style)

## Input generation

- [`createCellLigList()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/createCellLigList.md)
  : Read a cell-to-ligand compatibility matrix
- [`createCellRecList()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/createCellRecList.md)
  : Read a cell-to-receptor compatibility matrix
- [`createCellTypeDistr()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/createCellTypeDistr.md)
  : Read and normalize cell-type abundance estimates
- [`createInteractionDistr()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/createInteractionDistr.md)
  : Read ligand-receptor interaction probabilities
- [`Read_Lig_Rec_Interaction()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/Read_Lig_Rec_Interaction.md)
  : Read a ligand-receptor sign matrix
- [`genRandomCellTypeList()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/genRandomCellTypeList.md)
  : Sample cell-type labels for graph vertices

## Graph utilities

- [`EdgetoAdj()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/EdgetoAdj.md)
  : Convert an edge list to a sparse adjacency matrix
- [`EdgetoAdj_No_loop()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/EdgetoAdj_No_loop.md)
  : Convert an edge list to a sparse adjacency matrix without self-loops
- [`GSCC()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/GSCC.md)
  : Extract the giant strongly connected component

## Example data

- [`skcm_example`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/skcm_example.md)
  : SKCM melanoma example input for RaCInG

## Statistics and I/O

- [`wilcox_group_test()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/wilcox_group_test.md)
  : Run Wilcoxon tests across network features
- [`top_features_plot()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/top_features_plot.md)
  : Ranked bar plot of top up/down features from Wilcoxon results
- [`correlate_features_with_score()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/correlate_features_with_score.md)
  : Correlate network features with an external score
- [`correlation_plot()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/correlation_plot.md)
  : Rainfall plot of feature-score correlations
- [`Read_Sim_Output()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/Read_Sim_Output.md)
  : Read a RaCInG simulation output file
- [`compute_results_processing()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/compute_results_processing.md)
  : Convert raw simulation outputs into feature matrices
