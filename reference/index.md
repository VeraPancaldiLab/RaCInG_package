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
- [`calculateWedges()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/calculateWedges.md)
  : Calculate wedge features from a kernel
- [`computeTriangles()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/computeTriangles.md)
  : Compute triangle features from kernel matrices
- [`computeGSCC()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/computeGSCC.md)
  : Compute GSCC features from kernel matrices
- [`countWedges()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/countWedges.md)
  : Count wedges across Monte Carlo graph simulations
- [`countTrustTriangles()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/countTrustTriangles.md)
  : Count trust triangles across Monte Carlo graph simulations
- [`countCycleTriangles()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/countCycleTriangles.md)
  : Count cycle triangles across Monte Carlo graph simulations
- [`countDirect()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/countDirect.md)
  : Count direct edges across Monte Carlo graph simulations
- [`countGSCC()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/countGSCC.md)
  : Count GSCC contributions across Monte Carlo graph simulations
- [`model1()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/model1.md)
  : Generate a single RaCInG graph realization
- [`runSim()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/runSim.md)
  : Run Monte Carlo simulations for one or more patients
- [`getGSCCAnalytically()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/getGSCCAnalytically.md)
  : Legacy GSCC helper

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
- [`genRandomCellTypeDistr()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/genRandomCellTypeDistr.md)
  : Generate a random cell-type distribution
- [`genRandomLigRecDistr()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/genRandomLigRecDistr.md)
  : Generate a random ligand-receptor distribution
- [`genRandomCellLigands()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/genRandomCellLigands.md)
  : Generate a random cell-to-ligand compatibility matrix
- [`genRandomCellReceptors()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/genRandomCellReceptors.md)
  : Generate a random cell-to-receptor compatibility matrix
- [`genRandomCellTypeList()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/genRandomCellTypeList.md)
  : Sample cell-type labels for graph vertices
- [`generateUniformLRGraph()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/generateUniformLRGraph.md)
  : Generate a graph under a uniformized ligand-receptor baseline

## Graph utilities

- [`EdgetoAdj()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/EdgetoAdj.md)
  : Convert an edge list to a sparse adjacency matrix
- [`EdgetoAdj_No_loop()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/EdgetoAdj_No_loop.md)
  : Convert an edge list to a sparse adjacency matrix without self-loops
- [`Count_Types()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/Count_Types.md)
  : Count graph motifs by cell-type combination
- [`Trust_Triangles()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/Trust_Triangles.md)
  : Enumerate outward trust triangles
- [`Cycle_Triangles()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/Cycle_Triangles.md)
  : Enumerate directed cycle triangles
- [`Wedges()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/Wedges.md)
  : Enumerate wedges in a directed graph
- [`Find_Number_Trust_Triangles_Unique()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/Find_Number_Trust_Triangles_Unique.md)
  : Count unique trust triangles in a directed graph
- [`Find_Number_Triangles()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/Find_Number_Triangles.md)
  : Count triangles allowing multi-edges
- [`Find_Number_Triangles_Unique()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/Find_Number_Triangles_Unique.md)
  : Count unique triangles in a directed graph
- [`Find_Number_2Loops()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/Find_Number_2Loops.md)
  : Count reciprocal 2-loops allowing multi-edges
- [`Find_Number_2Loops_Unique()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/Find_Number_2Loops_Unique.md)
  : Count unique reciprocal 2-loops
- [`Find_Number_Wedges()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/Find_Number_Wedges.md)
  : Count wedges allowing multi-edges
- [`Find_Number_Wedges_Unique()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/Find_Number_Wedges_Unique.md)
  : Count unique wedges in a directed graph
- [`GSCC()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/GSCC.md)
  : Extract the giant strongly connected component
- [`IN()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/IN.md)
  : Compute the IN component of a directed graph
- [`OUT()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/OUT.md)
  : Compute the OUT component of a directed graph

## Example data

- [`skcm_example`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/skcm_example.md)
  : SKCM melanoma example input for RaCInG

## Statistics and I/O

- [`wilcox_group_test()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/wilcox_group_test.md)
  : Run Wilcoxon tests across network features
- [`volcano_plot()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/volcano_plot.md)
  : Create a volcano plot from Wilcoxon results
- [`correlate_features_with_score()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/correlate_features_with_score.md)
  : Correlate network features with an external score
- [`correlation_plot()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/correlation_plot.md)
  : Rainfall plot of feature-score correlations
- [`Read_Sim_Output()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/Read_Sim_Output.md)
  : Read a RaCInG simulation output file
- [`compute_results_processing()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/compute_results_processing.md)
  : Convert raw simulation outputs into feature matrices
