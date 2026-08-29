# Package index

## Main workflows

- [`prepare_input_files()`](https://mhurtado13.github.io/racing/reference/prepare_input_files.md)
  : Build RaCInG input files from raw count data
- [`compute_racing_kernel()`](https://mhurtado13.github.io/racing/reference/compute_racing_kernel.md)
  : Run the full kernel-based RaCInG workflow
- [`compute_racing_montecarlo()`](https://mhurtado13.github.io/racing/reference/compute_racing_montecarlo.md)
  : Run the full Monte Carlo RaCInG workflow

## Kernel and Monte Carlo methods

- [`compute_kernel()`](https://mhurtado13.github.io/racing/reference/compute_kernel.md)
  : Compute the RaCInG kernel for one or more patients
- [`compute_kernel_features()`](https://mhurtado13.github.io/racing/reference/compute_kernel_features.md)
  : Derive communication features from a kernel
- [`calculateDirect()`](https://mhurtado13.github.io/racing/reference/calculateDirect.md)
  : Calculate direct communication features from a kernel
- [`calculateWedges()`](https://mhurtado13.github.io/racing/reference/calculateWedges.md)
  : Calculate wedge features from a kernel
- [`computeTriangles()`](https://mhurtado13.github.io/racing/reference/computeTriangles.md)
  : Compute triangle features from kernel matrices
- [`computeGSCC()`](https://mhurtado13.github.io/racing/reference/computeGSCC.md)
  : Compute GSCC features from kernel matrices
- [`countWedges()`](https://mhurtado13.github.io/racing/reference/countWedges.md)
  : Count wedges across Monte Carlo graph simulations
- [`countTrustTriangles()`](https://mhurtado13.github.io/racing/reference/countTrustTriangles.md)
  : Count trust triangles across Monte Carlo graph simulations
- [`countCycleTriangles()`](https://mhurtado13.github.io/racing/reference/countCycleTriangles.md)
  : Count cycle triangles across Monte Carlo graph simulations
- [`countDirect()`](https://mhurtado13.github.io/racing/reference/countDirect.md)
  : Count direct edges across Monte Carlo graph simulations
- [`countGSCC()`](https://mhurtado13.github.io/racing/reference/countGSCC.md)
  : Count GSCC contributions across Monte Carlo graph simulations
- [`model1()`](https://mhurtado13.github.io/racing/reference/model1.md)
  : Generate a single RaCInG graph realization
- [`runSim()`](https://mhurtado13.github.io/racing/reference/runSim.md)
  : Run Monte Carlo simulations for one or more patients
- [`getGSCCAnalytically()`](https://mhurtado13.github.io/racing/reference/getGSCCAnalytically.md)
  : Legacy GSCC helper

## Input generation

- [`createCellLigList()`](https://mhurtado13.github.io/racing/reference/createCellLigList.md)
  : Read a cell-to-ligand compatibility matrix
- [`createCellRecList()`](https://mhurtado13.github.io/racing/reference/createCellRecList.md)
  : Read a cell-to-receptor compatibility matrix
- [`createCellTypeDistr()`](https://mhurtado13.github.io/racing/reference/createCellTypeDistr.md)
  : Read and normalize cell-type abundance estimates
- [`createInteractionDistr()`](https://mhurtado13.github.io/racing/reference/createInteractionDistr.md)
  : Read ligand-receptor interaction probabilities
- [`Read_Lig_Rec_Interaction()`](https://mhurtado13.github.io/racing/reference/Read_Lig_Rec_Interaction.md)
  : Read a ligand-receptor sign matrix
- [`genRandomCellTypeDistr()`](https://mhurtado13.github.io/racing/reference/genRandomCellTypeDistr.md)
  : Generate a random cell-type distribution
- [`genRandomLigRecDistr()`](https://mhurtado13.github.io/racing/reference/genRandomLigRecDistr.md)
  : Generate a random ligand-receptor distribution
- [`genRandomCellLigands()`](https://mhurtado13.github.io/racing/reference/genRandomCellLigands.md)
  : Generate a random cell-to-ligand compatibility matrix
- [`genRandomCellReceptors()`](https://mhurtado13.github.io/racing/reference/genRandomCellReceptors.md)
  : Generate a random cell-to-receptor compatibility matrix
- [`genRandomCellTypeList()`](https://mhurtado13.github.io/racing/reference/genRandomCellTypeList.md)
  : Sample cell-type labels for graph vertices
- [`generateUniformLRGraph()`](https://mhurtado13.github.io/racing/reference/generateUniformLRGraph.md)
  : Generate a graph under a uniformized ligand-receptor baseline

## Graph utilities

- [`EdgetoAdj()`](https://mhurtado13.github.io/racing/reference/EdgetoAdj.md)
  : Convert an edge list to a sparse adjacency matrix
- [`EdgetoAdj_No_loop()`](https://mhurtado13.github.io/racing/reference/EdgetoAdj_No_loop.md)
  : Convert an edge list to a sparse adjacency matrix without self-loops
- [`Count_Types()`](https://mhurtado13.github.io/racing/reference/Count_Types.md)
  : Count graph motifs by cell-type combination
- [`Trust_Triangles()`](https://mhurtado13.github.io/racing/reference/Trust_Triangles.md)
  : Enumerate outward trust triangles
- [`Cycle_Triangles()`](https://mhurtado13.github.io/racing/reference/Cycle_Triangles.md)
  : Enumerate directed cycle triangles
- [`Wedges()`](https://mhurtado13.github.io/racing/reference/Wedges.md)
  : Enumerate wedges in a directed graph
- [`Find_Number_Trust_Triangles_Unique()`](https://mhurtado13.github.io/racing/reference/Find_Number_Trust_Triangles_Unique.md)
  : Count unique trust triangles in a directed graph
- [`Find_Number_Triangles()`](https://mhurtado13.github.io/racing/reference/Find_Number_Triangles.md)
  : Count triangles allowing multi-edges
- [`Find_Number_Triangles_Unique()`](https://mhurtado13.github.io/racing/reference/Find_Number_Triangles_Unique.md)
  : Count unique triangles in a directed graph
- [`Find_Number_2Loops()`](https://mhurtado13.github.io/racing/reference/Find_Number_2Loops.md)
  : Count reciprocal 2-loops allowing multi-edges
- [`Find_Number_2Loops_Unique()`](https://mhurtado13.github.io/racing/reference/Find_Number_2Loops_Unique.md)
  : Count unique reciprocal 2-loops
- [`Find_Number_Wedges()`](https://mhurtado13.github.io/racing/reference/Find_Number_Wedges.md)
  : Count wedges allowing multi-edges
- [`Find_Number_Wedges_Unique()`](https://mhurtado13.github.io/racing/reference/Find_Number_Wedges_Unique.md)
  : Count unique wedges in a directed graph
- [`GSCC()`](https://mhurtado13.github.io/racing/reference/GSCC.md) :
  Extract the giant strongly connected component
- [`IN()`](https://mhurtado13.github.io/racing/reference/IN.md) :
  Compute the IN component of a directed graph
- [`OUT()`](https://mhurtado13.github.io/racing/reference/OUT.md) :
  Compute the OUT component of a directed graph

## Example data

- [`skcm_example`](https://mhurtado13.github.io/racing/reference/skcm_example.md)
  : SKCM melanoma example input for RaCInG

## Statistics and I/O

- [`wilcox_group_test()`](https://mhurtado13.github.io/racing/reference/wilcox_group_test.md)
  : Run Wilcoxon tests across network features
- [`volcano_plot()`](https://mhurtado13.github.io/racing/reference/volcano_plot.md)
  : Create a volcano plot from Wilcoxon results
- [`correlate_features_with_score()`](https://mhurtado13.github.io/racing/reference/correlate_features_with_score.md)
  : Correlate network features with an external score
- [`correlation_plot()`](https://mhurtado13.github.io/racing/reference/correlation_plot.md)
  : Rainfall plot of feature-score correlations
- [`Read_Sim_Output()`](https://mhurtado13.github.io/racing/reference/Read_Sim_Output.md)
  : Read a RaCInG simulation output file
- [`compute_results_processing()`](https://mhurtado13.github.io/racing/reference/compute_results_processing.md)
  : Convert raw simulation outputs into feature matrices
