#' Generate a random cell-type distribution
#'
#' @param cellTypeNo Number of cell types.
#'
#' @return A probability vector summing to 1.
#' @export
genRandomCellTypeDistr <- function(cellTypeNo) {
  # Generate a random probability distribution for cell types
  #
  # Arguments:
  #   cellTypeNo : number of cell types
  #
  # Returns:
  #   A numeric vector of length `cellTypeNo` representing the probability
  #   of each cell type. The sum of the probabilities is 1.
  
  probs <- stats::runif(cellTypeNo)  # generate random values for each cell type
  probs <- probs / sum(probs)  # normalize so total probability sums to 1
  
  return(probs)
}

#' Generate a random ligand-receptor distribution
#'
#' @param ligNo Number of ligand types.
#' @param recNo Number of receptor types.
#'
#' @return A normalized ligand-by-receptor probability matrix.
#' @export
genRandomLigRecDistr <- function(ligNo, recNo) {
  # Generate a random ligand-receptor probability distribution matrix
  #
  # Arguments:
  #   ligNo : number of ligand types (rows)
  #   recNo : number of receptor types (columns)
  #
  # Returns:
  #   A ligNo x recNo matrix where each element represents the probability
  #   of a specific ligand-receptor pair occurring. The sum of all elements is 1.
  
  mat <- matrix(
    stats::runif(ligNo * recNo),   # generate random values for each ligand-receptor pair
    nrow = ligNo,
    ncol = recNo
  )
  
  mat <- mat / sum(mat)      # normalize the matrix so that all probabilities sum to 1
  
  return(mat)
}

#' Generate a random cell-to-ligand compatibility matrix
#'
#' @param cellTypeNo Number of cell types.
#' @param ligNo Number of ligand types.
#'
#' @return A binary matrix describing ligand availability per cell type.
#' @export
genRandomCellLigands <- function(cellTypeNo, ligNo) {
  # Generate a random 0/1 matrix indicating which cell types can secrete which ligands
  # Arguments:
  #   cellTypeNo : number of cell types (rows)
  #   ligNo      : number of ligand types (columns)
  # Returns:
  #   A matrix of size cellTypeNo x ligNo
  #     - 1 means the cell type can secrete that ligand
  #     - 0 means the cell type cannot secrete that ligand
  
  mat <- matrix(
    sample(0:1, cellTypeNo * ligNo, replace = TRUE),  # randomly assign 0 or 1 for each cell-ligand pair
    nrow = cellTypeNo                                   # set number of rows = number of cell types
  )
  
  return(mat)
}

#' Generate a random cell-to-receptor compatibility matrix
#'
#' @param cellTypeNo Number of cell types.
#' @param recNo Number of receptor types.
#'
#' @return A binary matrix describing receptor availability per cell type.
#' @export
genRandomCellReceptors <- function(cellTypeNo, recNo) {
  # Generate a random 0/1 matrix indicating which cell types can express which receptors
  # Arguments:
  #   cellTypeNo : number of cell types (rows)
  #   recNo      : number of receptor types (columns)
  # Returns:
  #   A matrix of size cellTypeNo x recNo
  #     - 1 means the cell type can express that receptor
  #     - 0 means the cell type cannot express that receptor
  
  mat <- matrix(
    sample(0:1, cellTypeNo * recNo, replace = TRUE),  # randomly sample 0 or 1 for each cell-receptor pair
    nrow = cellTypeNo                                   # set number of rows = number of cell types
  )
  
  return(mat)
}

#' Sample cell-type labels for graph vertices
#'
#' @param Dcelltype Probability vector over cell types.
#' @param vertexNo Number of vertices to sample.
#'
#' @return An integer vector of sampled cell-type indices.
#' @export
genRandomCellTypeList <- function(Dcelltype, vertexNo) {
  # This function generates a random list of cell types for the graph vertices.
  # Each vertex is assigned a cell type based on the probabilities in Dcelltype.
  #
  # Parameters:
  #   Dcelltype : numeric vector
  #       The probability distribution of each cell type (should sum to 1).
  #   vertexNo : integer
  #       The number of vertices (cells) to generate.
  #
  # Returns:
  #   A numeric vector of length vertexNo where each entry represents
  #   the assigned cell type index (1-based) for that vertex.
  
  sample(
    seq_along(Dcelltype),  # Generates a sequence 1:length(Dcelltype)
    vertexNo,               # Number of samples to draw (number of vertices)
    replace = TRUE,         # Sampling with replacement so types can repeat
    prob = Dcelltype        # Weighted probabilities according to Dcelltype
  )
}

#' Sample an edge list from ligand-receptor probabilities
#'
#' @param Dligrec Ligand-by-receptor probability matrix.
#' @param vertextypelist Integer vector assigning a cell type to each vertex.
#' @param structurelig Cell-by-ligand compatibility matrix.
#' @param structurerec Cell-by-receptor compatibility matrix.
#' @param edgeNo Number of edges to generate (typically `round(avdeg * N)`).
#'
#' @return A list with `cell_connection` and `ligrec_type` matrices.
#' @keywords internal
genRandomEdgeList <- function(Dligrec, vertextypelist, structurelig, structurerec, edgeNo) {

  ligNo <- nrow(Dligrec)   # Number of ligand types
  recNo <- ncol(Dligrec)   # Number of receptor types

  # Flatten the ligand-receptor probability matrix to a vector
  distr <- as.vector(Dligrec)
  # Randomly sample edgeNo edges based on probabilities
  linearEdgeList <- sample(seq_along(distr), edgeNo, prob = distr, replace = TRUE)
  
  # Convert linear indices to row (ligand) and column (receptor) indices
  edge_indices <- arrayInd(linearEdgeList, .dim = dim(Dligrec))
  lig_indices <- edge_indices[, 1]  # Ligand indices (rows of Dligrec)
  rec_indices <- edge_indices[, 2]  # Receptor indices (columns of Dligrec)
  
  # Store the types of ligands and receptors for each edge
  edgetypelist <- cbind(lig_indices, rec_indices)
  # Initialize the final edge list with ligand and receptor positions
  edgelist <- edgetypelist
  
  # Assign ligands to compatible cells
  for (i in 1:ligNo) {
    count <- sum(edgetypelist[,1] == i)                  # Number of edges with ligand i
    acceptingCells <- which(structurelig[, i] == 1)      # Cell types that can produce ligand i
    connectionChoice <- which(vertextypelist %in% acceptingCells)  # Vertices of compatible cell types
    if (length(connectionChoice) > 0) {
      # Randomly assign these edges to compatible cells
      temp <- sample(connectionChoice, count, replace = TRUE)
    } else {
      # If no compatible cells exist, mark edges as invalid (-1)
      temp <- rep(-1, count)
    }
    # Update the first column of edgelist (ligand endpoints)
    edgelist[edgetypelist[,1] == i, 1] <- temp
  }
  
  # Assign receptors to compatible cells (similar procedure)
  for (i in 1:recNo) {
    count <- sum(edgetypelist[,2] == i)                  # Number of edges with receptor i
    acceptingCells <- which(structurerec[, i] == 1)      # Cell types that can accept receptor i
    connectionChoice <- which(vertextypelist %in% acceptingCells)  # Compatible vertices
    if (length(connectionChoice) > 0) {
      temp <- sample(connectionChoice, count, replace = TRUE)
    } else {
      temp <- rep(-1, count)                             # Mark invalid edges
    }
    # Update the second column of edgelist (receptor endpoints)
    edgelist[edgetypelist[,2] == i, 2] <- temp
  }
  
  # Remove edges where either ligand or receptor assignment failed (-1)
  keep <- !(edgelist[,1] == -1 | edgelist[,2] == -1)
  edgelist <- edgelist[keep, , drop = FALSE]
  edgetypelist <- edgetypelist[keep, , drop = FALSE]
  
  # Return the final edge list and the corresponding ligand-receptor types
  return(list(cell_connection = edgelist, ligrec_type = edgetypelist))
}

#' Generate a single RaCInG graph realization
#'
#' @param N Number of vertices (cells) in the graph.
#' @param avdeg Target average degree.
#' @param cellLigList Cell-by-ligand compatibility matrix.
#' @param cellRecList Cell-by-receptor compatibility matrix.
#' @param Dcelltype Cell-type abundance probabilities.
#' @param Dligrec Ligand-by-receptor probability matrix.
#' @param Signmatrix Optional ligand-receptor sign matrix.
#' @param genRandom Logical; if `TRUE`, generate random test inputs internally.
#'
#' @return A list with vertex labels, an edge list, and ligand-receptor types.
#' @export
model1 <- function(N, avdeg,
                   cellLigList = NULL, cellRecList = NULL,
                   Dcelltype = NULL, Dligrec = NULL, Signmatrix = NULL,
                   genRandom = TRUE) {
  
  if (genRandom) {
    cellTypeNo <- 11
    ligNo <- 10
    recNo <- 10
    M <- round(avdeg * N)
    Dcelltype <- genRandomCellTypeDistr(cellTypeNo)
    Dligrec <- genRandomLigRecDistr(ligNo, recNo)
    cellLigList <- genRandomCellLigands(cellTypeNo, ligNo)
    cellRecList <- genRandomCellReceptors(cellTypeNo, recNo)
  } else {
    cellTypeNo <- length(Dcelltype)
    ligNo <- nrow(Dligrec)
    recNo <- ncol(Dligrec)
    M <- round(avdeg * N)
  }
  
  # -------------------------------
  # Generate vertex list
  # -------------------------------
  # V is a vector of length N containing the **cell type of each vertex**.
  # Example: V[8] = 9 means vertex/cell #8 has type 9.
  V <- genRandomCellTypeList(Dcelltype, N)
  
  # -------------------------------
  # Generate edge list
  # -------------------------------
  # E is a matrix with two columns: ligand vertex index and receptor vertex index.
  # These are **indices of cells**, NOT their types.
  # For example, E[1,] = c(8, 11) means there is an edge from cell #8 to cell #11.
  # To get the types of these cells, you can use V[E[,1]] and V[E[,2]].
  edges <- genRandomEdgeList(Dligrec, V, cellLigList, cellRecList, edgeNo = M)
  E <- edges$cell_connection
  
  # -------------------------------
  # Ligand-receptor type list
  # -------------------------------
  # types is a matrix of the **ligand and receptor types** corresponding to each edge.
  # types[i,1] = ligand type, types[i,2] = receptor type for edge i.
  types <- edges$ligrec_type
  
  # Add sign info if provided
  if (!is.null(Signmatrix)) {
    interactions <- mapply(function(lig, rec) Signmatrix[lig, rec], types[,1], types[,2])
    types <- cbind(types, interactions)
  }
  
  return(list(V = V, E = E, types = types))
}

#' Generate a graph under a uniformized ligand-receptor baseline
#'
#' @param LRdistr Ligand-receptor tensor.
#' @param Lmatrix Cell-by-ligand compatibility matrix.
#' @param Rmatrix Cell-by-receptor compatibility matrix.
#' @param Cdistr Patient-by-cell-type abundance matrix.
#' @param cellTypes Character vector of cell-type labels.
#' @param patient Patient index to simulate.
#' @param N Number of cells in the generated graph.
#' @param avdeg Target average degree.
#'
#' @return A list containing the simulated graph and the uniform LR distribution used.
#' @export
generateUniformLRGraph <- function(LRdistr, Lmatrix, Rmatrix, Cdistr, cellTypes, patient = 1, N = 20, avdeg = 2) {
  # Set seed for reproducibility
  set.seed(1)
  
  # Create a new array to hold the uniform LR distribution
  # It has the same dimensions as the original LR distribution array
  LRdistrUniform <- array(0, dim = dim(LRdistr))
  
  # Calculate the value to assign to every non-zero entry
  # We want the non-zero probabilities to be uniform and sum to 1
  # sum(LRdistr[,,patient] != 0) counts how many non-zero entries there are
  normval <- 1 / sum(LRdistr[,,patient] != 0)
  
  # Copy the LR distribution for the selected patient (slice of the 3D array)
  copy <- LRdistr[,,patient]
  
  # Replace all non-zero entries with the uniform value
  # Zero entries remain zero
  copy[copy > 0] <- normval
  
  # Store the uniformized matrix back in the uniform LR array
  # Only for the selected patient
  LRdistrUniform[,,patient] <- copy
  
  # Generate graph with uniform LR distribution
  result <- model1(
    N = N,
    avdeg = avdeg,
    cellLigList = Lmatrix,
    cellRecList = Rmatrix,
    Dcelltype = Cdistr[patient, ],  # select patient distribution
    Dligrec = LRdistrUniform[,,patient],
    genRandom = FALSE
  )
  
  Vnorm <- result$V
  Enorm <- result$E
  typesNorm <- result$types
  
  return(list(
    V = Vnorm,
    E = Enorm,
    types = typesNorm,
    LRdistrUniform = LRdistrUniform[,,patient]
  ))
}
