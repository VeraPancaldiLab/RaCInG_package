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
#'
#' @return A list with vertex labels, an edge list, and ligand-receptor types.
#' @export
model1 <- function(N, avdeg, cellLigList, cellRecList, Dcelltype, Dligrec, Signmatrix = NULL) {

  M <- round(avdeg * N)

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
