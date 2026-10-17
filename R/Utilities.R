#' Convert an edge list to a sparse adjacency matrix
#'
#' @param E Two-column edge list with source and target vertex indices.
#' @param N Total number of vertices.
#'
#' @return A sparse adjacency matrix of class `dgCMatrix`.
#' @export
EdgetoAdj <- function(E, N) {
  # E: a two-column matrix where each row represents an edge (from, to)
  #    e.g., E[i,] = c(2,5) means an edge from vertex 2 to vertex 5
  # N: the total number of vertices in the graph
  # Returns: a sparse adjacency matrix (class "dgCMatrix") representing the graph
  # Note: If there are duplicate edges between the same vertices, the matrix entries
  #       will contain the sum of these edges (i.e., the number of parallel edges)
  
  # Construct the sparse adjacency matrix:
  # - i = row indices (source vertices)
  # - j = column indices (target vertices)
  # - x = 1 for each edge (will be summed automatically if duplicates exist)
  # - dims = dimensions of the adjacency matrix (N x N)
  Adj <- Matrix::sparseMatrix(
    i = E[,1], 
    j = E[,2], 
    x = rep(1, nrow(E)), 
    dims = c(N, N)
  )
  
  return(Adj)
}

#' Convert an edge list to a sparse adjacency matrix without self-loops
#'
#' @param E Two-column edge list with source and target vertex indices.
#' @param N Total number of vertices.
#'
#' @return A sparse adjacency matrix with diagonal entries removed.
#' @export
EdgetoAdj_No_loop <- function(E, N) {
  # E: two-column matrix where each row represents an edge (from, to)
  # N: total number of vertices
  # Returns: sparse adjacency matrix without self-loops (edges from a vertex to itself)
  
  # Step 1: Identify edges that are NOT self-loops
  # Keep only edges where the source and target are different
  keep <- E[,1] != E[,2]
  
  # Step 2: Filter the edge list to remove self-loops
  E_noloop <- E[keep,, drop = FALSE]
  
  # Step 3: Construct the sparse adjacency matrix
  # - i = source vertices
  # - j = target vertices
  # - x = 1 for each edge
  # - dims = N x N
  # This will sum duplicate edges automatically, if any
  Adj <- Matrix::sparseMatrix(
    i = E_noloop[,1], 
    j = E_noloop[,2], 
    x = rep(1, nrow(E_noloop)), 
    dims = c(N, N)
  )
  
  return(Adj)
}
