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

#' Count graph motifs by cell-type combination
#'
#' @param oblist Matrix of graph objects where each row contains vertex indices.
#' @param V Integer vector mapping vertices to cell types.
#' @param maxTypes Optional maximum number of types used to size the output array.
#'
#' @return An array counting how often each type combination occurs.
#' @export
Count_Types <- function(oblist, V, maxTypes = 0){

  # oblist: matrix of objects (each row = one object, e.g., a wedge or triangle)
  # V: vector mapping vertices to their cell types (V[vertex] = type)
  # maxTypes: optional, maximum number of cell types
  # Returns: multi-dimensional array that counts how many times each type combination occurs
  
  # Determine number of types to use for the array dimensions
  if (maxTypes > 0) {
    dim_val <- maxTypes
  } else {
    dim_val <- max(V)   # use the maximum type present in V
  }
  
  # Create an n-dimensional array to store counts
  # - ncol(oblist) = number of vertices per object (e.g., 3 for wedges)
  # - Each dimension = number of possible cell types
  counttensor <- array(0, dim = rep(dim_val, ncol(oblist)))
  
  # Replace vertex indices in each object with their corresponding types
  # Each row now shows the types for the vertices of that object
  typelist <- matrix(V[oblist], nrow = nrow(oblist), ncol = ncol(oblist))
  
  # Loop over all objects (seq_len avoids the 1:0 == c(1,0) pitfall when empty)
  for (i in seq_len(nrow(typelist))) {
    types <- typelist[i, ]  # types of the vertices for this object
    
    # Increment the count in the array at the position matching this type combination
    # - Example: if types = [1,3,2], then counttensor[1,3,2] increases by 1
    # - This way, we track how many times each combination of types occurs
    counttensor[matrix(types, nrow = 1)] <- counttensor[matrix(types, nrow = 1)] + 1
  }
  
  # Return the array with counts for all type combinations
  return(counttensor)
}

#' Poisson branching-process helper
#'
#' @param x Current estimate vector.
#' @param M Mean offspring matrix.
#' @param sens Length of `x`.
#'
#' @return A numeric vector used by `nleqslv::nleqslv()`.
#' @keywords internal
poiBPFunc <- function(x, M, sens) {
  # x: vector of current GSCC estimates
  # M: matrix of average group sizes
  # sens: size of x
  # Returns: vector of survival probabilities
  inter <- matrix(rep(x, each = sens), nrow = sens)
  result <- inter * M
  return(rep(1, sens) - x - exp(-rowSums(result)))
}

