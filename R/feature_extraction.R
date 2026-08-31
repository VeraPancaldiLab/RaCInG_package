#' Enumerate outward trust triangles
#'
#' @param Adj Adjacency matrix.
#'
#' @return A list with the triangle count and the vertex triplets found.
#' @export
Trust_Triangles <- function(Adj) {
  n <- nrow(Adj)

  # Precompute each vertex's out-neighbors once. Adj[v, ] on a sparse matrix
  # is O(nnz) per call, so looking it up fresh inside the nested loops below
  # (once per v, and again every time v appears as some other vertex's
  # neighbor) is the main cost for larger graphs - cache it instead.
  neighbors_of <- lapply(seq_len(n), function(x) which(Adj[x, ] != 0))

  per_vertex <- vector("list", n)
  NoTriangles <- 0L

  for (v in seq_len(n)) {
    neigh_v <- neighbors_of[[v]]
    rows_v <- vector("list", length(neigh_v))

    for (wi in seq_along(neigh_v)) {
      w <- neigh_v[wi]
      # Common neighbors u of v and w: v -> u, v -> w, and w -> u
      intersec <- intersect(neigh_v, neighbors_of[[w]])
      if (length(intersec) > 0) {
        rows_v[[wi]] <- cbind(intersec, v, w)
      }
    }

    rows_v <- rows_v[!vapply(rows_v, is.null, logical(1))]
    if (length(rows_v) > 0) {
      per_vertex[[v]] <- do.call(rbind, rows_v)
      NoTriangles <- NoTriangles + nrow(per_vertex[[v]])
    }
  }

  per_vertex <- per_vertex[!vapply(per_vertex, is.null, logical(1))]
  Triangle_list <- if (length(per_vertex) > 0) do.call(rbind, per_vertex) else matrix(nrow = 0, ncol = 3)
  dimnames(Triangle_list) <- NULL

  # Return:
  # - total number of triangles
  # - matrix where each row is a triangle [u, v, w] (0 rows if none found)
  return(list(NoTriangles = NoTriangles, Triangle_list = Triangle_list))
}

#' Enumerate directed cycle triangles
#'
#' @param Adj Adjacency matrix.
#'
#' @return A list with the cycle-triangle count and the vertex triplets found.
#' @export
Cycle_Triangles <- function(Adj) {
  n <- nrow(Adj)
  AdjT <- Matrix::t(Adj)

  # Precompute out- and in-neighbors once (see Trust_Triangles() for why).
  neighbors_of <- lapply(seq_len(n), function(x) which(Adj[x, ] != 0))
  backneighbors_of <- lapply(seq_len(n), function(x) which(AdjT[x, ] != 0))

  per_vertex <- vector("list", n)
  NoTriangles <- 0L

  for (v in seq_len(n)) {
    neigh_v <- neighbors_of[[v]]
    backneigh_v <- backneighbors_of[[v]]
    rows_v <- vector("list", length(neigh_v))

    for (wi in seq_along(neigh_v)) {
      w <- neigh_v[wi]
      # Nodes u with u -> v (backneigh_v) and w -> u (neighbors_of[[w]])
      # complete the cycle v -> w -> u -> v.
      intersec <- intersect(backneigh_v, neighbors_of[[w]])
      if (length(intersec) > 0) {
        rows_v[[wi]] <- cbind(intersec, v, w)
      }
    }

    rows_v <- rows_v[!vapply(rows_v, is.null, logical(1))]
    if (length(rows_v) > 0) {
      per_vertex[[v]] <- do.call(rbind, rows_v)
      NoTriangles <- NoTriangles + nrow(per_vertex[[v]])
    }
  }

  per_vertex <- per_vertex[!vapply(per_vertex, is.null, logical(1))]
  Triangle_list <- if (length(per_vertex) > 0) do.call(rbind, per_vertex) else matrix(nrow = 0, ncol = 3)
  dimnames(Triangle_list) <- NULL

  # Each triangle is counted 3 times (once from each vertex),
  # so divide by 3 to get the correct number of unique triangles
  # 1 → 2 → 3 → 1
  # v = 1: (1,2,3)
  # v = 2: (2,3,1)
  # v = 3: (3,1,2)
  return(list(NoTriangles = round(NoTriangles / 3), Triangle_list = Triangle_list))
}

#' Enumerate wedges in a directed graph
#'
#' @param Adj Adjacency matrix.
#'
#' @return A list with the wedge count and the vertex triplets found.
#' @export
Wedges <- function(Adj) {
  n <- nrow(Adj)  # Number of vertices in the graph

  # Precompute each vertex's out-neighbors once (see Trust_Triangles()).
  neighbors_of <- lapply(seq_len(n), function(x) which(Adj[x, ] != 0))

  per_vertex <- vector("list", n)
  NoWedges <- 0L

  # Loop over each vertex v
  for (v in seq_len(n)) {
    neigh_v <- neighbors_of[[v]]
    rows_v <- vector("list", length(neigh_v))

    # Loop over each neighbor w of v
    for (vi in seq_along(neigh_v)) {
      w <- neigh_v[vi]
      # Neighbors of w, excluding v to avoid trivial wedge loops
      neigh_w <- setdiff(neighbors_of[[w]], v)
      if (length(neigh_w) > 0) {
        # Each u in neigh_w forms a wedge (v -> w -> u)
        rows_v[[vi]] <- cbind(v, w, neigh_w)
      }
    }

    rows_v <- rows_v[!vapply(rows_v, is.null, logical(1))]
    if (length(rows_v) > 0) {
      per_vertex[[v]] <- do.call(rbind, rows_v)
      NoWedges <- NoWedges + nrow(per_vertex[[v]])
    }
  }

  per_vertex <- per_vertex[!vapply(per_vertex, is.null, logical(1))]
  Wedge_list <- if (length(per_vertex) > 0) do.call(rbind, per_vertex) else matrix(nrow = 0, ncol = 3)
  dimnames(Wedge_list) <- NULL

  # Return the total count and the full list as a matrix (0 rows if none found)
  return(list(NoWedges = NoWedges, Wedge_list = Wedge_list))
}

#' Find strongly connected components with Tarjan's algorithm
#'
#' @param Adj Adjacency matrix.
#'
#' @return A list of strongly connected components.
#' @keywords internal
TarjanIterative <- function(Adj) {
  N <- nrow(Adj)
  disc <- rep(NA_integer_, N)     # discovery index
  lowlink <- rep(NA_integer_, N)
  onstack <- rep(FALSE, N)
  stack <- integer(0)             # vertices currently on the SCC stack
  groups <- list()
  timer <- 1L

  # Cache each vertex's out-neighbors once, since Adj[v, ] lookups are
  # repeated every time we resume a suspended call frame for v.
  neighbors_of <- lapply(seq_len(N), function(v) which(Adj[v, ] != 0))

  for (start in seq_len(N)) {
    if (!is.na(disc[start])) next

    # Explicit call stack emulating recursive strongconnect(start). Each
    # frame tracks its node and how many of its neighbors have been visited
    # so far; the position is advanced *before* recursing into a neighbor,
    # so resuming a frame naturally continues from the next neighbor.
    disc[start] <- timer; lowlink[start] <- timer; timer <- timer + 1L
    stack <- c(stack, start); onstack[start] <- TRUE
    call_stack <- list(list(node = start, pos = 1L))

    while (length(call_stack) > 0) {
      frame <- call_stack[[length(call_stack)]]
      v <- frame$node
      neigh <- neighbors_of[[v]]

      if (frame$pos <= length(neigh)) {
        w <- neigh[frame$pos]
        call_stack[[length(call_stack)]]$pos <- frame$pos + 1L

        if (is.na(disc[w])) {
          disc[w] <- timer; lowlink[w] <- timer; timer <- timer + 1L
          stack <- c(stack, w); onstack[w] <- TRUE
          call_stack[[length(call_stack) + 1]] <- list(node = w, pos = 1L)
        } else if (onstack[w]) {
          lowlink[v] <- min(lowlink[v], disc[w])
        }
      } else {
        # All neighbors of v visited: pop v's frame and propagate its
        # lowlink to the parent frame (if any) exactly once.
        call_stack[[length(call_stack)]] <- NULL
        if (length(call_stack) > 0) {
          parent <- call_stack[[length(call_stack)]]$node
          lowlink[parent] <- min(lowlink[parent], lowlink[v])
        }

        if (lowlink[v] == disc[v]) {
          com <- integer(0)
          repeat {
            w <- stack[length(stack)]
            stack <- stack[-length(stack)]
            onstack[w] <- FALSE
            com <- c(com, w)
            if (w == v) break
          }
          groups[[length(groups) + 1]] <- com
        }
      }
    }
  }

  return(groups)
}

#' Extract the giant strongly connected component
#'
#' @param Adj Adjacency matrix.
#'
#' @return Integer vector of vertex indices in the GSCC.
#' @export
GSCC <- function(Adj) {
  SCCs <- TarjanIterative(Adj)
  sizes <- sapply(SCCs, length)
  giant_idx <- which.max(sizes)
  return(SCCs[[giant_idx]])
}