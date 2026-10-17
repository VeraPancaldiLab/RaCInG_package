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