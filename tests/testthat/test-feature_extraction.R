sp <- function(i, j, n) Matrix::sparseMatrix(i = i, j = j, x = rep(1, length(i)), dims = c(n, n))

test_that("TarjanIterative()/GSCC() find strongly connected components correctly", {
  # Regression test: the previous iterative implementation used `stack <- ...`
  # (instead of `<<-`) inside a closure, silently desyncing its SCC stack from
  # the tracked state, and `seq(i + 1, length(neigh))` produced a descending
  # sequence for zero-out-degree vertices - both crashed or corrupted output
  # on ordinary graphs.

  # A 3-cycle: one SCC {1, 2, 3}
  adj <- sp(c(1, 2, 3), c(2, 3, 1), 3)
  sccs <- lapply(TarjanIterative(adj), sort)
  expect_equal(sccs[[1]], c(1L, 2L, 3L))

  # A DAG (no cycles): every vertex is its own SCC
  adj_dag <- sp(c(1, 2), c(2, 3), 3)
  sccs_dag <- TarjanIterative(adj_dag)
  expect_equal(length(sccs_dag), 3)
  expect_true(all(lengths(sccs_dag) == 1))

  # A graph containing an isolated (zero out-degree) vertex must not error
  adj_isolated <- sp(1, 2, 5)
  expect_silent(GSCC(adj_isolated))
})

test_that("Wedges()/Trust_Triangles()/Cycle_Triangles() match brute-force counts", {
  # Small, hand-checkable directed graph:
  # 1 -> 2 -> 3 -> 1 (a 3-cycle) plus 1 -> 4
  adj <- sp(c(1, 2, 3, 1), c(2, 3, 1, 4), 4)

  w <- Wedges(adj)
  # Wedges v -> w -> u, u != v: (1,2,3),(2,3,1),(3,1,2),(3,1,4)
  expect_equal(w$NoWedges, 4)
  expect_equal(nrow(w$Wedge_list), 4)

  tt <- Trust_Triangles(adj)
  # Trust triangles: v -> w and a common out-neighbor u of both.
  # v=1,w=2: neigh(1)={2,4}, neigh(2)={3} -> no common neighbor
  # v=2,w=3: neigh(2)={3}, neigh(3)={1} -> no common neighbor
  # v=3,w=1: neigh(3)={1}, neigh(1)={2,4} -> no common neighbor
  expect_equal(tt$NoTriangles, 0)

  ct <- Cycle_Triangles(adj)
  # Exactly one 3-cycle: 1 -> 2 -> 3 -> 1
  expect_equal(ct$NoTriangles, 1)
  expect_equal(nrow(ct$Triangle_list), 3)  # counted once per vertex on the cycle
})

test_that("Count_Types() handles an empty (zero-row) object list", {
  # Regression test: do.call(rbind, list()) returns NULL, and 1:nrow(NULL)
  # used to produce c(1, 0) instead of an empty sequence.
  empty <- matrix(nrow = 0, ncol = 3)
  V <- c(1L, 2L, 1L)
  expect_silent(res <- Count_Types(empty, V, maxTypes = 2))
  expect_equal(sum(res), 0)
})
