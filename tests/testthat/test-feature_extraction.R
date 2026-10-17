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
