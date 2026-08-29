test_that("EdgetoAdj_No_loop() removes self-loops but keeps other edges", {
  E <- matrix(c(1, 1, 2, 1, 2, 3), ncol = 2, byrow = TRUE)  # 1->1 (loop), 2->1, 2->3
  adj <- EdgetoAdj_No_loop(E, 3)
  expect_equal(adj[1, 1], 0)
  expect_equal(adj[2, 1], 1)
  expect_equal(adj[2, 3], 1)
  expect_equal(sum(adj), 2)
})

test_that("EdgetoAdj() keeps self-loops", {
  E <- matrix(c(1, 1, 2, 1), ncol = 2, byrow = TRUE)
  adj <- EdgetoAdj(E, 2)
  expect_equal(adj[1, 1], 1)
  expect_equal(sum(adj), 2)
})

test_that("Count_Types() tabulates type combinations correctly", {
  oblist <- matrix(c(1, 2, 3, 2, 3, 1), nrow = 2, byrow = TRUE)  # rows are vertex triplets
  V <- c(1, 2, 1)  # vertex 1 -> type 1, vertex 2 -> type 2, vertex 3 -> type 1
  res <- Count_Types(oblist, V, maxTypes = 2)
  # row 1: vertices (1,2,3) -> types (1,2,1); row 2: vertices (2,3,1) -> types (2,1,1)
  expect_equal(res[1, 2, 1], 1)
  expect_equal(res[2, 1, 1], 1)
  expect_equal(sum(res), 2)
})
