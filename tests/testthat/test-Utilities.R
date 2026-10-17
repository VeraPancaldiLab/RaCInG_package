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
