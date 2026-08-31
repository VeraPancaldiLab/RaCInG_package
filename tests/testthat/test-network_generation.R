test_that("model1() generates the requested number of edges", {
  # Regression test: genRandomEdgeList() used to ignore the caller's edge
  # count entirely and always sample nrow(Dligrec) * ncol(Dligrec) edges
  # instead of round(avdeg * N), silently corrupting graph density for any
  # dataset where those two numbers differ (i.e. essentially always).
  set.seed(1)
  cellnames <- c("A", "B", "C")
  lig <- matrix(1, nrow = 3, ncol = 5, dimnames = list(cellnames, paste0("L", 1:5)))
  rec <- matrix(1, nrow = 3, ncol = 7, dimnames = list(cellnames, paste0("R", 1:7)))
  Dcell <- c(0.5, 0.3, 0.2)
  Dconn <- matrix(runif(5 * 7), nrow = 5)
  Dconn <- Dconn / sum(Dconn)

  N <- 100
  av <- 4
  g <- model1(N, av, lig, rec, Dcell, Dconn)

  # ligNo * recNo == 35 here, round(av * N) == 400: if the bug regresses,
  # nrow(E) collapses to 35 regardless of N/av.
  expect_equal(nrow(g$E), round(av * N))
})

test_that("genRandomEdgeList() honors an explicit edgeNo", {
  set.seed(2)
  Dligrec <- matrix(runif(9), nrow = 3)
  Dligrec <- Dligrec / sum(Dligrec)
  structurelig <- matrix(1, nrow = 4, ncol = 3)
  structurerec <- matrix(1, nrow = 4, ncol = 3)
  V <- sample(1:4, 50, replace = TRUE)

  edges <- genRandomEdgeList(Dligrec, V, structurelig, structurerec, edgeNo = 123)
  # All cells are compatible with everything here, so no edges are dropped.
  expect_equal(nrow(edges$cell_connection), 123)
  expect_equal(nrow(edges$ligrec_type), 123)
})
