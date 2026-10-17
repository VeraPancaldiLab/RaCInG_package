# Hand-checkable example: 2 cell types (A, B), 2 ligands, 2 receptors, 1 patient.
# A can only send L1 and receive R1; B can only send L2 and receive R2.
# The only nonzero interactions are L1-R2 and L2-R1, so A only talks to B and
# vice versa (A->A and B->B are structurally impossible).
make_toy_kernel_inputs <- function() {
  liglist <- matrix(c(1, 0, 0, 1), nrow = 2, byrow = TRUE, dimnames = list(c("A", "B"), c("L1", "L2")))
  reclist <- matrix(c(1, 0, 0, 1), nrow = 2, byrow = TRUE, dimnames = list(c("A", "B"), c("R1", "R2")))
  Cmatrix <- matrix(c(0.3, 0.7), nrow = 1, dimnames = list("P1", c("A", "B")))
  LRmatrix <- array(0, dim = c(2, 2, 1), dimnames = list(c("L1", "L2"), c("R1", "R2"), "P1"))
  LRmatrix[1, 2, 1] <- 1  # L1-R2
  LRmatrix[2, 1, 1] <- 1  # L2-R1
  list(liglist = liglist, reclist = reclist, Cmatrix = Cmatrix, LRmatrix = LRmatrix)
}

test_that("compute_kernel() returns the unweighted structural kernel", {
  # Regression test: compute_kernel() used to bake Cvec into both kernel
  # indices; abundances are applied once per node by calculateDirect() instead.
  d <- make_toy_kernel_inputs()
  k <- compute_kernel(d$liglist, d$reclist, d$Cmatrix, d$LRmatrix, normalize = FALSE)

  # kernel[A,B] = (1/weightlig_L1) * (1/weightrec_R2) = (1/0.3) * (1/0.7)
  expected <- (1 / 0.3) * (1 / 0.7)
  expect_equal(k[1, 2, 1], expected, tolerance = 1e-8)
  expect_equal(k[2, 1, 1], expected, tolerance = 1e-8)
  # No structural path for A->A or B->B
  expect_equal(k[1, 1, 1], 0)
  expect_equal(k[2, 2, 1], 0)
})

test_that("calculateDirect() weights raw scores by abundance exactly once per node", {
  d <- make_toy_kernel_inputs()
  kernel <- compute_kernel(d$liglist, d$reclist, d$Cmatrix, d$LRmatrix, normalize = FALSE)

  df <- calculateDirect(kernel, cells = c("A", "B"), Dcell = d$Cmatrix, bundle = FALSE)
  # Dir_A_B = Dcell[A] * Dcell[B] * kernel[A,B] = 0.3 * 0.7 * (1/0.21) == 1
  expect_equal(df$Dir_A_B[1], 1, tolerance = 1e-8)
  expect_equal(df$Dir_B_A[1], 1, tolerance = 1e-8)
})

test_that("kernel direct features run raw and normalized without NAs", {
  data(skcm_example, envir = environment())
  for (nm in c(TRUE, FALSE)) {
    res <- compute_racing_kernel(input_data = skcm_example, communication_type = "D", norm = nm)
    expect_false(any(is.na(res$features)), info = paste("norm =", nm))
  }
})

test_that("kernel method rejects removed or Monte-Carlo-only feature types", {
  data(skcm_example, envir = environment())
  for (ct in c("W", "TT", "GSCC")) {
    expect_error(compute_racing_kernel(input_data = skcm_example, communication_type = ct), "only computes direct", info = ct)
  }
  expect_error(compute_kernel_features(array(0, c(2, 2, 1)), celltypes = c("A", "B"), communication_type = c("D", "W"),
                                       Dcell = matrix(0.5, 1, 2)), "only computes direct")
})

test_that("Monte Carlo computes D and GSCC and rejects removed types", {
  data(skcm_example, envir = environment())
  inp <- skcm_example
  rownames(inp$Cmatrix) <- paste0("P", seq_len(nrow(inp$Cmatrix)))
  set.seed(1)
  res <- compute_racing_montecarlo(input_data = inp, output_folder = tempfile(), communication_type = c("D", "GSCC"),
                                   Ncells = 200, Ngraphs = 2, Ndegree = 3, nPatients = 2)
  expect_equal(names(res$output), c("D", "GSCC"))
  expect_equal(rownames(res$output$GSCC$mean), c("P1", "P2"))
  expect_true(all(startsWith(colnames(res$output$GSCC$mean), "GSCC_")))
  for (ct in c("W", "TT", "CT")) {
    expect_error(compute_racing_montecarlo(input_data = inp, output_folder = tempfile(), communication_type = ct), "wedges and triangles were removed", info = ct)
  }
})
