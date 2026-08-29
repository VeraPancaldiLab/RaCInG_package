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
  # indices, which is correct for Direct communication but silently
  # double-counts shared/middle nodes for Wedges, Triangles, and (since the
  # branching-process equations are nonlinear) GSCC in both norm modes.
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

test_that("calculateWedges() does not double-weight the shared middle node", {
  d <- make_toy_kernel_inputs()
  kernel <- compute_kernel(d$liglist, d$reclist, d$Cmatrix, d$LRmatrix, normalize = FALSE)

  df <- calculateWedges(kernel, cells = c("A", "B"), Dcell = d$Cmatrix, bundle = FALSE)
  # W_B_A_B: i=B,j=A,k=B -> Dcell[B]*Dcell[A]*Dcell[B] * kernel[B,A]*kernel[A,B]
  Dcell <- d$Cmatrix[1, ]
  kAB <- kernel[1, 2, 1]
  kBA <- kernel[2, 1, 1]
  expected <- Dcell["B"] * Dcell["A"] * Dcell["B"] * kBA * kAB
  expect_equal(unname(df$W_B_A_B[1]), unname(expected), tolerance = 1e-8)
})

test_that("kernel features run for every communication type without NAs", {
  data(skcm_example, envir = environment())
  for (ct in c("D", "W", "TT", "GSCC")) {
    for (nm in c(TRUE, FALSE)) {
      res <- compute_racing_kernel(input_data = skcm_example, communication_type = ct, norm = nm)
      expect_false(any(is.na(res$features)), info = paste(ct, nm))
    }
  }
})
