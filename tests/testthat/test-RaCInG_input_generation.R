toy_inputs <- function() {
  genes <- c("L1", "R1", "L2", "R2", "G5")
  profile <- data.frame(A = c(100, 0, 0, 100, 1), C = c(0, 100, 100, 0, 1), row.names = genes)
  net <- data.frame(source_genesymbol = c("L1", "L2"), target_genesymbol = c("R1", "R2"))
  set.seed(1)
  samples <- paste0("Pt", 1:4, ".baseline")
  counts <- matrix(stats::rexp(5 * 4, 1 / 50), 5, 4, dimnames = list(genes, samples))
  deconv <- matrix(c(.5, .5, .3, .7, .8, .2, .1, .9), 4, 2, byrow = TRUE, dimnames = list(samples, c("A", "C")))
  list(profile = profile, net = net, counts = counts, deconv = deconv)
}

test_that("prepare_input_files() works without a trailing slash and keeps dotted sample names", {
  skip_if_not_installed("ADImpute")
  skip_if_not_installed("multideconv")
  t <- toy_inputs()
  out_dir <- tempfile() # no trailing slash, as tempdir()/tempfile() return
  inp <- suppressWarnings(prepare_input_files(t$counts, output_folder = out_dir, deconv = t$deconv, cc_network = t$net,
                                              cell_expr_profile = t$profile, file_name = "toy", already_normalized = TRUE))
  expect_equal(inp$celltypes, c("A", "C"))
  lr <- utils::read.csv(file.path(out_dir, "LRmatrix_toy.csv"), row.names = 1, check.names = FALSE)
  expect_equal(rownames(lr), colnames(t$counts))
})

test_that("compute_racing_kernel() accepts TPM counts via already_normalized", {
  skip_if_not_installed("ADImpute")
  skip_if_not_installed("multideconv")
  t <- toy_inputs()
  res <- suppressWarnings(compute_racing_kernel(t$counts, output_folder = tempfile(), deconv = t$deconv, cc_network = t$net,
                                                cell_expr_profile = t$profile, file_name = "toy", communication_type = "D",
                                                already_normalized = TRUE))
  expect_equal(nrow(res$features), ncol(t$counts))
})

test_that("workflow wrappers forward already_normalized and expr_threshold to prepare_input_files()", {
  captured <- NULL
  local_mocked_bindings(prepare_input_files = function(...) { captured <<- list(...); stop("captured") })
  expect_error(compute_racing_kernel(counts = matrix(1), already_normalized = TRUE, expr_threshold = 3), "captured")
  expect_true(captured$already_normalized); expect_equal(captured$expr_threshold, 3)
  expect_error(compute_racing_montecarlo(counts = matrix(1), already_normalized = TRUE, expr_threshold = 4), "captured")
  expect_true(captured$already_normalized); expect_equal(captured$expr_threshold, 4)
})

test_that("workflow wrappers default output_folder to Results/", {
  expect_equal(formals(compute_racing_kernel)$output_folder, "Results/")
  expect_equal(formals(compute_racing_montecarlo)$output_folder, "Results/")
  expect_equal(formals(runSim)$output_folder, "Results/")
})

test_that("prepare_input_files() stops when deconv rows don't match counts columns in order", {
  skip_if_not_installed("ADImpute")
  skip_if_not_installed("multideconv")
  t <- toy_inputs()
  run <- function(d) suppressWarnings(prepare_input_files(t$counts, output_folder = tempfile(), deconv = d, cc_network = t$net,
                                                          cell_expr_profile = t$profile, file_name = "toy", already_normalized = TRUE))
  expect_error(run(t$deconv[4:1, ]), "different order")
  expect_error(run(t$deconv[1:3, ]), "exactly the samples")
  expect_error(run(unname(t$deconv)), "exactly the samples")
})

test_that("generateInput() stops when Lmatrix and Rmatrix list different cell types", {
  dir <- tempfile(); dir.create(dir)
  utils::write.csv(matrix(1, 2, 1, dimnames = list(c("A", "B"), "L1")), file.path(dir, "Lmatrix_x.csv"))
  utils::write.csv(matrix(1, 2, 1, dimnames = list(c("A", "C"), "R1")), file.path(dir, "Rmatrix_x.csv"))
  expect_error(generateInput("x", dir), "same cell types")
})

test_that("prepare_input_files() keeps receive-only cell types", {
  skip_if_not_installed("ADImpute")
  skip_if_not_installed("multideconv")
  genes <- c("L1", "R1", "L2", "R2", "G5")
  # A sends L1 to B's R1; B expresses no ligand (receive-only); C talks to itself via L2-R2
  profile <- data.frame(A = c(100, 0, 0, 0, 1), B = c(0, 100, 0, 0, 1), C = c(0, 0, 100, 100, 1), row.names = genes)
  net <- data.frame(source_genesymbol = c("L1", "L2"), target_genesymbol = c("R1", "R2"))
  counts <- matrix(50 + seq_len(20), 5, 4, dimnames = list(genes, paste0("S", 1:4)))
  deconv <- matrix(c(.5, .3, .2, .2, .2, .6, .4, .4, .2, .1, .1, .8), 4, 3, byrow = TRUE,
                   dimnames = list(paste0("S", 1:4), c("A", "B", "C")))
  inp <- suppressWarnings(prepare_input_files(counts, output_folder = tempfile(), deconv = deconv, cc_network = net,
                                              cell_expr_profile = profile, file_name = "recv", already_normalized = TRUE))
  expect_equal(inp$celltypes, c("A", "B", "C"))
  expect_equal(unname(inp$Lmatrix[2, ]), c(0, 0)) # B sends nothing
  k <- compute_kernel(inp$Lmatrix, inp$Rmatrix, inp$Cmatrix, inp$LRmatrix)
  expect_true(all(k[1, 2, ] > 0)) # A -> B communication is kept
})

test_that("compute_racing_kernel() labels rows with the patients it used (nPatients with counts)", {
  skip_if_not_installed("ADImpute")
  skip_if_not_installed("multideconv")
  t <- toy_inputs()
  res <- suppressWarnings(compute_racing_kernel(t$counts, output_folder = tempfile(), deconv = t$deconv, cc_network = t$net,
                                                cell_expr_profile = t$profile, file_name = "toy", nPatients = 2,
                                                already_normalized = TRUE))
  expect_equal(rownames(res$features), colnames(t$counts)[1:2])
})

test_that("gene and cell-type names are kept unchanged", {
  skip_if_not_installed("ADImpute")
  skip_if_not_installed("multideconv")
  t <- toy_inputs()
  rownames(t$profile)[1] <- rownames(t$counts)[1] <- "HLA-A"
  t$net$source_genesymbol[1] <- "HLA-A"
  inp <- suppressWarnings(prepare_input_files(t$counts, output_folder = tempfile(), deconv = t$deconv, cc_network = t$net,
                                              cell_expr_profile = t$profile, file_name = "names", already_normalized = TRUE))
  expect_true("HLA-A" %in% inp$ligands)
  expect_equal(colnames(calculateDirect(array(1, c(2, 2, 1)), cells = c("B", "CD8+ T"), Dcell = matrix(0.5, 1, 2))),
               c("Dir_B_B", "Dir_B_CD8+ T", "Dir_CD8+ T_CD8+ T"))
})

test_that("user-supplied cell-type names outside multideconv's nomenclature give a warning", {
  skip_if_not_installed("multideconv")
  expect_warning(.warn_nonstandard_celltypes(c("A", "B.cells")), "not in multideconv's standardized nomenclature: A\\.")
  expect_no_warning(.warn_nonstandard_celltypes(c("B.cells", "B.cells_Subgroup.2", "Quantiseq_TIL10_CD8.cells", "Plasma_Subgroup.1")))
})
