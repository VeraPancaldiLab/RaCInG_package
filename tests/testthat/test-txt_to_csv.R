test_that("compute_results_processing() does not scramble cell types whose names contain underscores", {
  # Regression test: remove_direction used to group/relabel columns by
  # strsplit(colname, "_"), which silently corrupts merging whenever a
  # cell-type name itself contains "_" (e.g. multideconv's method-tagged
  # names like "CBSX_CBSX.HNSCC.scRNAseq_B.cells") -- there is no way to
  # tell, from the string alone, which underscores separate the two cell
  # types in a pair from underscores that are part of one cell type's own
  # name. Fixed by grouping on the underlying integer index grid instead.
  set.seed(1)
  cellnames <- c("A_1", "B_2", "C_3")
  lig <- matrix(1, nrow = 3, ncol = 4, dimnames = list(cellnames, paste0("L", 1:4)))
  rec <- matrix(1, nrow = 3, ncol = 4, dimnames = list(cellnames, paste0("R", 1:4)))
  Cmatrix <- matrix(c(0.5, 0.3, 0.2), nrow = 1, dimnames = list("Pat1", cellnames))
  LRmatrix <- array(runif(4 * 4 * 1), dim = c(4, 4, 1))
  LRmatrix <- LRmatrix / sum(LRmatrix)

  out_dir <- tempfile()
  dir.create(out_dir)
  runSim(lig, rec, Cmatrix, LRmatrix, cellnames, communication_type = "D",
         N = 20, itNo = 2, av = 3, output_folder = out_dir, file.name = "underscore_test")

  res <- compute_results_processing(
    celltypes = cellnames, patient_names = "Pat1", triangle_type = "D",
    remove_direction = TRUE, normalized = FALSE,
    sim_raw_file = file.path(out_dir, "underscore_test.out"),
    output_folder = out_dir, file.name = "underscore_test.csv"
  )

  expected <- paste0("Dir_", c(
    "A_1_A_1", "A_1_B_2", "A_1_C_3", "B_2_B_2", "B_2_C_3", "C_3_C_3"
  ))
  expect_setequal(colnames(res$mean), expected)
  expect_false(any(duplicated(colnames(res$mean))))
})
