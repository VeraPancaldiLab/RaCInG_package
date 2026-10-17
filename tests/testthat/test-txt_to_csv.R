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
    celltypes = cellnames, patient_names = "Pat1", communication_type = "D",
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

sim_toy <- function() {
  set.seed(2)
  cellnames <- c("A", "B", "C")
  lig <- matrix(1, nrow = 3, ncol = 4, dimnames = list(cellnames, paste0("L", 1:4)))
  rec <- matrix(1, nrow = 3, ncol = 4, dimnames = list(cellnames, paste0("R", 1:4)))
  Cmatrix <- matrix(c(0.5, 0.3, 0.2, 0.2, 0.5, 0.3, 0.3, 0.3, 0.4), nrow = 3, byrow = TRUE,
                    dimnames = list(c("Pat1", "Pat2", "Pat3"), cellnames))
  LRmatrix <- array(runif(4 * 4 * 3), dim = c(4, 4, 3))
  list(cellnames = cellnames, lig = lig, rec = rec, Cmatrix = Cmatrix, LRmatrix = LRmatrix)
}

test_that("runSim() header counts only the simulated patients and results keep their names", {
  t <- sim_toy()
  out_dir <- tempfile()
  runSim(t$lig, t$rec, t$Cmatrix, t$LRmatrix, t$cellnames, communication_type = "D",
         N = 30, itNo = 2, av = 3, output_folder = out_dir, file.name = "subset", patient_idx = 2)
  sim_file <- file.path(out_dir, "subset.out")
  expect_equal(as.integer(strsplit(readLines(sim_file)[2], ",")[[1]][1]), 1L)
  expect_equal(Read_Sim_Output(sim_file)$patients, 2L)

  res <- compute_results_processing(t$cellnames, rownames(t$Cmatrix), "D", normalized = FALSE,
                                    sim_raw_file = sim_file, output_folder = out_dir, file.name = "subset.csv")
  expect_equal(rownames(res$mean), "Pat2")

  # Files written before the header fix state nrow(Cmatrix) instead; still read correctly.
  old_lines <- readLines(sim_file)
  old_lines[2] <- sub("^1,", "3,", old_lines[2])
  old_file <- file.path(out_dir, "old_header.out")
  writeLines(old_lines, old_file)
  old <- Read_Sim_Output(old_file)
  expect_equal(old$patients, 2L)
  expect_equal(dim(old$data_avg)[3], 1L)
})

test_that("runSim() writes into Results/ by default, not the filesystem root", {
  t <- sim_toy()
  wd <- tempfile(); dir.create(wd)
  old_wd <- setwd(wd); on.exit(setwd(old_wd), add = TRUE)
  runSim(t$lig, t$rec, t$Cmatrix, t$LRmatrix, t$cellnames, communication_type = "D", N = 30, itNo = 2, av = 3, pats = 1)
  expect_true(file.exists(file.path(wd, "Results", "RaCInG_input.out")))
})

test_that("countAllTypes() Direct counts match countDirect() on the same graphs", {
  t <- sim_toy()
  set.seed(3); a <- countDirect(t$Cmatrix[1, ], t$LRmatrix[, , 1], t$lig, t$rec, t$cellnames, N = 50, av = 4, itNo = 3)
  set.seed(3); b <- countAllTypes(t$Cmatrix[1, ], t$LRmatrix[, , 1], t$lig, t$rec, t$cellnames, N = 50, av = 4, itNo = 3, types = "D")$D
  expect_equal(b, a)
})

test_that("normalized Monte Carlo features report 0/0 as 0", {
  t <- sim_toy()
  t$Cmatrix[, "C"] <- 0 # C never appears, so every pair involving C is 0 in both simulations
  t$Cmatrix <- t$Cmatrix / rowSums(t$Cmatrix)
  out_dir <- tempfile()
  set.seed(4)
  for (nm in c(FALSE, TRUE)) {
    runSim(t$lig, t$rec, t$Cmatrix, t$LRmatrix, t$cellnames, communication_type = "D",
           N = 50, itNo = 2, av = 3, output_folder = out_dir, file.name = "zero", norm = nm)
  }
  res <- compute_results_processing(t$cellnames, rownames(t$Cmatrix), "D", normalized = TRUE,
                                    sim_raw_file = file.path(out_dir, "zero.out"),
                                    sim_norm_file = file.path(out_dir, "zero_norm.out"),
                                    output_folder = out_dir, file.name = "zero.csv")
  expect_true(all(res$mean[, c("Dir_A_C", "Dir_B_C", "Dir_C_C")] == 0))
  expect_true(all(res$sd[, c("Dir_A_C", "Dir_B_C", "Dir_C_C")] == 0))
})

test_that("compute_racing_montecarlo() works when Cmatrix has no rownames (e.g. skcm_example)", {
  t <- sim_toy()
  Cmat <- unname(t$Cmatrix)
  input <- list(Lmatrix = t$lig, Rmatrix = t$rec, Cmatrix = Cmat, LRmatrix = t$LRmatrix,
                celltypes = t$cellnames, ligands = colnames(t$lig), receptors = colnames(t$rec))
  out_dir <- tempfile()
  set.seed(4)
  res <- compute_racing_montecarlo(input_data = input, output_folder = out_dir, pt_idx = 2,
                                   Ncells = 30, Ngraphs = 2, Ndegree = 3)
  expect_equal(rownames(res$output$mean), "Patient_2")

  kres <- compute_racing_kernel(input_data = input, pt_idx = 2)
  expect_equal(rownames(kres$features), "Patient_2")
})

test_that("compute_results_processing() writes <stem>_sd.csv and only normalizes with a baseline file", {
  t <- sim_toy()
  out_dir <- tempfile()
  runSim(t$lig, t$rec, t$Cmatrix, t$LRmatrix, t$cellnames, communication_type = "D",
         N = 30, itNo = 2, av = 3, output_folder = out_dir, file.name = "stem")
  res <- compute_results_processing(t$cellnames, rownames(t$Cmatrix), "D",
                                    sim_raw_file = file.path(out_dir, "stem.out"),
                                    output_folder = out_dir, file.name = "stem.csv")
  expect_true(file.exists(file.path(out_dir, "stem_sd.csv")))
  expect_false(file.exists(file.path(out_dir, "stem.csv_sd.csv")))
  expect_equal(nrow(res$mean), 3L)
})

test_that("compute_racing_montecarlo()$input can be passed back as input_data", {
  t <- sim_toy()
  input <- list(Lmatrix = t$lig, Rmatrix = t$rec, Cmatrix = t$Cmatrix, LRmatrix = t$LRmatrix,
                celltypes = t$cellnames, ligands = colnames(t$lig), receptors = colnames(t$rec))
  out_dir <- tempfile()
  set.seed(5)
  res <- compute_racing_montecarlo(input_data = input, output_folder = out_dir, pt_idx = 1,
                                   Ncells = 30, Ngraphs = 2, Ndegree = 3)
  expect_named(res$input, c("LRmatrix", "Lmatrix", "Rmatrix", "Cmatrix", "celltypes", "ligands", "receptors"),
               ignore.order = TRUE)
  expect_equal(compute_racing_kernel(input_data = res$input)$features,
               compute_racing_kernel(input_data = input)$features)
})
