#' Compute the RaCInG kernel for one or more patients
#'
#' The kernel is the *unweighted* structural kernel: cell-type abundances are
#' not baked into it. `calculateDirect()`, `calculateWedges()`,
#' `computeTriangles()`, and `computeGSCC()` each multiply in the abundance of
#' every node they use, exactly once per node, when computing raw
#' (unnormalized) features from this kernel.
#'
#' @param liglist Cell-by-ligand compatibility matrix.
#' @param reclist Cell-by-receptor compatibility matrix.
#' @param Cmatrix Patient-by-cell-type abundance matrix.
#' @param LRmatrix Ligand-by-receptor-by-patient interaction tensor.
#' @param normalize Logical; if `TRUE`, also compute a second kernel
#'   (`kernel_norm`) under a NULL DISTRIBUTION where every ligand-receptor
#'   pair that is structurally possible (nonzero in `LRmatrix`) is given the
#'   same uniform strength (`1 / number of active pairs for that patient`),
#'   instead of its real measured expression-derived strength. `kernel_norm`
#'   therefore reflects only network topology and cell-type abundance --
#'   "how much communication would you expect between these two cell types if
#'   every possible ligand-receptor pair were equally active" -- with no
#'   information about which pairs are actually more or less expressed.
#'   Downstream functions (`calculateDirect()`, `calculateWedges()`, etc.) use
#'   this as a baseline to compute a composition-independent enrichment score
#'   (`kernel / kernel_norm`) instead of the raw, abundance-weighted score --
#'   the abundance weighting cancels out of that ratio algebraically, since
#'   both `kernel` and `kernel_norm` are built from the same per-patient
#'   `lig_weight`/`rec_weight` terms. Use `normalize = FALSE` (default) when
#'   you want absolute communication strength (abundance and real LR
#'   expression both matter); use `normalize = TRUE` when you want a score
#'   that isolates specificity/enrichment of a cell-type pair's communication
#'   relative to what topology alone would predict, independent of how common
#'   those cell types are or how strong LR expression is overall.
#'
#' @return Either a 3D kernel array `[sender x receiver x patient]`, or (when
#'   `normalize = TRUE`) a list with `kernel` and `kernel_norm`.
#' @export
compute_kernel <- function(liglist, reclist, Cmatrix, LRmatrix, normalize = FALSE) {
  # liglist:   [cell type × ligand] binary/weighted matrix mapping cell types to ligands
  # reclist:   [cell type × receptor] binary/weighted matrix mapping cell types to receptors
  # Cmatrix:   [patient × cell type] cell type abundances for each patient
  # LRmatrix:  [ligand × receptor × patient] ligand-receptor interaction strengths for each patient
  # normalize: if TRUE, also compute the kernel under a uniformized (null-distribution) LRmatrix --
  #            see @param normalize above for what this represents and why you'd want it.
  #
  # Return:
  #  - if normalize == FALSE: kernel array [sender x receiver x patient]
  #  - if normalize == TRUE: list(kernel = kernel, kernel_norm = normalized_kernel)
  #
  # NOTE: the returned kernel is the *unweighted* structural kernel (cell-type
  # abundances are NOT baked in). Each feature-extraction function multiplies
  # in the abundance of every node it uses, exactly once per node, mirroring
  # the original Python implementation.

  # Ensure matrix format
  Cmatrix <- if (is.null(dim(Cmatrix))) matrix(Cmatrix, nrow = 1) else Cmatrix
  n_patients <- nrow(Cmatrix)
  n_celltypes <- ncol(Cmatrix)

  # Prepare normalized LRmatrix copy if requested (do not overwrite original LRmatrix)
  if (normalize) {
    LRmatrix_norm <- LRmatrix
    normvec <- 1 / apply(LRmatrix_norm > 0, 3, sum)  # inverse of positive interactions per patient
    for (p in 1:n_patients) {
      copy <- LRmatrix_norm[, , p]
      copy[copy > 0] <- normvec[p]           # scale positive interactions
      LRmatrix_norm[, , p] <- copy
    }
  }

  # Output: [sender x receiver x patient]
  kernel <- array(0, dim = c(n_celltypes, n_celltypes, n_patients))
  if (normalize) kernel_norm <- array(0, dim = c(n_celltypes, n_celltypes, n_patients))

  for (p in 1:n_patients) {

    Cvec <- Cmatrix[p, ]  # cell abundances

    # ---- ligand weights ----
    # weightlig[i] = total abundance-weighted ligand i expression across cell types
    lig_norm <- colSums(sweep(liglist, 1, Cvec, "*"))
    lig_norm[lig_norm == 0] <- 1  # avoid division by zero
    # lig_weight[A, i] = P(ligand i is sent by cell type A)
    lig_weight <- sweep(liglist, 2, lig_norm, "/")

    # ---- receptor weights ----
    rec_norm <- colSums(sweep(reclist, 1, Cvec, "*"))
    rec_norm[rec_norm == 0] <- 1
    # rec_weight[B, j] = P(receptor j is received by cell type B)
    rec_weight <- sweep(reclist, 2, rec_norm, "/")

    # ---- aggregate over ligand-receptor pairs ----
    # kernel[A, B] = sum over (i,j) of P(A sends i) * LRmatrix[i,j] * P(B receives j)
    kernel[, , p] <- lig_weight %*% LRmatrix[, , p] %*% t(rec_weight)

    if (normalize) {
      kernel_norm[, , p] <- lig_weight %*% LRmatrix_norm[, , p] %*% t(rec_weight)
    }
  }

  if (normalize) {
    return(list(kernel = kernel, kernel_norm = kernel_norm))
  } else {
    return(kernel)
  }
}

#' Calculate direct communication features from a kernel
#'
#' @param kernel Kernel array returned by [compute_kernel()].
#' @param unifKernel Optional normalized baseline kernel.
#' @param cells Character vector of cell-type names.
#' @param Dcell Patient-by-cell-type abundance matrix, used to weight raw
#'   (unnormalized) scores. Ignored when `unifKernel` is supplied, since the
#'   abundance weights cancel out of the ratio.
#' @param bundle Logical; if `TRUE`, combine reciprocal directions.
#'
#' @return A patient-by-feature data frame of direct communication scores.
#' @export
calculateDirect <- function(kernel, unifKernel = NULL, cells, Dcell = NULL, bundle = TRUE) {
  n_celltypes <- length(cells)      # number of cell types
  n_patients <- dim(kernel)[3]      # number of patients

  df_list <- list()                 # list to store each column before converting to data frame

  # Loop over all cell type pairs (sender i, receiver j).
  # If bundle = TRUE, only the upper triangle is used to avoid double counting.
  for (i in 1:n_celltypes) {
    for (j in if (bundle) i:n_celltypes else 1:n_celltypes) {

      colname <- paste0("Dir_", cells[i], "_", cells[j])

      if (!is.null(unifKernel)) {
        # normalized direct score: abundance weights cancel out of the ratio
        value <- if (bundle) {
          (kernel[i, j, ] + kernel[j, i, ]) / (unifKernel[i, j, ] + unifKernel[j, i, ])
        } else {
          kernel[i, j, ] / unifKernel[i, j, ]
        }
        value[is.nan(value)] <- 0 # no possible ligand-receptor pathway between i and j (0/0): no communication, not undefined
      } else {
        # raw direct score: weight explicitly by the abundance of each node
        value <- if (bundle) {
          Dcell[, i] * Dcell[, j] * (kernel[i, j, ] + kernel[j, i, ])
        } else {
          Dcell[, i] * Dcell[, j] * kernel[i, j, ]
        }
      }

      df_list[[colname]] <- value   # add the vector of patient scores as a new column
    }
  }

  df <- as.data.frame(df_list)
  rownames(df) <- paste0("Patient_", 1:n_patients)

  return(df)
}

#' Calculate wedge features from a kernel
#'
#' @param kernel Kernel array returned by [compute_kernel()].
#' @param unifKernel Optional normalized baseline kernel.
#' @param cells Character vector of cell-type names.
#' @param Dcell Patient-by-cell-type abundance matrix, used to weight raw
#'   (unnormalized) scores. Ignored when `unifKernel` is supplied.
#' @param bundle Logical; if `TRUE`, combine directionally equivalent wedges.
#'
#' @return A patient-by-feature data frame of wedge scores.
#' @export
calculateWedges <- function(kernel, unifKernel = NULL, cells, Dcell = NULL, bundle = TRUE) {
  n_celltypes <- length(cells)
  n_patients <- dim(kernel)[3]

  df_list <- list()   # list to collect columns before converting to data frame

  # Loop over all triplets of cell types (i, j, k): i = sender, j = intermediate, k = receiver
  for (i in 1:n_celltypes) {
    for (j in 1:n_celltypes) {
      for (k in if (bundle) j:n_celltypes else 1:n_celltypes) {

        colname <- paste0("W_", cells[i], "_", cells[j], "_", cells[k])

        if (!is.null(unifKernel)) {
          # normalized wedge score: abundance weights cancel out of the ratio
          value <- if (bundle) {
            (kernel[i,j,] * kernel[j,k,] + kernel[k,j,] * kernel[j,i,]) /
              (unifKernel[i,j,] * unifKernel[j,k,] + unifKernel[k,j,] * unifKernel[j,i,])
          } else {
            (kernel[i,j,] * kernel[j,k,]) /
              (unifKernel[i,j,] * unifKernel[j,k,])
          }
          value[is.nan(value)] <- 0 # no possible ligand-receptor pathway through this triplet (0/0): no communication, not undefined
        } else {
          # raw wedge: weight explicitly by the abundance of each of the 3 nodes
          value <- if (bundle) {
            Dcell[,i] * Dcell[,j] * Dcell[,k] * (kernel[i,j,] * kernel[j,k,] + kernel[k,j,] * kernel[j,i,])
          } else {
            Dcell[,i] * Dcell[,j] * Dcell[,k] * (kernel[i,j,] * kernel[j,k,])
          }
        }

        df_list[[colname]] <- value   # store column for this triplet
      }
    }
  }

  df <- as.data.frame(df_list)
  rownames(df) <- paste0("Patient_", 1:n_patients)

  return(df)
}

#' Legacy GSCC helper
#'
#' This helper depended on project-specific files from the original development
#' workflow. The packaged interface now recommends using [computeGSCC()] directly.
#'
#' @param cancer Legacy dataset identifier.
#' @param lab Scaling factor.
#' @param norm Logical; kept for backward compatibility.
#' @param test Logical; kept for backward compatibility.
#'
#' @return This function stops with a message directing users to [computeGSCC()].
#' @export
getGSCCAnalytically <- function(cancer, lab = 1, norm = TRUE, test = FALSE) {
  stop(
    "`getGSCCAnalytically()` is a legacy helper. Use `computeGSCC()` with an explicit kernel and abundance matrix instead.",
    call. = FALSE
  )
}

#' Compute GSCC features from kernel matrices
#'
#' @param kernel Kernel array from [compute_kernel()].
#' @param Dcell Patient-by-cell-type abundance matrix.
#' @param cell_names Character vector of cell-type labels.
#' @param patient_names Character vector of patient names.
#' @param unifKernel Optional normalized baseline kernel.
#' @param norm Logical; if `TRUE`, divide by the baseline GSCC values.
#' @param lab Scaling factor for interaction strengths.
#'
#' @return A data frame with GSCC feature values per patient.
#' @export
computeGSCC <- function(kernel, Dcell, cell_names, patient_names,
                        unifKernel = NULL, norm = FALSE, lab = 1) {
  # kernel: [sender x receiver x patient] matrix from compute_kernel
  # Dcell: [patient x cell type] cell abundances (can be used to scale if needed)
  # cell_names: vector of cell type names
  # patient_names: vector of patient names
  # unifKernel: optional uniform kernel (same dimensions as kernel) for normalization
  # norm: boolean, whether to normalize by uniform kernel
  # lab: scaling factor for interaction strengths

  n_patients <- nrow(Dcell)
  n_cells <- length(cell_names)

  # Matrices to store GSCC per patient and per cell type
  GSCCsizes <- matrix(0, nrow = n_patients, ncol = n_cells + 1)  # last column = total GSCC
  GSCCsizesN <- if (!is.null(unifKernel)) matrix(0, nrow = n_patients, ncol = n_cells + 1) else NULL

  # Loop over patients
  for (p in 1:n_patients) {
    q <- Dcell[p, ]  # cell abundances vector for this patient

    # Construct mu matrices for Poisson branching process
    muP <- matrix(0, nrow = n_cells, ncol = n_cells)
    muM <- matrix(0, nrow = n_cells, ncol = n_cells)
    for (i in 1:n_cells) {
      for (j in 1:n_cells) {
        muP[i, j] <- lab * kernel[j, i, p] * q[j]  # incoming contribution
        muM[i, j] <- lab * kernel[i, j, p] * q[j]  # outgoing contribution
      }
    }

    # Solve the branching process equations using nleqslv
    solP <- nleqslv::nleqslv(rep(1, n_cells), poiBPFunc, M = muP, sens = n_cells)
    solM <- nleqslv::nleqslv(rep(1, n_cells), poiBPFunc, M = muM, sens = n_cells)

    x <- solP$x
    y <- solM$x

    # GSCC per cell type and total
    GSCCsizes[p, 1:n_cells] <- x * y * q
    GSCCsizes[p, n_cells + 1] <- sum(x * y * q)

    # If uniform kernel is provided, calculate same for normalization
    if (!is.null(unifKernel)) {
      muP_N <- matrix(0, nrow = n_cells, ncol = n_cells)
      muM_N <- matrix(0, nrow = n_cells, ncol = n_cells)
      for (i in 1:n_cells) {
        for (j in 1:n_cells) {
          muP_N[i, j] <- lab * unifKernel[j, i, p] * q[j]
          muM_N[i, j] <- lab * unifKernel[i, j, p] * q[j]
        }
      }
      solP_N <- nleqslv::nleqslv(rep(1, n_cells), poiBPFunc, M = muP_N, sens = n_cells)
      solM_N <- nleqslv::nleqslv(rep(1, n_cells), poiBPFunc, M = muM_N, sens = n_cells)
      xN <- solP_N$x
      yN <- solM_N$x
      GSCCsizesN[p, 1:n_cells] <- xN * yN * q
      GSCCsizesN[p, n_cells + 1] <- sum(xN * yN * q)
    }
  }

  # Normalize by uniform kernel if requested
  if (norm && !is.null(unifKernel)) {
    normvals <- GSCCsizes / GSCCsizesN
    normvals[is.nan(normvals)] <- 1  # handle division by zero
    df <- as.data.frame(normvals)
  } else {
    df <- as.data.frame(GSCCsizes)
  }

  # Set column and row names
  colnames(df) <- c(paste0("GSCC_", cell_names), "GSCC_total")
  rownames(df) <- patient_names

  return(df)
}

#' Compute triangle features from kernel matrices
#'
#' @param kernel Kernel array from [compute_kernel()].
#' @param cell_names Character vector of cell-type names.
#' @param patient_names Character vector of patient names.
#' @param Dcell Patient-by-cell-type abundance matrix, used to weight raw
#'   (unnormalized) scores. Ignored when `unifKernel`/`norm` is used.
#' @param unifKernel Optional normalized baseline kernel.
#' @param norm Logical; if `TRUE`, divide by the baseline triangle scores.
#' @param bundle Logical; if `TRUE`, aggregate all directions into a single
#'   "Tr" triangle score per triple. If `FALSE`, two scores are returned per
#'   triple: "TT" (trust/transitive triangle, `i->j->k` and `i->k`) and,
#'   for `i<=j<=k` only, "CT" (cycle triangle, `i->j->k->i`).
#'
#' @return A patient-by-feature data frame of triangle scores.
#' @export
computeTriangles <- function(kernel, cell_names, patient_names, Dcell = NULL,
                             unifKernel = NULL, norm = FALSE, bundle = TRUE) {
  n_patients <- dim(kernel)[3]   # number of patients
  n_cells <- length(cell_names)   # number of cell types

  df_list <- list()  # will store triangle values per patient

  weight_or_ratio <- function(num_fun, i, j, k) {
    # num_fun(K) computes the triangle numerator from a kernel array K
    if (norm && !is.null(unifKernel)) {
      value <- num_fun(kernel) / num_fun(unifKernel)
      value[is.nan(value)] <- 1
    } else {
      value <- Dcell[, i] * Dcell[, j] * Dcell[, k] * num_fun(kernel)
    }
    value
  }

  for (i in 1:n_cells) {
    for (j in if (bundle) i:n_cells else 1:n_cells) {
      for (k in if (bundle) j:n_cells else 1:n_cells) {

        if (bundle) {
          # Sum all 8 edge-direction combinations of the 3 undirected edges
          # (2 cyclic + 6 transitive) to remove directionality entirely.
          triangle_sum <- function(K) {
            K[i,j,]*K[j,k,]*K[k,i,] + K[i,j,]*K[j,k,]*K[i,k,] +
            K[i,j,]*K[k,j,]*K[i,k,] + K[i,j,]*K[k,j,]*K[k,i,] +
            K[j,i,]*K[k,j,]*K[k,i,] + K[j,i,]*K[j,k,]*K[k,i,] +
            K[j,i,]*K[j,k,]*K[i,k,] + K[j,i,]*K[k,j,]*K[i,k,]
          }
          df_list[[paste0("Tr_", cell_names[i], "_", cell_names[j], "_", cell_names[k])]] <-
            weight_or_ratio(triangle_sum, i, j, k)

        } else {
          # Trust triangle (transitive): i->j, j->k, i->k
          trust_triangle <- function(K) K[i,j,] * K[j,k,] * K[i,k,]
          df_list[[paste0("TT_", cell_names[i], "_", cell_names[j], "_", cell_names[k])]] <-
            weight_or_ratio(trust_triangle, i, j, k)

          # Cycle triangle: i->j->k->i (only counted once per unordered triple)
          if (i <= j && j <= k) {
            cycle_triangle <- function(K) K[i,j,] * K[j,k,] * K[k,i,]
            df_list[[paste0("CT_", cell_names[i], "_", cell_names[j], "_", cell_names[k])]] <-
              weight_or_ratio(cycle_triangle, i, j, k)
          }
        }
      }
    }
  }

  # Convert the list into a dataframe: rows = patients, columns = triangles
  df <- as.data.frame(df_list)
  rownames(df) <- patient_names

  return(df)
}

#' Derive communication features from a kernel
#'
#' @param kernel Kernel array from [compute_kernel()].
#' @param unifKernel Optional normalized baseline kernel.
#' @param celltypes Character vector of cell-type labels.
#' @param communication_type Feature family to compute (`"D"`, `"W"`, `"TT"`, or `"GSCC"`).
#' @param bundle Logical; if `TRUE`, merge directionally equivalent features where appropriate.
#' @param patient_names Optional patient labels.
#' @param Dcell Patient-by-cell-type abundance matrix. Always required for
#'   `"GSCC"`; required for `"D"`, `"W"`, `"TT"` only when `unifKernel` is not
#'   supplied (raw, unnormalized features).
#' @param norm Logical; if `TRUE`, compute normalized features when a baseline is supplied.
#' @param patient_idx Optional patient index subset.
#'
#' @return A data frame of feature values for the selected patients.
#' @export
compute_kernel_features <- function(kernel, unifKernel = NULL, celltypes, communication_type = "D", bundle = TRUE,
                                    patient_names = NULL, Dcell = NULL, norm = FALSE, patient_idx = NULL) {

  # Subset patients if requested
  if (!is.null(patient_idx)) {
    kernel <- kernel[,,patient_idx, drop = FALSE]
    if (!is.null(unifKernel)) unifKernel <- unifKernel[,,patient_idx, drop = FALSE]
    if (!is.null(patient_names)) patient_names <- patient_names[patient_idx]
    if (!is.null(Dcell)) Dcell <- Dcell[patient_idx, , drop = FALSE]
  }

  comm <- toupper(communication_type)

  if (comm %in% c("D", "W", "TT") && is.null(Dcell) && is.null(unifKernel)) {
    stop("Dcell (cell abundance matrix) is required to compute raw (unnormalized) ", comm, " features.")
  }

  if (comm == "D") {
    df <- calculateDirect(kernel = kernel, unifKernel = unifKernel, cells = celltypes, Dcell = Dcell, bundle = bundle)

  } else if (comm == "W") {
    df <- calculateWedges(kernel = kernel, unifKernel = unifKernel, cells = celltypes, Dcell = Dcell, bundle = bundle)

  } else if (comm == "TT") {
    # trust/cycle triangles (directed): use bundle = FALSE in computeTriangles
    df <- computeTriangles(kernel = kernel, cell_names = celltypes, Dcell = Dcell,
                           patient_names = if (!is.null(patient_names)) patient_names else paste0("Patient_", seq_len(dim(kernel)[3])),
                           unifKernel = unifKernel, norm = norm, bundle = FALSE)

  } else if (comm == "GSCC") {
    if (is.null(Dcell)) stop("Dcell (cell abundance matrix) is required for GSCC computation.")
    # computeGSCC expects patient_names and Dcell rows == patients
    df <- computeGSCC(kernel = kernel, Dcell = Dcell,
                      cell_names = celltypes,
                      patient_names = if (!is.null(patient_names)) patient_names else paste0("Patient_", seq_len(nrow(Dcell))),
                      unifKernel = unifKernel, norm = norm)

  } else {
    stop("Unsupported communication_type. Use D, W, TT, or GSCC.")
  }

  return(df)

}

#' Run the full kernel-based RaCInG workflow
#'
#' @param counts Gene-by-sample count matrix. Required when `input_data` is not
#'   supplied; ignored otherwise.
#' @param output_folder Directory used to write and read intermediate input files.
#' @param deconv Optional patient-by-cell-type abundance matrix. If omitted, it
#'   is computed via `multideconv::compute.deconvolution()` followed by
#'   `multideconv::compute.deconvolution.analysis()` (which identifies and
#'   collapses correlated cell-type subgroups) and
#'   `multideconv::standardize_celltype_colnames()`. See [prepare_input_files()].
#' @param cc_network Optional ligand-receptor prior network. If omitted, it is
#'   retrieved via `liana::get_curated_omni()`. See [prepare_input_files()].
#' @param fun_LR Function used to combine ligand and receptor expression values.
#' @param cell_expr_profile Optional gene-by-cell-type expression profile
#'   matrix. If omitted, it is estimated from `counts` and `deconv` via
#'   per-gene non-negative least squares. See [prepare_input_files()].
#' @param source,target Column names to use as ligand and receptor identifiers in `cc_network`.
#' @param signed Logical; if `TRUE`, also try to load a sign matrix.
#' @param deconv_method Deconvolution method(s) used when `deconv` is not supplied.
#' @param cbsx.name,cbsx.token Optional credentials for the deconvolution workflow.
#' @param file_name File stem used for intermediate input files.
#' @param nPatients Number of patients to process, or `"all"`.
#' @param communication_type Feature family to compute.
#' @param norm Logical; if `TRUE`, compute a normalized baseline kernel.
#' @param pt_idx Optional single patient index to process.
#' @param remove_direction Logical; if `TRUE`, merge directionally equivalent features.
#' @param input_data Optional named list of pre-computed input matrices as returned
#'   by [prepare_input_files()].
#'   Must contain `Lmatrix`, `Rmatrix`, `Cmatrix`, `LRmatrix`, `celltypes`,
#'   `ligands`, and `receptors`.
#'   When supplied, the `counts` argument and all preprocessing parameters
#'   (`deconv`, `cc_network`, etc.) are ignored.
#'
#' @return A list with the kernel arrays and the derived feature matrix.
#' @export
compute_racing_kernel = function(counts = NULL, output_folder = "~/Documents/racing/vignettes/", deconv = NULL, cc_network = NULL, fun_LR = min, 
                                 cell_expr_profile = NULL, source = "source_genesymbol", target = "target_genesymbol", signed = FALSE,
                                 deconv_method = "Quantiseq", cbsx.name = NULL, cbsx.token = NULL, file_name = NULL, nPatients = "all", 
                                 communication_type = "W", norm = TRUE, pt_idx = NULL, remove_direction = TRUE,
                                 input_data = NULL) {

  if (is.null(file_name)) {
    file_name <- "RaCInG_input"
  }

  if (!is.null(input_data)) {
    cat("Using pre-computed input matrices; skipping input generation.\n")
    res <- input_data
  } else {
    if (is.null(counts)) {
      stop("Either 'counts' or 'input_data' must be provided.", call. = FALSE)
    }
    res <- prepare_input_files(counts, output_folder = output_folder, deconv = deconv, cc_network = cc_network, fun_LR = fun_LR, 
                               cell_expr_profile = cell_expr_profile, source = source, target = target,
                               deconv_method = deconv_method, cbsx.name = cbsx.name, cbsx.token = cbsx.token, file_name = file_name,
                               signed = signed)
  }

  Lmatrix   <- res$Lmatrix
  Rmatrix   <- res$Rmatrix
  Cmatrix    <- res$Cmatrix
  LRmatrix   <- res$LRmatrix
  cellTypes <- res$celltypes
  ligs      <- res$ligands
  recs      <- res$receptors

  if(nPatients == "all"){
    nPatients = nrow(Cmatrix)
  }                        

  if(!is.null(pt_idx)){
    cat("Because patient index is provided, nPatients argument will be ignored.\n")
  }else{
    cat("Computing kernel for ", ifelse(nPatients == "all", "all", nPatients), " patients\n")
    Cmatrix  <- Cmatrix[1:nPatients, , drop = FALSE]
    LRmatrix <- LRmatrix[,,1:nPatients, drop = FALSE]
  }

  # -----------------------------
  # Calculate kernel
  # -----------------------------
  cat("Calculating kernel...\n")
  
  res <- compute_kernel(Lmatrix, Rmatrix, Cmatrix, LRmatrix, normalize = norm)
  kernel_raw <- if (isTRUE(norm)) res$kernel else res
  kernel_norm <- if (isTRUE(norm)) res$kernel_norm else NULL

  # -----------------------------
  # Calculate features
  # -----------------------------
  cat("Calculating features...\n")
  patient_names <- if (!is.null(counts)) colnames(counts) else rownames(Cmatrix)
  features <- compute_kernel_features(
    kernel = kernel_raw,
    unifKernel = kernel_norm,
    celltypes = cellTypes,
    communication_type = communication_type,
    bundle = remove_direction,
    patient_names = patient_names,
    Dcell = Cmatrix,
    norm = norm,
    patient_idx = pt_idx
  )

  return(list(kernel = kernel_raw, kernel_norm = kernel_norm, features = features))

}