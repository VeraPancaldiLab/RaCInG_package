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
#' @param patient_names Optional character vector of patient labels, matching
#'   the patient order in `kernel`. Defaults to `Patient_1`, `Patient_2`, ...
#'   when omitted, matching [compute_racing_montecarlo()]'s output when it
#'   is also given real names.
#'
#' @return A patient-by-feature data frame of direct communication scores.
#'   Cell-type pairs with no possible ligand-receptor pathway in any patient
#'   (zero for every patient) are dropped rather than returned as all-zero
#'   columns.
#' @export
calculateDirect <- function(kernel, unifKernel = NULL, cells, Dcell = NULL, bundle = TRUE, patient_names = NULL) {
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
  rownames(df) <- if (!is.null(patient_names)) patient_names else paste0("Patient_", 1:n_patients)

  # Drop cell-type pairs with no possible ligand-receptor pathway in any
  # patient (zero for every patient -- both branches above fill the
  # structurally-impossible NaN case with 0, so a plain all-zero check is a
  # correct and sufficient test here regardless of raw vs normalized).
  df <- df[, colSums(df != 0) > 0, drop = FALSE]

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
#' @param patient_names Optional character vector of patient labels, matching
#'   the patient order in `kernel`. Defaults to `Patient_1`, `Patient_2`, ...
#'   when omitted.
#'
#' @return A patient-by-feature data frame of wedge scores. Triplets with no
#'   possible ligand-receptor pathway in any patient (zero for every
#'   patient) are dropped rather than returned as all-zero columns.
#' @export
calculateWedges <- function(kernel, unifKernel = NULL, cells, Dcell = NULL, bundle = TRUE, patient_names = NULL) {
  n_celltypes <- length(cells)
  n_patients <- dim(kernel)[3]

  # Build the full (i, j, k) index grid vectorized, up front, instead of a
  # triple nested for-loop with a per-triplet df_list[[colname]] <- value
  # assignment -- that pattern grows an R named list one entry at a time,
  # which for realistic cell-type counts (e.g. 106) means n^2*(n+1)/2 (~600K,
  # bundled) or n^3 (~1.2M, unbundled) individual list insertions plus a
  # final as.data.frame() on all of them: this was the actual crash/hang,
  # not the underlying arithmetic. A single vectorized computation over the
  # whole grid does the same math with none of that per-element overhead
  # (mirrors the fix already applied to compute_results_processing() on the
  # Monte Carlo side for the same O(n^3) list-growth pattern).
  if (bundle) {
    grid <- do.call(rbind, lapply(seq_len(n_celltypes), function(jj) {
      expand.grid(i = seq_len(n_celltypes), j = jj, k = jj:n_celltypes)
    }))
  } else {
    grid <- expand.grid(i = seq_len(n_celltypes), j = seq_len(n_celltypes), k = seq_len(n_celltypes))
  }
  i <- grid$i; j <- grid$j; k <- grid$k

  # Flatten kernel's [sender, receiver, patient] array to [(sender*receiver), patient]
  # so kernel[a, b, ] for every (a, b) pair in the grid can be pulled out in
  # one indexed lookup per patient column, instead of one R-level `kernel[i,j,]`
  # subset per triplet.
  kernel_mat <- matrix(kernel, nrow = n_celltypes * n_celltypes, ncol = n_patients)
  idx2 <- function(a, b) a + (b - 1L) * n_celltypes
  Kij <- kernel_mat[idx2(i, j), , drop = FALSE]
  Kjk <- kernel_mat[idx2(j, k), , drop = FALSE]
  if (bundle) {
    Kkj0 <- kernel_mat[idx2(k, j), , drop = FALSE]
    Kji0 <- kernel_mat[idx2(j, i), , drop = FALSE]
    struct_core <- Kij * Kjk + Kkj0 * Kji0
  } else {
    struct_core <- Kij * Kjk
  }

  # Drop triplets with no possible ligand-receptor pathway in any patient
  # (zero structural kernel product for every patient) before doing any
  # further (abundance-weighting or ratio) computation on them -- unifKernel
  # is a rescaling of the same nonzero pattern as kernel (normalization only
  # rescales already-positive entries, see compute_kernel()), so checking
  # kernel's pattern alone is correct and sufficient for both raw and
  # normalized output. This both shrinks the returned feature count and cuts
  # the working-set size for the (possibly very large, e.g. ~600K rows at
  # 106 cell types) computation below.
  possible <- rowSums(struct_core != 0) > 0
  if (!all(possible)) {
    i <- i[possible]; j <- j[possible]; k <- k[possible]
    Kij <- Kij[possible, , drop = FALSE]
    Kjk <- Kjk[possible, , drop = FALSE]
    if (bundle) { Kkj0 <- Kkj0[possible, , drop = FALSE]; Kji0 <- Kji0[possible, , drop = FALSE] }
  }
  rm(struct_core)

  if (!is.null(unifKernel)) {
    unif_mat <- matrix(unifKernel, nrow = n_celltypes * n_celltypes, ncol = n_patients)
    UKij <- unif_mat[idx2(i, j), , drop = FALSE]
    UKjk <- unif_mat[idx2(j, k), , drop = FALSE]

    if (bundle) {
      UKkj <- unif_mat[idx2(k, j), , drop = FALSE]
      UKji <- unif_mat[idx2(j, i), , drop = FALSE]
      num <- Kij * Kjk + Kkj0 * Kji0
      den <- UKij * UKjk + UKkj * UKji
    } else {
      num <- Kij * Kjk
      den <- UKij * UKjk
    }
    value <- num / den
    value[is.nan(value)] <- 0 # no possible ligand-receptor pathway through this triplet (0/0): no communication, not undefined
  } else {
    # Weight by abundance per patient via a small per-patient loop
    # (n_patients is typically <=~150) rather than materializing three
    # separate full [n_combo x n_patients] copies of Dcell (Di, Dj, Dk) --
    # for n_combo in the hundreds of thousands to over a million (large
    # cell-type counts), those three extra full-size matrices were the
    # single biggest avoidable memory cost in this function; Dmat[p, i] for
    # one patient p is just a cheap n_combo-length index/recycle, not a copy.
    Dmat <- as.matrix(Dcell)  # [n_patients x n_celltypes]
    core <- if (bundle) Kij * Kjk + Kkj0 * Kji0 else Kij * Kjk
    value <- core
    for (p in seq_len(n_patients)) {
      value[, p] <- Dmat[p, i] * Dmat[p, j] * Dmat[p, k] * core[, p]
    }
  }

  out <- t(value)  # [n_patients x n_combo]
  colnames(out) <- paste0("W_", cells[i], "_", cells[j], "_", cells[k])
  df <- as.data.frame(out)
  rownames(df) <- if (!is.null(patient_names)) patient_names else paste0("Patient_", 1:n_patients)

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
#' @return A patient-by-feature data frame of triangle scores. Triples with
#'   no possible ligand-receptor pathway in any patient (zero for every
#'   patient) are dropped rather than returned as all-zero columns.
#' @export
computeTriangles <- function(kernel, cell_names, patient_names, Dcell = NULL,
                             unifKernel = NULL, norm = FALSE, bundle = TRUE) {
  n_patients <- dim(kernel)[3]   # number of patients
  n_cells <- length(cell_names)   # number of cell types

  # Build the (i, j, k) index grid(s) vectorized, up front, instead of a
  # triple nested for-loop with a per-triple df_list[[colname]] <- value
  # assignment -- for realistic cell-type counts this pattern reaches
  # n(n+1)(n+2)/6 (bundled) or a full n^3 (unbundled -- what "TT" always
  # requests via compute_kernel_features()) individual list insertions, which
  # was the actual crash/hang for larger cohorts (e.g. n=106 gives ~1.2M
  # triples), not the underlying arithmetic. See the equivalent fix in
  # calculateWedges() for the same pattern.
  kernel_mat <- matrix(kernel, nrow = n_cells * n_cells, ncol = n_patients)
  unif_mat <- if (!is.null(unifKernel)) matrix(unifKernel, nrow = n_cells * n_cells, ncol = n_patients) else NULL
  idx2 <- function(a, b) a + (b - 1L) * n_cells
  get_val <- function(mat, a, b) mat[idx2(a, b), , drop = FALSE]

  # Weight by abundance per patient via a small per-patient loop
  # (n_patients is typically <=~150) rather than materializing three
  # separate full [n_combo x n_patients] copies of Dcell (Di, Dj, Dk) --
  # for n_combo in the hundreds of thousands to over a million (large
  # cell-type counts), those three extra full-size matrices were the
  # single biggest avoidable memory cost in this function; Dmat[p, idx] for
  # one patient p is just a cheap n_combo-length index/recycle, not a copy.
  weight_or_ratio_vec <- function(num_mat, den_mat, i_idx, j_idx, k_idx) {
    if (norm && !is.null(unifKernel)) {
      value <- num_mat / den_mat
      value[is.nan(value)] <- 1
      value
    } else {
      Dmat <- as.matrix(Dcell)
      value <- num_mat
      for (p in seq_len(n_patients)) {
        value[, p] <- Dmat[p, i_idx] * Dmat[p, j_idx] * Dmat[p, k_idx] * num_mat[, p]
      }
      value
    }
  }

  if (bundle) {
    grid <- do.call(rbind, lapply(seq_len(n_cells), function(ii) {
      do.call(rbind, lapply(ii:n_cells, function(jj) {
        expand.grid(i = ii, j = jj, k = jj:n_cells)
      }))
    }))
    i <- grid$i; j <- grid$j; k <- grid$k

    # Sum all 8 edge-direction combinations of the 3 undirected edges
    # (2 cyclic + 6 transitive) to remove directionality entirely.
    Kij <- get_val(kernel_mat, i, j); Kjk <- get_val(kernel_mat, j, k); Kki <- get_val(kernel_mat, k, i)
    Kik <- get_val(kernel_mat, i, k); Kkj <- get_val(kernel_mat, k, j); Kji <- get_val(kernel_mat, j, i)
    num <- Kij*Kjk*Kki + Kij*Kjk*Kik + Kij*Kkj*Kik + Kij*Kkj*Kki +
           Kji*Kkj*Kki + Kji*Kjk*Kki + Kji*Kjk*Kik + Kji*Kkj*Kik

    # Drop triples with no possible ligand-receptor pathway in any patient
    # (zero structural kernel product for every patient), before the
    # (possibly ratio or abundance-weighting) computation below -- unifKernel
    # shares kernel's nonzero pattern (see calculateWedges() for why), so
    # checking the raw structural product alone is correct here too.
    possible <- rowSums(num != 0) > 0
    if (!all(possible)) {
      i <- i[possible]; j <- j[possible]; k <- k[possible]
      Kij <- Kij[possible,,drop=FALSE]; Kjk <- Kjk[possible,,drop=FALSE]; Kki <- Kki[possible,,drop=FALSE]
      Kik <- Kik[possible,,drop=FALSE]; Kkj <- Kkj[possible,,drop=FALSE]; Kji <- Kji[possible,,drop=FALSE]
      num <- num[possible,,drop=FALSE]
    }

    den <- NULL
    if (norm && !is.null(unifKernel)) {
      UKij <- get_val(unif_mat, i, j); UKjk <- get_val(unif_mat, j, k); UKki <- get_val(unif_mat, k, i)
      UKik <- get_val(unif_mat, i, k); UKkj <- get_val(unif_mat, k, j); UKji <- get_val(unif_mat, j, i)
      den <- UKij*UKjk*UKki + UKij*UKjk*UKik + UKij*UKkj*UKik + UKij*UKkj*UKki +
             UKji*UKkj*UKki + UKji*UKjk*UKki + UKji*UKjk*UKik + UKji*UKkj*UKik
    }

    value <- weight_or_ratio_vec(num, den, i, j, k)
    out <- t(value)
    colnames(out) <- paste0("Tr_", cell_names[i], "_", cell_names[j], "_", cell_names[k])

  } else {
    grid <- expand.grid(i = seq_len(n_cells), j = seq_len(n_cells), k = seq_len(n_cells))
    i_full <- grid$i; j_full <- grid$j; k_full <- grid$k

    # Kij/Kjk (and their unifKernel counterparts) are shared by both TT and
    # CT below, so computed once and reused; Kik/Kki are each used by only
    # one of the two and freed (rm + gc) right after, since at up to ~1.2M
    # combinations x n_patients each of these full-size matrices is itself a
    # substantial allocation.
    Kij <- get_val(kernel_mat, i_full, j_full)
    Kjk <- get_val(kernel_mat, j_full, k_full)
    UKij <- if (norm && !is.null(unifKernel)) get_val(unif_mat, i_full, j_full) else NULL
    UKjk <- if (norm && !is.null(unifKernel)) get_val(unif_mat, j_full, k_full) else NULL

    # --- Trust triangle (transitive): i->j, j->k, i->k -- every (i,j,k) triple ---
    Kik <- get_val(kernel_mat, i_full, k_full)
    num_tt <- Kij * Kjk * Kik
    rm(Kik)

    # Drop triples with no possible pathway in any patient (zero structural
    # kernel product for every patient), before the (possibly ratio or
    # abundance-weighting) computation -- unifKernel shares kernel's nonzero
    # pattern (see calculateWedges()), so checking the raw structural
    # product alone is correct here too.
    possible_tt <- rowSums(num_tt != 0) > 0
    i_tt <- i_full[possible_tt]; j_tt <- j_full[possible_tt]; k_tt <- k_full[possible_tt]
    num_tt <- num_tt[possible_tt, , drop = FALSE]
    den_tt <- if (norm && !is.null(unifKernel)) {
      (UKij * UKjk * get_val(unif_mat, i_full, k_full))[possible_tt, , drop = FALSE]
    } else NULL
    value_tt <- weight_or_ratio_vec(num_tt, den_tt, i_tt, j_tt, k_tt)
    rm(num_tt, den_tt); gc(FALSE)
    colnames_tt <- paste0("TT_", cell_names[i_tt], "_", cell_names[j_tt], "_", cell_names[k_tt])

    # --- Cycle triangle: i->j->k->i (only counted once per unordered triple, i<=j<=k) ---
    keep_ct <- i_full <= j_full & j_full <= k_full
    i_ct0 <- i_full[keep_ct]; j_ct0 <- j_full[keep_ct]; k_ct0 <- k_full[keep_ct]
    Kki <- get_val(kernel_mat, k_full, i_full)
    num_ct <- (Kij * Kjk * Kki)[keep_ct, , drop = FALSE]
    den_ct_full <- if (norm && !is.null(unifKernel)) {
      (UKij * UKjk * get_val(unif_mat, k_full, i_full))[keep_ct, , drop = FALSE]
    } else NULL
    rm(Kki, Kij, Kjk, UKij, UKjk)

    possible_ct <- rowSums(num_ct != 0) > 0
    i_ct <- i_ct0[possible_ct]; j_ct <- j_ct0[possible_ct]; k_ct <- k_ct0[possible_ct]
    num_ct <- num_ct[possible_ct, , drop = FALSE]
    den_ct <- if (!is.null(den_ct_full)) den_ct_full[possible_ct, , drop = FALSE] else NULL
    value_ct <- weight_or_ratio_vec(num_ct, den_ct, i_ct, j_ct, k_ct)
    rm(num_ct, den_ct); gc(FALSE)
    colnames_ct <- paste0("CT_", cell_names[i_ct], "_", cell_names[j_ct], "_", cell_names[k_ct])

    out <- t(rbind(value_tt, value_ct))
    colnames(out) <- c(colnames_tt, colnames_ct)
  }

  # Convert to a dataframe: rows = patients, columns = triangles
  df <- as.data.frame(out)
  rownames(df) <- patient_names

  return(df)
}

#' Derive communication features from a kernel
#'
#' @param kernel Kernel array from [compute_kernel()].
#' @param unifKernel Optional normalized baseline kernel.
#' @param celltypes Character vector of cell-type labels.
#' @param communication_type Feature family to compute (`"D"`, `"W"`, `"TT"`,
#'   or `"GSCC"`; `"TT"` returns both trust- and cycle-triangle columns
#'   together, see [computeTriangles()]), or a character vector of several of
#'   these. Since the kernel itself (`kernel`/`unifKernel`) is passed in
#'   already computed, requesting multiple types here costs nothing extra to
#'   compute per type beyond that one kernel -- no re-derivation happens for
#'   any of them. With more than one type, the return value is a named list
#'   (one data frame per type) instead of a single data frame.
#' @param bundle Logical; if `TRUE`, merge directionally equivalent features where appropriate.
#' @param patient_names Optional patient labels.
#' @param Dcell Patient-by-cell-type abundance matrix. Always required for
#'   `"GSCC"`; required for `"D"`, `"W"`, `"TT"` only when `unifKernel` is not
#'   supplied (raw, unnormalized features).
#' @param norm Logical; if `TRUE`, compute normalized features when a baseline is supplied.
#' @param patient_idx Optional patient index subset.
#'
#' @return A data frame of feature values for the selected patients, or (when
#'   `communication_type` has more than one entry) a named list of such data
#'   frames, one per requested type.
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

  .compute_one_type <- function(comm) {
    comm <- toupper(comm)

    if (comm %in% c("D", "W", "TT") && is.null(Dcell) && is.null(unifKernel)) {
      stop("Dcell (cell abundance matrix) is required to compute raw (unnormalized) ", comm, " features.")
    }

    if (comm == "D") {
      calculateDirect(kernel = kernel, unifKernel = unifKernel, cells = celltypes, Dcell = Dcell, bundle = bundle, patient_names = patient_names)

    } else if (comm == "W") {
      calculateWedges(kernel = kernel, unifKernel = unifKernel, cells = celltypes, Dcell = Dcell, bundle = bundle, patient_names = patient_names)

    } else if (comm == "TT") {
      # trust/cycle triangles (directed): use bundle = FALSE in computeTriangles
      computeTriangles(kernel = kernel, cell_names = celltypes, Dcell = Dcell,
                       patient_names = if (!is.null(patient_names)) patient_names else paste0("Patient_", seq_len(dim(kernel)[3])),
                       unifKernel = unifKernel, norm = norm, bundle = FALSE)

    } else if (comm == "GSCC") {
      if (is.null(Dcell)) stop("Dcell (cell abundance matrix) is required for GSCC computation.")
      # computeGSCC expects patient_names and Dcell rows == patients
      computeGSCC(kernel = kernel, Dcell = Dcell,
                 cell_names = celltypes,
                 patient_names = if (!is.null(patient_names)) patient_names else paste0("Patient_", seq_len(nrow(Dcell))),
                 unifKernel = unifKernel, norm = norm)

    } else {
      stop("Unsupported communication_type. Use D, W, TT, or GSCC.")
    }
  }

  if (length(communication_type) > 1) {
    out <- lapply(communication_type, .compute_one_type)
    names(out) <- communication_type
    return(out)
  }

  return(.compute_one_type(communication_type))
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
#' @param communication_type Feature family to compute (`"D"`, `"W"`, `"TT"`,
#'   or `"GSCC"`), or a character vector of several of these (e.g.
#'   `c("D", "W", "TT", "GSCC")`). The kernel (`compute_kernel()`) is always
#'   computed exactly once regardless of how many types are requested --
#'   feature extraction from an already-computed kernel is cheap, so there is
#'   no need to call `compute_racing_kernel()` again just to get a different
#'   `communication_type`; request them all in one call instead. With more
#'   than one type, `features` in the return value is a named list (one data
#'   frame per type) instead of a single data frame.
#' @param norm Logical; if `TRUE`, also compute a normalized baseline kernel and
#'   express features as an enrichment ratio over it (isolates specificity from
#'   abundance/topology). Default `FALSE` returns the raw, abundance-weighted
#'   communication magnitude.
#' @param pt_idx Optional single patient index to process.
#' @param remove_direction Logical; if `TRUE`, merge directionally equivalent features.
#' @param input_data Optional named list of pre-computed input matrices as returned
#'   by [prepare_input_files()].
#'   Must contain `Lmatrix`, `Rmatrix`, `Cmatrix`, `LRmatrix`, `celltypes`,
#'   `ligands`, and `receptors`.
#'   When supplied, the `counts` argument and all preprocessing parameters
#'   (`deconv`, `cc_network`, etc.) are ignored.
#'
#' @return A list with the kernel arrays and the derived feature matrix (or,
#'   for multiple `communication_type` entries, a named list of feature
#'   matrices) in `features`.
#' @export
compute_racing_kernel = function(counts = NULL, output_folder = "~/Documents/racing/vignettes/", deconv = NULL, cc_network = NULL, fun_LR = min, 
                                 cell_expr_profile = NULL, source = "source_genesymbol", target = "target_genesymbol", signed = FALSE,
                                 deconv_method = "Quantiseq", cbsx.name = NULL, cbsx.token = NULL, file_name = NULL, nPatients = "all", 
                                 communication_type = "W", norm = FALSE, pt_idx = NULL, remove_direction = TRUE,
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