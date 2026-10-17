#' Compute the RaCInG kernel for one or more patients
#'
#' The kernel is the *unweighted* structural kernel: cell-type abundances are
#' not baked into it. [calculateDirect()] multiplies in the abundance of both
#' cell types of a pair when computing raw (unnormalized) features from it.
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
#'   [calculateDirect()] uses this as a baseline to compute a composition-independent enrichment score
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
#' @examples
#' data(skcm_example)
#' k <- compute_kernel(skcm_example$Lmatrix, skcm_example$Rmatrix,
#'                     skcm_example$Cmatrix, skcm_example$LRmatrix)
#' dim(k)
#' feats <- compute_kernel_features(k, celltypes = skcm_example$celltypes,
#'                                  Dcell = skcm_example$Cmatrix,
#'                                  patient_names = rownames(skcm_example$Cmatrix))
#' head(feats[, 1:4])
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
  # abundances are NOT baked in). calculateDirect() multiplies in the abundance
  # of both cell types, mirroring the original Python implementation.

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
#'   when omitted.
#' @param zero_threshold Drop a cell-type pair once its fraction of
#'   zero-valued patients reaches this threshold. Default `1` only drops
#'   pairs that are zero for every patient (the original behavior); lower it
#'   (e.g. `0.9`) to also drop merely zero-inflated pairs.
#'
#' @return A patient-by-feature data frame of direct communication scores.
#'   Cell-type pairs with no possible ligand-receptor pathway in any patient
#'   (zero for every patient) are dropped rather than returned as all-zero
#'   columns.
#' @export
calculateDirect <- function(kernel, unifKernel = NULL, cells, Dcell = NULL, bundle = TRUE, patient_names = NULL, zero_threshold = 1) {
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

  df <- as.data.frame(df_list, check.names = FALSE)
  rownames(df) <- if (!is.null(patient_names)) patient_names else paste0("Patient_", 1:n_patients)

  # Drop cell-type pairs whose fraction of zero-valued patients reaches
  # zero_threshold -- default 1 keeps only the original behavior (drop a pair
  # only if it's zero for every single patient, a structurally-impossible
  # pathway; both branches above fill the impossible NaN case with 0, so this
  # is a correct and sufficient test regardless of raw vs normalized). Lower
  # zero_threshold (e.g. 0.9) to also drop merely zero-inflated pairs, which
  # this all-zero check alone does not catch.
  df <- df[, colMeans(df == 0) < zero_threshold, drop = FALSE]

  return(df)
}

#' Derive communication features from a kernel
#'
#' @param kernel Kernel array from [compute_kernel()].
#' @param unifKernel Optional normalized baseline kernel. When supplied,
#'   features are the enrichment ratio over it (see [calculateDirect()]).
#' @param celltypes Character vector of cell-type labels.
#' @param communication_type Feature family to compute. The kernel method only
#'   computes direct communication (`"D"`); GSCC is available with
#'   [compute_racing_montecarlo()].
#' @param bundle Logical; if `TRUE`, merge the two directions of each cell-type pair.
#' @param patient_names Optional patient labels.
#' @param Dcell Patient-by-cell-type abundance matrix. Required when
#'   `unifKernel` is not supplied (raw, unnormalized features).
#' @param patient_idx Optional patient index subset.
#' @param zero_threshold Drop a feature once its fraction of zero-valued
#'   patients reaches this threshold. Default `1` only drops features that
#'   are zero for every patient; lower it (e.g. `0.9`) to also drop merely
#'   zero-inflated features. Passed through to [calculateDirect()].
#'
#' @return A patient-by-feature data frame of direct communication scores for
#'   the selected patients.
#' @export
compute_kernel_features <- function(kernel, unifKernel = NULL, celltypes, communication_type = "D", bundle = TRUE,
                                    patient_names = NULL, Dcell = NULL, patient_idx = NULL, zero_threshold = 1) {
  .check_kernel_communication_type(communication_type)

  # Subset patients if requested
  if (!is.null(patient_idx)) {
    kernel <- kernel[,,patient_idx, drop = FALSE]
    if (!is.null(unifKernel)) unifKernel <- unifKernel[,,patient_idx, drop = FALSE]
    if (!is.null(patient_names)) patient_names <- patient_names[patient_idx]
    if (!is.null(Dcell)) Dcell <- Dcell[patient_idx, , drop = FALSE]
  }

  if (is.null(Dcell) && is.null(unifKernel)) {
    stop("Dcell (cell abundance matrix) is required to compute raw (unnormalized) direct features.")
  }

  calculateDirect(kernel = kernel, unifKernel = unifKernel, cells = celltypes, Dcell = Dcell, bundle = bundle,
                  patient_names = patient_names, zero_threshold = zero_threshold)
}

#' Check that only direct communication is requested from the kernel method
#'
#' @param communication_type Requested feature family.
#'
#' @return Invisibly `TRUE`; errors for anything other than `"D"`.
#' @keywords internal
.check_kernel_communication_type <- function(communication_type) {
  if (!identical(toupper(communication_type), "D")) {
    stop("The kernel method only computes direct communication (communication_type = \"D\"). ",
         "GSCC is available with compute_racing_montecarlo(); wedges and triangles were removed from RaCInG.",
         call. = FALSE)
  }
  invisible(TRUE)
}

#' Run the full kernel-based RaCInG workflow
#'
#' @param counts Gene-by-sample count matrix. Required when `input_data` is not
#'   supplied; ignored otherwise.
#' @param output_folder Directory used to write and read intermediate input files.
#' @param already_normalized Logical; if `TRUE`, `counts` is taken as linear-scale
#'   TPM and not renormalized. See [prepare_input_files()].
#' @param deconv Optional patient-by-cell-type abundance matrix. If omitted, it
#'   is computed via `multideconv::compute.deconvolution()` followed by
#'   `multideconv::compute.deconvolution.analysis()` (which identifies and
#'   collapses correlated cell-type subgroups; requires multideconv >= 0.2.0).
#'   Columns keep multideconv's names: subgroups (`"B.cells_Subgroup.1"`)
#'   and unclustered features (`"Quantiseq_TIL10_B.cells"`). See [prepare_input_files()].
#' @param cc_network Optional ligand-receptor prior network. If omitted, it is
#'   retrieved via `liana::get_curated_omni()`. See [prepare_input_files()].
#' @param fun_LR Function combining a ligand's and receptor's expression
#'   values (as a length-2 vector) into one interaction strength. Default
#'   `min`: the weaker-expressed partner limits the interaction. Any function
#'   collapsing 2 values to 1 works, e.g. `prod` for the product of the two
#'   instead.
#' @param cell_expr_profile Optional gene-by-cell-type expression profile
#'   matrix. If omitted, it is estimated from `counts` and `deconv` via
#'   per-gene non-negative least squares. See [prepare_input_files()].
#' @param expr_threshold Minimum expression (in linear TPM, not log-scale)
#'   a ligand or receptor must reach in its sending/receiving cell type to be
#'   kept. See [prepare_input_files()].
#' @param source,target Column names to use as ligand and receptor identifiers in `cc_network`.
#' @param signed Logical; if `TRUE`, also try to load a sign matrix.
#' @param deconv_method Deconvolution method(s) used when `deconv` is not supplied.
#' @param cbsx.name,cbsx.token Optional credentials for the deconvolution workflow.
#' @param file_name File stem used for intermediate input files.
#' @param nPatients Number of patients to process, or `"all"`.
#' @param communication_type Feature family to compute. The kernel method only
#'   computes direct communication (`"D"`, the default); GSCC is available with
#'   [compute_racing_montecarlo()].
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
#' @param zero_threshold Drop a feature once its fraction of zero-valued
#'   patients reaches this threshold. Default `1` only drops features that
#'   are zero for every patient (the original behavior); lower it (e.g.
#'   `0.9`) to also drop merely zero-inflated features. See
#'   [compute_kernel_features()].
#'
#' @return A list with three elements:
#'   \describe{
#'     \item{`kernel`}{The raw kernel array `[sender x receiver x patient]`
#'       from [compute_kernel()].}
#'     \item{`kernel_norm`}{The uniformized baseline kernel when `norm = TRUE`,
#'       otherwise `NULL`.}
#'     \item{`features`}{A patient-by-feature data frame of direct
#'       communication scores, with columns named `"Dir_A_B"` and rows named
#'       by patient.}
#'   }
#' @export
#' @examples
#' data(skcm_example)
#' res <- compute_racing_kernel(input_data = skcm_example)
#' dim(res$features)
#' head(res$features[, 1:4])
#'
#' # Enrichment ratio over a uniform ligand-receptor baseline
#' res_norm <- compute_racing_kernel(input_data = skcm_example, norm = TRUE)
#' head(res_norm$features[, 1:4])
compute_racing_kernel = function(counts = NULL, output_folder = "Results/", deconv = NULL, cc_network = NULL, fun_LR = min,
                                 cell_expr_profile = NULL, source = "source_genesymbol", target = "target_genesymbol", signed = FALSE,
                                 deconv_method = "Quantiseq", cbsx.name = NULL, cbsx.token = NULL, file_name = NULL, nPatients = "all",
                                 communication_type = "D", norm = FALSE, pt_idx = NULL, remove_direction = TRUE,
                                 input_data = NULL, zero_threshold = 1, expr_threshold = 10, already_normalized = FALSE) {

  .check_kernel_communication_type(communication_type)

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
                               signed = signed, already_normalized = already_normalized, expr_threshold = expr_threshold)
  }

  Lmatrix   <- res$Lmatrix
  Rmatrix   <- res$Rmatrix
  Cmatrix    <- res$Cmatrix
  LRmatrix   <- res$LRmatrix
  cellTypes <- res$celltypes
  # Name unnamed patients by their Cmatrix row, so subset runs keep the right labels
  if (is.null(rownames(Cmatrix))) rownames(Cmatrix) <- paste0("Patient_", seq_len(nrow(Cmatrix)))

  if (nPatients == "all") {
    nPatients <- nrow(Cmatrix)
  }

  if(!is.null(pt_idx)){
    cat("Because patient index is provided, nPatients argument will be ignored.\n")
  }else{
    cat("Computing kernel for", nPatients, "patients\n")
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
  patient_names <- rownames(Cmatrix) # the rows actually used (subset by nPatients above)
  features <- compute_kernel_features(
    kernel = kernel_raw,
    unifKernel = kernel_norm,
    celltypes = cellTypes,
    communication_type = communication_type,
    bundle = remove_direction,
    patient_names = patient_names,
    Dcell = Cmatrix,
    patient_idx = pt_idx,
    zero_threshold = zero_threshold
  )

  return(list(kernel = kernel_raw, kernel_norm = kernel_norm, features = features))

}