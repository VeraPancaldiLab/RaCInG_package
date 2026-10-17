#' Sort a character vector and return the permutation index
#'
#' Internal helper used to keep labels aligned after alphabetical reordering.
#'
#' @param stringlist Character vector to sort.
#'
#' @return A list with the sorted vector (`sortlist`) and its permutation index (`perm`).
#' @keywords internal
sortPermute <- function(stringlist) {
  ord <- order(stringlist)
  sortlist <- stringlist[ord]
  perm <- ord
  return(list(sortlist = sortlist, perm = perm))
}

#' Read a cell-to-ligand compatibility matrix
#'
#' @param filename Path to a CSV file with cell types in rows and ligands in columns.
#'
#' @return A list containing the matrix plus ligand and cell-type labels.
#' @export
createCellLigList <- function(filename) {
  if (!file.exists(filename)) {
    print(filename)
    stop("ERROR: The file with cell-ligand interactions does not exist...")
  }
  
  df <- utils::read.csv(filename, row.names = 1, check.names = FALSE)
  
  # Store the column names (ligands) and row names (cell types)
  ligands <- colnames(df)
  celltypes <- rownames(df)
  
  # Convert the data frame into a numeric matrix while removing rownames
  CellLigList = df %>%
    tibble::rownames_to_column("Celltype") %>%
    dplyr::select(-Celltype) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), as.numeric)) %>%
    as.matrix() 
  
  # Sort the cell types alphabetically and get the permutation indices
  sp <- sortPermute(celltypes)
  celltypes <- sp$sortlist # sorted cell type names
  perm <- sp$perm # permutation to reorder rows
  
  CellLigList <- CellLigList[perm, , drop = FALSE]
  
  return(list(CellLigList = CellLigList,
              ligands = ligands,
              celltypes = celltypes))
}

#' Read a cell-to-receptor compatibility matrix
#'
#' @param filename Path to a CSV file with cell types in rows and receptors in columns.
#'
#' @return A list containing the matrix plus receptor and cell-type labels.
#' @export
createCellRecList <- function(filename) {
  if (!file.exists(filename)) {
    stop("ERROR: The file with cell-receptor interactions does not exist...")
  }
  
  df <- utils::read.csv(filename, row.names = 1, check.names = FALSE)
  
  # Store the column names (receptors) and row names (cell types)
  receptors <- colnames(df)
  celltypes <- rownames(df)

  # Convert the data frame into a numeric matrix while removing rownames
  CellRecList = df %>%
    tibble::rownames_to_column("Celltype") %>%  # move rownames into a column
    dplyr::select(-Celltype) %>%                # remove the rowname column
    dplyr::mutate(dplyr::across(dplyr::everything(), as.numeric)) %>%  # convert all columns to numeric
    as.matrix()   

  # Sort the cell types alphabetically and get the permutation indices
  sp <- sortPermute(celltypes)
  celltypes <- sp$sortlist
  perm <- sp$perm
  
  CellRecList <- CellRecList[perm, , drop = FALSE]
  
  return(list(CellRecList = CellRecList,
              receptors = receptors,
              celltypes = celltypes))
}

#' Read and normalize cell-type abundance estimates
#'
#' @param cells Character vector of expected cell types.
#' @param filename Path to a CSV file containing patient-by-cell-type abundances.
#'
#' @return A list with the normalized `Dtypes` matrix plus labels.
#' @export
createCellTypeDistr <- function(cells, filename) {
  if (!file.exists(filename)) {
    stop("Cell type distribution file does not exist...")
  }
  
  # Read CSV file, using the first column as rownames (datasets)
  df <- utils::read.csv(filename, row.names = 1, check.names = FALSE)
  
  # Extract column names as cell type names and row names as dataset/patient names
  celltypes <- colnames(df)
  datasets <- rownames(df)

  # Convert the data frame to a numeric matrix, keeping dataset names as rownames
  # so downstream functions (e.g. compute_racing_montecarlo) can recover patient IDs
  Dtypes <- as.matrix(df)
  storage.mode(Dtypes) <- "numeric"
  rownames(Dtypes) <- datasets

  # Normalize rows so each dataset sums to 1 (relative fractions)
  row_sums <- rowSums(Dtypes)
  normalized <- Dtypes / row_sums

  # Sort cell types alphabetically
  sp <- sortPermute(celltypes)
  celltypes <- sp$sortlist
  perm <- sp$perm
  Dtypes <- normalized[, perm, drop = FALSE]

  # Restrict to the cell types that actually appear in Lmatrix/Rmatrix.
  # ccc_table's expression threshold can leave some cell types with no
  # ligand/receptor pair at all (they neither send nor receive anything above
  # threshold), so `cells` (from Lmatrix) can be a strict subset of the
  # deconvolution's full cell-type vocabulary. Those cell types contribute
  # nothing to the network, so they're dropped from Cmatrix here rather than
  # left in as columns with no matching Lmatrix/Rmatrix row.
  if (!identical(cells, celltypes)) {
    if (!all(cells %in% celltypes)) {
      stop("Cell types in Lmatrix/Rmatrix are not a subset of Cmatrix's cell types -- mismatched input files.")
    }
    keep <- celltypes %in% cells
    Dtypes <- Dtypes[, keep, drop = FALSE]
    celltypes <- celltypes[keep]
    Dtypes <- Dtypes / rowSums(Dtypes) # re-normalize so each dataset still sums to 1 after dropping cell types with no active ligand/receptor
  }

  return(list(Dtypes = Dtypes,
              celltypes = celltypes,
              datasets = datasets))
}

#' Read ligand-receptor interaction probabilities
#'
#' @param filename Path to a CSV file containing patient-by-interaction weights.
#' @param ligands Character vector of ligand names.
#' @param receptors Character vector of receptor names.
#'
#' @return A 3D array with dimensions ligand × receptor × patient.
#' @export
createInteractionDistr <- function(filename, ligands, receptors) {
  if (!file.exists(filename)) {
    stop("Interaction distribution file does not exist...")
  }
  
  df <- utils::read.csv(filename, row.names = 1, check.names = FALSE)
  
  # Store column names (interaction names like LIG_REC)
  interactions <- colnames(df)
  
  # Convert data frame to numeric matrix, replacing NA values with 0
  distr <- df %>%
    tibble::rownames_to_column("Datasets") %>%
    dplyr::select(-Datasets) %>%
    dplyr::mutate(dplyr::across(
      dplyr::everything(),                      # apply to all columns
      ~ tidyr::replace_na(as.numeric(.), 0)     # convert to numeric and replace NA with 0
    )) %>%
    as.matrix() 
  
  # Split interaction names "LIG_REC" into ligand and receptor parts
  parts <- strsplit(interactions, "_")
  rowWordIndex <- sapply(parts, `[`, 1)
  colWordIndex <- sapply(parts, `[`, 2)
  
  # Convert ligand/receptor names into numeric indices matching the provided ligand/receptor lists
  row <- match(rowWordIndex, ligands)
  col <- match(colWordIndex, receptors)
  
  # Create a 3D tensor to hold ligand-receptor probabilities
  # Dimensions: ligands x receptors x datasets
  DconnectionTensor <- array(0, dim = c(length(ligands), length(receptors), nrow(distr)))
  
  # Fill the tensor: for each dataset, assign probability values to the correct ligand-receptor pairs
  for (i in 1:nrow(distr)) {
    DconnectionTensor[cbind(row, col, i)] <- distr[i, ]
  }
  
  # Normalize each dataset slice so that the sum of all probabilities equals 1
  # This ensures that each dataset represents a valid probability distribution over all ligand-receptor pairs
  normsum <- apply(DconnectionTensor, 3, sum)
  
  for (i in 1:length(normsum)) {
    if (normsum[i] != 0) { # avoids division by 0
      DconnectionTensor[, , i] <- DconnectionTensor[, , i] / normsum[i]
    }
  }
  
  # Return the normalized 3D tensor
  return(DconnectionTensor)
}

#' Read a ligand-receptor sign matrix
#'
#' @param filename Path to a delimited file containing ligand-receptor signs or weights.
#'
#' @return A list with the numeric sign matrix and its ligand/receptor labels.
#' @export
Read_Lig_Rec_Interaction <- function(filename) {
  if (!file.exists(filename)) {
    stop("Interaction sign file does not exist...")
  }
  
  df <- data.table::fread(filename, header = TRUE)
  
  receptor_names <- colnames(df)[-1]
  ligand_names <- df[[1]]
  
  sign_matrix <- as.matrix(df[, -1])
  storage.mode(sign_matrix) <- "numeric"
  
  return(list(sign_matrix = sign_matrix,
              ligand_names = ligand_names,
              receptor_names = receptor_names))
}

#' Load RaCInG input matrices from disk
#'
#' Internal helper called by [prepare_input_files()] to read back the generated
#' CSV files and assemble the normalised matrices and 3-D tensor. Normalisation
#' happens in the readers it calls: [createCellTypeDistr()] rescales each
#' patient's cell-type abundances to sum to 1, and [createInteractionDistr()]
#' rescales each patient's ligand-receptor weights to sum to 1.
#' Users should call [prepare_input_files()] instead.
#'
#' @param file_name File stem used for the `Lmatrix_`, `Rmatrix_`, `Cmatrix_`, and
#'   `LRmatrix_` CSV files.
#' @param output_folder Directory containing the exported input files.
#' @param read_signs Logical; if `TRUE`, attempts to read `Sign_matrix_<file_name>.csv`
#'   from `output_folder`.
#'
#' @return A named list containing the input matrices and their labels.
#' @keywords internal
generateInput <- function(file_name, output_folder, read_signs = FALSE) {
  
  # -----------------------------
  # Read cell-to-ligand compatibility information
  # The function returns:
  #   CellLigList: numeric matrix indicating which cell type can interact with which ligand
  #   ligands: vector of ligand names
  #   celltypes: vector of cell type names (alphabetically sorted)
  # -----------------------------
  lig_path <- file.path(output_folder, paste0("Lmatrix_", file_name, ".csv"))
  lig_data <- createCellLigList(lig_path)
  
  # -----------------------------
  # Read cell-to-receptor compatibility information
  # The function returns:
  #   CellRecList: numeric matrix indicating which cell type can interact with which receptor
  #   receptors: vector of receptor names
  #   celltypes: vector of cell type names (alphabetically sorted)
  # -----------------------------
  rec_path <- file.path(output_folder, paste0("Rmatrix_", file_name, ".csv"))
  rec_data <- createCellRecList(rec_path)

  # Rmatrix rows are used with Lmatrix's cell-type labels below, so both files
  # must list the same cell types (both readers sort them alphabetically).
  if (!identical(lig_data$celltypes, rec_data$celltypes)) {
    stop("Lmatrix and Rmatrix must list the same cell types. Only in Lmatrix: ",
         paste(setdiff(lig_data$celltypes, rec_data$celltypes), collapse = ", "),
         "; only in Rmatrix: ", paste(setdiff(rec_data$celltypes, lig_data$celltypes), collapse = ", "), ".",
         call. = FALSE)
  }
  
  # -----------------------------
  # Read cell type distribution per dataset/patient
  # The function returns:
  #   Dtypes: matrix with relative abundance of each cell type per dataset (rows normalized)
  #   celltypes: vector of cell types
  #   datasets: vector of dataset/patient names
  # -----------------------------
  dist_path <- file.path(output_folder, paste0("Cmatrix_", file_name, ".csv"))
  dist_data <- createCellTypeDistr(
    cells = lig_data$celltypes,
    filename = dist_path
  )
  
  # -----------------------------
  # Read ligand-receptor interaction weights
  # The function returns:
  #   DconnectionTensor: 3D array (ligands x receptors x datasets) of interaction probabilities
  # -----------------------------
  lr_path <- file.path(output_folder, paste0("LRmatrix_", file_name, ".csv"))
  DconnectionTensor <- createInteractionDistr(
    filename = lr_path,
    ligands = lig_data$ligands,
    receptors = rec_data$receptors
  )
  
  # -----------------------------
  # Optionally read the sign of interactions (+1 stimulating, -1 inhibiting, 0 unknown)
  # If no sign matrix is supplied, fall back to zeros (unknown sign).
  # -----------------------------
  if (isTRUE(read_signs)) {
    sign_path <- file.path(output_folder, paste0("Sign_matrix_", file_name, ".csv"))
    if (!file.exists(sign_path)) {
      stop(
        "`read_signs = TRUE` requires a sign matrix file at ", sign_path,
        call. = FALSE
      )
    }
    sign_data <- Read_Lig_Rec_Interaction(sign_path)
    Sign_matrix <- sign_data$sign_matrix
  } else {
    Sign_matrix <- matrix(
      0,
      nrow = dim(DconnectionTensor)[1],
      ncol = dim(DconnectionTensor)[2]
    )
  }
  
  # -----------------------------
  # Return all input data in a list for use in the network generation model
  # -----------------------------
  return(list(
    Lmatrix = lig_data$CellLigList,       # cell-ligand compatibility matrix
    Rmatrix = rec_data$CellRecList,       # cell-receptor compatibility matrix
    Cmatrix = dist_data$Dtypes,                # normalized cell type distributions
    LRmatrix = DconnectionTensor,    # ligand-receptor interaction tensor
    celltypes = lig_data$celltypes,           # vector of cell type names
    ligands = lig_data$ligands,               # vector of ligand names
    receptors = rec_data$receptors,           # vector of receptor names
    Sign_matrix = Sign_matrix                  # matrix of interaction signs
  ))
}

#' Warn about cell-type names outside multideconv's nomenclature
#'
#' A name is accepted if, after dropping a `"_Subgroup.<n>"` suffix, it is a
#' multideconv cell type or ends in `"_<cell type>"` (e.g.
#' `"Quantiseq_TIL10_B.cells"`).
#'
#' @param celltypes Character vector of cell-type (column) names.
#'
#' @return Invisibly, the non-standard names.
#' @keywords internal
.warn_nonstandard_celltypes <- function(celltypes) {
  vocab <- multideconv::get_cell_type_nomenclature()
  base <- sub("_Subgroup\\.[0-9]+$", "", celltypes)
  ok <- vapply(base, function(x) any(x == vocab | endsWith(x, paste0("_", vocab))), logical(1))
  bad <- celltypes[!ok]
  if (length(bad) > 0) {
    warning("These cell-type names in `deconv` are not in multideconv's standardized nomenclature: ",
            paste(bad, collapse = ", "), ". Please name cell types with ",
            "multideconv::get_cell_type_nomenclature() (multideconv::standardize_celltype_colnames() ",
            "can help), so they are consistent across methods. Ignore this for cell types the ",
            "nomenclature does not cover.", call. = FALSE)
  }
  invisible(bad)
}

#' Estimate per-cell-type average expression profiles from one composition
#'
#' For each gene, solves a non-negative least squares regression of bulk TPM
#' expression across samples against per-sample cell-type fractions, giving
#' one expression value per gene per cell type:
#' `bulk_expr[gene, sample] ~= sum_celltype fractions[sample, celltype] * profile[gene, celltype]`.
#' Valid only when every column of `deconv` already belongs to one coherent
#' composition (fractions summing to ~1 per sample). When `deconv` instead
#' spans multiple deconvolution methods/signatures concatenated together (as
#' `multideconv::compute.deconvolution()` produces before subgroup analysis),
#' use [.estimate_expression_profiles_subgroups()] instead -- fitting NNLS
#' directly across method/signature blocks here would violate the mixture
#' equation this regression assumes.
#'
#' @param counts_tpm Gene-by-sample TPM matrix.
#' @param deconv Sample-by-cell-type fraction matrix.
#'
#' @return A gene-by-cell-type data frame of estimated expression profiles.
#' @keywords internal
.estimate_expression_profiles_samples <- function(counts_tpm, deconv) {
  samples <- intersect(colnames(counts_tpm), rownames(deconv))
  X <- as.matrix(deconv[samples, , drop = FALSE])
  Y <- as.matrix(counts_tpm[, samples, drop = FALSE])

  profile <- matrix(0, nrow = nrow(Y), ncol = ncol(X),
                    dimnames = list(rownames(Y), colnames(X)))
  for (g in seq_len(nrow(Y))) {
    profile[g, ] <- stats::coef(nnls::nnls(X, Y[g, ]))
  }

  as.data.frame(profile)
}

#' Estimate per-cell-type expression profiles across multiple deconvolution
#' methods/signatures
#'
#' `multideconv::compute.deconvolution()`'s raw output concatenates one block
#' of columns per method-signature combination (named
#' `"{method}_{signature}_{celltype}"`). Each block is, on its own, a valid
#' composition -- multideconv itself verifies every block's standardized
#' cell-type columns sum to 1 per sample -- but
#' `multideconv::compute.deconvolution.analysis()`'s merged "Deconvolution
#' matrix" is not: it stitches together representatives from different blocks
#' per cell type, so fitting [.estimate_expression_profiles_samples()]'s NNLS
#' mixture equation directly on it is invalid. This instead fits one profile
#' per method-signature block (always a valid composition), then merges the
#' resulting *profiles* -- never the raw fractions -- according to whichever
#' subgroups `subgroups` (already computed by the caller via
#' `multideconv::compute.deconvolution.analysis()`) forms across those blocks.
#'
#' Requires multideconv >= 0.2.0, whose subgroups come from a single
#' complete-linkage clustering (members are always raw `deconv_wide` columns,
#' never other subgroups) and keep feature names unchanged, so subgroup members
#' and surviving raw features match `colnames(deconv_wide)` as-is. Each
#' subgroup's profile is the median of its members' profiles, mirroring how
#' multideconv builds the subgroup fraction itself.
#'
#' @param counts_tpm Gene-by-sample TPM matrix.
#' @param deconv_wide Sample-by-`{method}_{signature}_{celltype}` raw fraction
#'   matrix -- `multideconv::compute.deconvolution()`'s direct output, before
#'   subgroup analysis.
#' @param subgroups The list returned by
#'   `multideconv::compute.deconvolution.analysis(deconv_wide)`.
#'
#' @return A gene-by-cell-type data frame of estimated expression profiles,
#'   one column per feature surviving in `subgroups[["Deconvolution matrix"]]`.
#' @keywords internal
.estimate_expression_profiles_subgroups <- function(counts_tpm, deconv_wide, subgroups) {
  block_id <- sub("^([^_]+_[^_]+)_.+$", "\\1", colnames(deconv_wide))

  profiles_by_block <- lapply(split(colnames(deconv_wide), block_id), function(cols) {
    as.matrix(.estimate_expression_profiles_samples(counts_tpm, deconv_wide[, cols, drop = FALSE]))
  })

  common_genes <- Reduce(intersect, lapply(profiles_by_block, rownames))
  all_profiles <- do.call(cbind, lapply(profiles_by_block, function(p) p[common_genes, , drop = FALSE]))

  dt <- subgroups[["Deconvolution matrix"]]
  # One flat subgroup-name -> member-columns map across all cell types.
  subgroup_lookup <- unlist(unname(subgroups[["Deconvolution subgroups composition"]]), recursive = FALSE)

  merged <- sapply(colnames(dt), function(feature) {
    raw_cols <- if (feature %in% names(subgroup_lookup)) {
      subgroup_lookup[[feature]]
    } else {
      feature
    }
    missing_cols <- setdiff(raw_cols, colnames(all_profiles))
    if (length(missing_cols) > 0) {
      stop("Deconvolution feature '", feature, "' could not be traced back to its ",
           "method-signature column(s) (", paste(missing_cols, collapse = ", "), ") -- ",
           "naming mismatch between compute.deconvolution.analysis() output and ",
           "deconv_wide. RaCInG requires multideconv >= 0.2.0.")
    }
    apply(as.matrix(all_profiles[, raw_cols, drop = FALSE]), 1, stats::median)
  })
  dimnames(merged) <- list(common_genes, colnames(dt))

  stopifnot(ncol(merged) == ncol(dt))
  as.data.frame(merged)
}

#' Build RaCInG input files from raw count data
#'
#' This function combines the preprocessing and input-loading steps into a single
#' call. It generates the `L`, `R`, `C`, and `LR` CSV files from raw counts,
#' then reads them back to produce the normalised matrices and 3-D tensor required
#' by the kernel and Monte Carlo workflows.
#'
#' @param counts Gene-by-sample count matrix.
#' @param output_folder Directory where the generated `L`, `R`, `C`, and `LR` files are written.
#' @param deconv Optional patient-by-cell-type abundance matrix. If omitted, it
#'   is computed via `multideconv::compute.deconvolution()` followed by
#'   `multideconv::compute.deconvolution.analysis()` (which identifies and
#'   collapses correlated cell-type subgroups; requires multideconv >= 0.2.0).
#'   Columns keep multideconv's names: subgroups (`"B.cells_Subgroup.1"`)
#'   and unclustered features (`"Quantiseq_TIL10_B.cells"`). When you supply
#'   `deconv` yourself, its rows must be exactly `colnames(counts)`, in the same
#'   order, and its columns should be named with multideconv's standardized
#'   cell-type nomenclature (`multideconv::get_cell_type_nomenclature()`); a
#'   warning lists any that aren't.
#' @param cc_network Optional ligand-receptor prior network. If omitted, it is
#'   retrieved via `liana::get_curated_omni()` (a curated, OmniPath-backed
#'   consensus of CellPhoneDB/CellChatDB/ICELLNET/connectomeDB2020/CellTalkDB)
#'   the first time this is called, then cached indefinitely as an RDS file
#'   under `tools::R_user_dir("RaCInG", "cache")` and loaded from there on every
#'   subsequent call -- avoiding a live network round-trip (and its associated
#'   transient-failure risk) on every invocation. Delete that cache file to
#'   force a refresh from OmniPath.
#' @param fun_LR Function combining a ligand's and receptor's expression
#'   values (as a length-2 vector) into one interaction strength. Default
#'   `min`: the weaker-expressed partner limits the interaction. Any function
#'   collapsing 2 values to 1 works, e.g. `prod` for the product of the two
#'   instead.
#' @param cell_expr_profile Optional gene-by-cell-type expression profile
#'   matrix. If omitted, it is estimated via per-gene non-negative least
#'   squares -- per method-signature block, then merged according to the
#'   multideconv subgroup analysis, when `deconv` is auto-computed (see
#'   `.estimate_expression_profiles_subgroups()`), or directly from the
#'   supplied `deconv` otherwise (see `.estimate_expression_profiles_samples()`).
#' @param source,target Column names to use as ligand and receptor identifiers when `cc_network` is supplied.
#' @param deconv_method Deconvolution method(s) passed to `multideconv::compute.deconvolution()`.
#' @param cbsx.name,cbsx.token Optional credentials forwarded to the deconvolution workflow.
#' @param file_name File stem used when exporting the generated CSV files.
#' @param signed Logical; if `TRUE`, also try to load a sign matrix from `output_folder`.
#'
#' @param already_normalized Logical; if `TRUE`, `counts` is treated as
#'   already TPM-normalized, **linear scale, not logged**, and the internal
#'   `ADImpute::NormalizeTPM()` gene-length/library-size normalization is
#'   skipped. Use this for datasets where only pre-normalized expression is
#'   available (e.g. public cohorts distributed as log2(TPM+1) matrices with
#'   no raw counts -- convert those back to linear scale, e.g. `2^x - 1`,
#'   before passing them here) -- applying `NormalizeTPM()` to already-
#'   normalized data would re-apply its gene-length correction on top of the
#'   original normalization, distorting relative gene expression levels.
#'   Emits a warning, since `expr_threshold` below is only meaningful if this
#'   linear-scale assumption holds -- comparing a log-scale value against it
#'   would silently apply the wrong cutoff.
#'
#' @param expr_threshold Minimum expression value, in linear-scale TPM (in
#'   `cell_expr_profile` units), a ligand or receptor must reach in its
#'   sending/receiving cell type for that ligand-receptor-celltype triple to
#'   be kept in `CC_table`.
#'
#' @return A named list with the processed input matrices and their labels:
#'   `Lmatrix`, `Rmatrix`, `Cmatrix` (normalised), `LRmatrix` (3-D tensor),
#'   `celltypes`, `ligands`, `receptors`, `Sign_matrix`, and `CC_table`.
#' @export
prepare_input_files <- function(counts, output_folder = "Results/", deconv = NULL, cc_network = NULL, fun_LR = min,
                                cell_expr_profile = NULL, source = "source_genesymbol", target = "target_genesymbol",
                                deconv_method = "Quantiseq", cbsx.name = NULL, cbsx.token = NULL, file_name = NULL,
                                signed = FALSE, already_normalized = FALSE, expr_threshold = 10){

  if (is.null(file_name)) {
    file_name <- "RaCInG_input"
  }

  .check_installed_packages(if (already_normalized) "multideconv" else c("ADImpute", "multideconv"))
  if (already_normalized) {
    warning(
      "already_normalized = TRUE: 'counts' is assumed to already be linear-scale TPM. ",
      "expr_threshold (default ", expr_threshold, ") is compared against this scale -- ",
      "if 'counts' is actually log-transformed (e.g. log2(TPM+1)), that comparison is ",
      "wrong; convert back to linear TPM first or set expr_threshold accordingly.",
      call. = FALSE
    )
    counts.tpm <- as.matrix(counts)
  } else {
    counts.tpm <- ADImpute::NormalizeTPM(counts)
  }
  counts.log.tpm <- log2(counts.tpm + 1)

  cat("Calculating deconvolution estimates...\n")

  ## C-matrix
  # deconv_wide/subgroups stay NULL on the user-supplied-deconv path below --
  # there's no known method/signature block structure to fit expression
  # profiles per block, so that path falls back to treating the whole matrix
  # as one composition (see the expression-profile step further down).
  deconv_wide <- NULL
  subgroups <- NULL
  if (is.null(deconv)) {
    deconv_wide <- multideconv::compute.deconvolution(
      counts.tpm,
      normalized = FALSE,
      methods = deconv_method,
      credentials.mail = cbsx.name,
      credentials.token = cbsx.token,
      file_name = file_name
    )

    ## Same samples as counts, in counts' column order (the LR matrix follows that order)
    if (!setequal(rownames(deconv_wide), colnames(counts))) {
      stop("Patient names in deconvolution estimates do not match those in the counts matrix.")
    }
    deconv_wide <- deconv_wide[colnames(counts), , drop = FALSE]

    cat("Identifying cell type subgroups...\n")
    subgroups <- multideconv::compute.deconvolution.analysis(deconv_wide)
    deconv <- subgroups[["Deconvolution matrix"]]
  } else {
    cat("Using provided deconvolution estimates...\n")
    # Cmatrix rows follow rownames(deconv) and the LR matrix follows colnames(counts),
    # and they are later paired by position -- so they must be identical, in order.
    if (!identical(rownames(deconv), colnames(counts))) {
      if (!is.null(rownames(deconv)) && setequal(rownames(deconv), colnames(counts))) {
        stop("rownames(deconv) contain the same samples as colnames(counts) but in a different order; ",
             "reorder it first, e.g. deconv[colnames(counts), , drop = FALSE].", call. = FALSE)
      }
      stop("rownames(deconv) must be exactly the samples in colnames(counts), in the same order ",
           "(missing from deconv: ", length(setdiff(colnames(counts), rownames(deconv))),
           ", not in counts: ", length(setdiff(rownames(deconv), colnames(counts))), ").", call. = FALSE)
    }
    .warn_nonstandard_celltypes(colnames(deconv))
  }

  ## Cell type expression profiles
  cat("\nEstimating cell type expression profiles...\n")
  if (is.null(cell_expr_profile)) {
    .check_installed_packages("nnls")
    if (!is.null(deconv_wide)) {
      # Fit one profile per method-signature block (each a valid composition),
      # then merge according to whichever subgroups multideconv formed --
      # never fit NNLS directly on the merged Deconvolution matrix (`deconv`),
      # which stitches together representatives from different blocks and is
      # not itself a valid composition. See .estimate_expression_profiles_subgroups().
      cell_expr_profile <- .estimate_expression_profiles_subgroups(counts.tpm, deconv_wide, subgroups)
    } else {
      # A user-supplied deconv has no known block structure -- treat it as one
      # composition, as before.
      cell_expr_profile <- .estimate_expression_profiles_samples(counts.tpm, deconv)
    }
  }

  ## Prior knowledge network CC
  cat("Processing cell-cell interaction network...\n")
  if (is.null(cc_network)) {
    .check_installed_packages(c("OmnipathR", "liana"))
    cache_dir <- tools::R_user_dir("RaCInG", "cache")
    cache_file <- file.path(cache_dir, "curated_omni_network.rds")
    if (file.exists(cache_file)) {
      cc_network <- readRDS(cache_file)
    } else {
      # get_curated_omni()'s underlying calls hit OmnipathR's own on-disk cache for
      # most of the network data, but at least one internal step (parameter
      # validation against OmnipathR::resources()) always makes a live request
      # regardless -- so this whole result is still subject to transient network
      # failures on every fresh R session unless cached at this level too.
      cc_network <- liana::get_curated_omni() %>%
        dplyr::filter(category_intercell_source %in% c("cell_surface_ligand", "ligand") &
                      category_intercell_target %in% c("receptor", "adhesion")) %>%
        liana::decomplexify(columns = c("source_genesymbol", "target_genesymbol"))
      dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
      saveRDS(cc_network, cache_file)
    }
  } else {
    cc_network <- cc_network %>%
      dplyr::mutate(source_genesymbol = .data[[source]], target_genesymbol = .data[[target]])
  }
  
  ## Subset expression for LR pairs (keep only interactions present in cell_expr_profile)
  keep <- integer()
  for (i in seq_len(nrow(cc_network))) {
    ligand  <- cc_network$source_genesymbol[i]
    receptor <- cc_network$target_genesymbol[i]
    if (ligand %in% rownames(cell_expr_profile) & receptor %in% rownames(cell_expr_profile)) {
      keep <- c(keep, i)
    }
  }
  cc_interations_sub <- cc_network[keep, , drop = FALSE]
  
  ## Build CC table (Sender x each LR pair x Receiver)
  ## Vectorized and processed one sender at a time: the previous implementation
  ## built one data.frame per (cell_type, LR pair) via nested
  ## do.call(rbind, lapply(...)) calls, i.e. O(n_celltypes * n_LR_pairs) tiny
  ## data.frames rbound together. With many deconvolution features (e.g.
  ## combining several deconvolution methods) n_celltypes can reach 100+, and
  ## with a curated LR network of several thousand pairs this produced
  ## hundreds of thousands of rbind pieces and tens of millions of rows, slow
  ## enough (tens of minutes) and memory-heavy enough to crash on typical
  ## machines. A single vectorized pass across the full sender x LR x receiver
  ## grid fixed the speed but still peaks at O(n_celltypes^2 * n_LR_pairs)
  ## elements transiently (~15GB observed with 122 cell types, i.e. the
  ## machine's entire RAM) before the expression-threshold filter shrinks it
  ## back down. Building one sender's full receiver grid at a time and
  ## filtering immediately keeps peak memory at O(n_celltypes * n_LR_pairs)
  ## instead, independent of how many senders there are.
  celltypes_all <- colnames(cell_expr_profile)
  n_ct <- length(celltypes_all)
  n_lr <- nrow(cc_interations_sub)
  cell_expr_mat <- as.matrix(cell_expr_profile)
  L_all <- cc_interations_sub$source_genesymbol
  R_all <- cc_interations_sub$target_genesymbol
  L_row <- match(L_all, rownames(cell_expr_mat))
  R_row <- match(R_all, rownames(cell_expr_mat))

  lr_idx_rep       <- rep(seq_len(n_lr), each = n_ct)
  receiver_idx_rep <- rep(seq_len(n_ct), times = n_lr)
  R_row_rep        <- R_row[lr_idx_rep]

  ccc_table <- do.call(rbind, lapply(seq_len(n_ct), function(s) {
    Expr_Ligand   <- cell_expr_mat[L_row[lr_idx_rep], s]
    Expr_Receptor <- cell_expr_mat[cbind(R_row_rep, receiver_idx_rep)]
    keep <- Expr_Ligand >= expr_threshold & Expr_Receptor >= expr_threshold ## Keep only pairs expressed above threshold
    data.frame(
      Sender = rep(celltypes_all[s], sum(keep)),
      Ligand = L_all[lr_idx_rep[keep]],
      Expr_Ligand = Expr_Ligand[keep],
      Receiver = celltypes_all[receiver_idx_rep[keep]],
      Receptor = R_all[lr_idx_rep[keep]],
      Expr_Receptor = Expr_Receptor[keep]
    )
  }))

  ## Use the actual matched (Ligand, Receptor) pairs present in ccc_table --
  ## previously this zipped unique(Ligand) with unique(Receptor) positionally,
  ## silently mismatching real pairs whenever the two vectors had different
  ## lengths (recycling warning: "number of columns of result is not a
  ## multiple of vector length").
  lr_pairs_df <- unique(ccc_table[, c("Ligand", "Receptor")])
  ligands <- lr_pairs_df$Ligand
  receptors <- lr_pairs_df$Receptor
  LR_pairs <- paste0(ligands, "_", receptors)

  ## LR-matrix: compute fun(L, R) per sample using counts (rows=genes, cols=samples): Default min(L,R)
  cat("Computing LR matrix...\n")
  pos_lr <- match(c(rbind(ligands, receptors)), rownames(counts.log.tpm))
  pos_lr <- matrix(pos_lr, nrow = 2)
  
  LR_matrix <- apply(pos_lr, 2, function(idx) {
    if (any(is.na(idx))) {
      return(rep(NA, ncol(counts.log.tpm)))
    }
    apply(counts.log.tpm[idx, , drop = FALSE], 2, function(x) fun_LR(x))
  })
  colnames(LR_matrix) <- LR_pairs
  
  ## L-matrix: ligands vs cell types compatibility (from filtered cc table)
  cat("Computing L-matrix...\n")
  # Every cell type that sends *or* receives at least one interaction; a
  # receive-only cell type keeps an all-zero Lmatrix row.
  celltypes <- unique(c(ccc_table$Sender, ccc_table$Receiver))
  ligs <- unique(ccc_table$Ligand)
  Lmatrix <- matrix(0, nrow = length(celltypes), ncol = length(ligs),
                    dimnames = list(celltypes, ligs))
  for (ct in celltypes) {
    tmp <- ccc_table %>% dplyr::filter(Sender == ct) %>% dplyr::select(Ligand) %>% unique()
    Lmatrix[ct, tmp$Ligand] <- 1
  }
  
  ## R-matrix: receptors vs cell types compatibility
  cat("Computing R-matrix...\n")
  recs <- unique(ccc_table$Receptor)
  Rmatrix <- matrix(0, nrow = length(celltypes), ncol = length(recs),
                    dimnames = list(celltypes, recs))
  for (ct in celltypes) {
    tmp <- ccc_table %>% dplyr::filter(Receiver == ct) %>% dplyr::select(Receptor) %>% unique()
    Rmatrix[ct, tmp$Receptor] <- 1
  }
  
  ## Output files
  cat("Exporting input files for RaCiNG...\n")
  if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)
  utils::write.csv(Lmatrix, file.path(output_folder, paste0("Lmatrix_", file_name, ".csv")))
  utils::write.csv(Rmatrix, file.path(output_folder, paste0("Rmatrix_", file_name, ".csv")))
  utils::write.csv(deconv, file.path(output_folder, paste0("Cmatrix_", file_name, ".csv")))
  utils::write.csv(LR_matrix, file.path(output_folder, paste0("LRmatrix_", file_name, ".csv")))
  utils::write.csv(ccc_table, file.path(output_folder, paste0("CC_table_", file_name, ".csv")), row.names = FALSE)
  
  cat("Input files generated and saved to:", output_folder, "\n")

  ## Read the written CSVs back through generateInput to obtain the
  ## normalised matrices and 3-D tensor expected by downstream functions.
  cat("Processing input matrices...\n")
  processed <- generateInput(file_name, output_folder = output_folder,
                             read_signs = signed)

  ## Append the CC table to the processed result
  processed$CC_table <- ccc_table

  return(processed)
}