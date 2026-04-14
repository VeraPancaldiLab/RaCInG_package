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
  
  df <- utils::read.csv(filename, row.names = 1)
  
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
  
  df <- utils::read.csv(filename, row.names = 1)
  
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
  df <- utils::read.csv(filename, row.names = 1)
  
  # Extract column names as cell type names and row names as dataset names
  celltypes <- colnames(df)
  datasets <- rownames(df)
  
  # Convert the data frame to a numeric matrix and remove rownames
  Dtypes = df %>%
    tibble::rownames_to_column("Datasets") %>%
    dplyr::select(-Datasets) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), as.numeric)) %>%
    as.matrix()
  
  # Normalize rows so each dataset sums to 1 (relative fractions) --> ############### Estimates from deconv should be always proportions (this should not be necessary unless there is MCP or XCell)
  row_sums <- rowSums(Dtypes)
  normalized <- Dtypes / row_sums
  
  # Sort cell types alphabetically
  sp <- sortPermute(celltypes)
  celltypes <- sp$sortlist
  perm <- sp$perm
  Dtypes <- normalized[, perm, drop = FALSE]
  
  # Merge M1 and M2 if needed --> This will break for other deconvolution results (need to be fix)
  if (length(cells) != length(celltypes)) {
    idx_M1 <- which(celltypes == "M1")
    idx_M2 <- which(celltypes == "M2")
    
    if (length(idx_M1) > 0 && length(idx_M2) > 0) {
      Dtypes[, idx_M1] <- Dtypes[, idx_M1] + Dtypes[, idx_M2]
      Dtypes <- Dtypes[, -idx_M2, drop = FALSE]
      celltypes[idx_M1] <- "M"
      celltypes <- celltypes[-idx_M2]
    }
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
  
  df <- utils::read.csv(filename, row.names = 1)
  
  # Store column names (interaction names like LIG_REC) and row names (datasets)
  interactions <- colnames(df)
  datasets <- rownames(df)
  
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
#' CSV files and assemble the normalised matrices and 3-D tensor.
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
  lig_path = paste0(output_folder, "Lmatrix_", file_name, ".csv")
  lig_data <- createCellLigList(lig_path)
  
  # -----------------------------
  # Read cell-to-receptor compatibility information
  # The function returns:
  #   CellRecList: numeric matrix indicating which cell type can interact with which receptor
  #   receptors: vector of receptor names
  #   celltypes: vector of cell type names (alphabetically sorted)
  # -----------------------------
  rec_path = paste0(output_folder, "Rmatrix_", file_name, ".csv")
  rec_data <- createCellRecList(rec_path)
  
  # -----------------------------
  # Read cell type distribution per dataset/patient
  # The function returns:
  #   Dtypes: matrix with relative abundance of each cell type per dataset (rows normalized)
  #   celltypes: vector of cell types
  #   datasets: vector of dataset/patient names
  # -----------------------------
  dist_path = paste0(output_folder, "Cmatrix_", file_name, ".csv")
  dist_data <- createCellTypeDistr(
    cells = lig_data$celltypes,
    filename = dist_path
  )
  
  # -----------------------------
  # Read ligand-receptor interaction weights
  # The function returns:
  #   DconnectionTensor: 3D array (ligands x receptors x datasets) of interaction probabilities
  # -----------------------------
  lr_path = paste0(output_folder, "LRmatrix_", file_name, ".csv")
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

#' Build RaCInG input files from raw count data
#'
#' This function combines the preprocessing and input-loading steps into a single
#' call. It generates the `L`, `R`, `C`, and `LR` CSV files from raw counts,
#' then reads them back to produce the normalised matrices and 3-D tensor required
#' by the kernel and Monte Carlo workflows.
#'
#' @param counts Gene-by-sample count matrix.
#' @param output_folder Directory where the generated `L`, `R`, `C`, and `LR` files are written.
#' @param deconv Optional deconvolution matrix. If omitted, the function will try to compute it.
#' @param cc_network Optional ligand-receptor prior network.
#' @param fun_LR Function used to combine ligand and receptor expression values.
#' @param cell_expr_profile Optional cell-type expression profile matrix.
#' @param source,target Column names to use as ligand and receptor identifiers when `cc_network` is supplied.
#' @param deconv_method Deconvolution method passed to `multideconv::compute.deconvolution()`.
#' @param cbsx.name,cbsx.token Optional credentials forwarded to the deconvolution workflow.
#' @param file_name File stem used when exporting the generated CSV files.
#' @param signed Logical; if `TRUE`, also try to load a sign matrix from `output_folder`.
#'
#' @return A named list with the processed input matrices and their labels:
#'   `Lmatrix`, `Rmatrix`, `Cmatrix` (normalised), `LRmatrix` (3-D tensor),
#'   `celltypes`, `ligands`, `receptors`, `Sign_matrix`, and `CC_table`.
#' @export
prepare_input_files <- function(counts, output_folder = "Results/", deconv = NULL, cc_network = NULL, fun_LR = min, 
                                cell_expr_profile = NULL, source = "source_genesymbol", target = "target_genesymbol",
                                deconv_method = "Quantiseq", cbsx.name = NULL, cbsx.token = NULL, file_name = NULL,
                                signed = FALSE){
  
  if (is.null(file_name)) {
    file_name <- "RaCInG_input"
  }

  .check_installed_packages(c("ADImpute", "multideconv"))
  counts.tpm <- ADImpute::NormalizeTPM(counts)
  counts.log.tpm <- log2(counts.tpm + 1)

  cat("Calculating deconvolution estimates...\n")
  ## C-matrix
  if (is.null(deconv)) {
    deconv <- multideconv::compute.deconvolution(
      counts.tpm,
      normalized = FALSE,
      methods = deconv_method,
      credentials.mail = cbsx.name,
      credentials.token = cbsx.token,
      file_name = file_name
    )
  } else {
    cat("Using provided deconvolution estimates...\n")
  }
  
  ## Cell type expression profiles
  cat("\nEstimating cell type expression profiles...\n")
  if (is.null(cell_expr_profile)) {
    cell_expr_profile <- multideconv::estimate_expression_profiles(counts, deconv)
    cell_expr_profile = rbind(sapply(cell_expr_profile, function(x) colMeans(x))) %>% as.data.frame()
  }

  ## Verify if patients names match between files
  if(!all(rownames(deconv) %in% colnames(counts))){
    stop("Patient names in deconvolution estimates do not match those in the counts matrix.")
  }  

  ## Prior knowledge network CC
  cat("Processing cell-cell interaction network...\n")
  if (is.null(cc_network)) {
    .check_installed_packages(c("OmnipathR", "liana"))
    cc_network <- OmnipathR::import_intercell_network(high_confidence = TRUE) %>%
      dplyr::filter(category_intercell_source %in% c("cell_surface_ligand", "ligand") &
                    category_intercell_target %in% c("receptor", "adhesion")) %>%
      liana::decomplexify(columns = c("source_genesymbol", "target_genesymbol"))
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
  ccc_table <- do.call(rbind, lapply(colnames(cell_expr_profile), function(cell_type) {
    do.call(rbind, lapply(seq_len(nrow(cc_interations_sub)), function(LR_pos) {
      L <- cc_interations_sub$source_genesymbol[LR_pos]
      R <- cc_interations_sub$target_genesymbol[LR_pos]
      data.frame(
        Sender = cell_type,
        Ligand = L,
        Expr_Ligand = as.numeric(cell_expr_profile[L, cell_type]),
        Receiver = colnames(cell_expr_profile),
        Receptor = R,
        Expr_Receptor = as.numeric(cell_expr_profile[R, ])
      )
    }))
  })) %>%
    dplyr::filter(Expr_Ligand >= 10 & Expr_Receptor >= 10) ## Keep only pairs expressed above threshold (TPM >= 10)
  
  receptors <- unique(ccc_table$Receptor)
  ligands <- unique(ccc_table$Ligand)
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
  rownames(LR_matrix) <- gsub(".", "-", rownames(LR_matrix), fixed = TRUE)
  
  ## L-matrix: ligands vs cell types compatibility (from filtered cc table)
  cat("Computing L-matrix...\n")
  celltypes <- unique(ccc_table$Sender)
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