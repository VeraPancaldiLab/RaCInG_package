#' Count wedges across Monte Carlo graph simulations
#'
#' @param Dcell Cell-type abundance vector for one patient.
#' @param Dconn Ligand-receptor probability matrix for one patient.
#' @param lig Cell-by-ligand compatibility matrix.
#' @param rec Cell-by-receptor compatibility matrix.
#' @param cellnames Character vector of cell-type names.
#' @param N Number of cells per simulated graph.
#' @param av Target average degree.
#' @param itNo Number of Monte Carlo iterations.
#'
#' @return A list of average and standard-deviation wedge counts.
#' @export
countWedges <- function(Dcell, Dconn, lig, rec, cellnames, N, av, itNo) {
  # Number of cell types
  nCells <- length(cellnames)
  
  # Create an array to hold counts of wedges by type for each Monte Carlo iteration
  triangletensor <- array(0, dim = c(nCells, nCells, nCells, itNo))
  
  # Vector to store total number of wedges in each iteration
  trianglecount <- numeric(itNo)
  
  for (i in 1:itNo) {
    # Generate a random graph for this iteration
    graph <- model1(N, av, lig, rec, Dcell, Dconn, genRandom = FALSE)
    V <- graph$V      # Vector of cell types per vertex
    E <- graph$E      # Edge list
    
    # Build adjacency matrix from edge list, excluding self-loops
    Adj <- EdgetoAdj_No_loop(E, length(V))

    # --- CALL THE EXISTING WEDGES FUNCTION ---
    wedge_result <- Wedges(Adj)   # <- this uses already defined function

    # Store total number of wedges
    trianglecount[i] <- wedge_result$NoWedges

    # Count wedges by cell type combination
    triangletensor[,,,i] <- Count_Types(wedge_result$Wedge_list, V, maxTypes = nCells)
  }
  
  # Compute averages and standard deviations across Monte Carlo iterations
  av_triag <- apply(triangletensor, c(1,2,3), mean)
  std_triag <- apply(triangletensor, c(1,2,3), stats::sd)
  av_count <- mean(trianglecount)
  std_count <- stats::sd(trianglecount)
  
  return(list(
    av_triag = av_triag, 
    std_triag = std_triag, 
    av_count = av_count, 
    std_count = std_count
  ))
}

#' Count trust triangles across Monte Carlo graph simulations
#'
#' @inheritParams countWedges
#'
#' @return A list of average and standard-deviation trust-triangle counts.
#' @export
countTrustTriangles <- function(Dcell, Dconn, lig, rec, cellnames, N, av, itNo) {
  # Number of cell types
  nCells <- length(cellnames)
  
  # Create a 4D array to store counts of trust triangles by type
  # Dimensions: [type of vertex 1, type of vertex 2, type of vertex 3, Monte Carlo iteration]
  triangletensor <- array(0, dim = c(nCells, nCells, nCells, itNo))
  
  # Vector to store total number of trust triangles in each iteration
  trianglecount <- numeric(itNo)
  
  # Loop over Monte Carlo iterations
  for (i in 1:itNo) {
    # Generate a random graph for this iteration
    graph <- model1(N, av, lig, rec, Dcell, Dconn, genRandom = FALSE)
    V <- graph$V  # Vector assigning cell type to each vertex
    E <- graph$E  # Edge list
    
    # Build adjacency matrix from edge list, excluding self-loops
    Adj <- EdgetoAdj_No_loop(E, length(V))

    # --- CALL EXISTING FUNCTION TO FIND TRUST TRIANGLES ---
    tri_result <- Trust_Triangles(Adj)

    # Store total number of trust triangles in this iteration
    trianglecount[i] <- tri_result$NoTriangles

    # Count the number of triangles by cell type combination
    triangletensor[,,,i] <- Count_Types(tri_result$Triangle_list, V, maxTypes = nCells)
  }
  
  # Compute average triangle counts per cell-type combination across iterations
  av_triag <- apply(triangletensor, c(1,2,3), mean)
  
  # Compute standard deviation of triangle counts per combination across iterations
  std_triag <- apply(triangletensor, c(1,2,3), stats::sd)
  
  # Compute average total number of trust triangles per iteration
  av_count <- mean(trianglecount)
  
  # Compute standard deviation of total triangles across iterations
  std_count <- stats::sd(trianglecount)
  
  # Return all results as a list
  return(list(
    av_triag = av_triag,      # Average triangle counts by type combination
    std_triag = std_triag,    # Standard deviation by type combination
    av_count = av_count,      # Average total trust triangles
    std_count = std_count     # Std of total trust triangles
  ))
}

#' Count cycle triangles across Monte Carlo graph simulations
#'
#' @inheritParams countWedges
#'
#' @return A list of average and standard-deviation cycle-triangle counts.
#' @export
countCycleTriangles <- function(Dcell, Dconn, lig, rec, cellnames, N, av, itNo) {
  # Number of different cell types
  nCells <- length(cellnames)
  
  # 4D array to store counts of cycle triangles by type for each iteration
  # Dimensions: [vertex1 type, vertex2 type, vertex3 type, Monte Carlo iteration]
  triangletensor <- array(0, dim = c(nCells, nCells, nCells, itNo))
  
  # Vector to store total number of cycle triangles in each iteration
  trianglecount <- numeric(itNo)
  
  # Loop over Monte Carlo iterations
  for (i in 1:itNo) {
    # Generate random graph for this iteration
    graph <- model1(N, av, lig, rec, Dcell, Dconn, genRandom = FALSE)
    V <- graph$V  # Vector: cell type for each vertex
    E <- graph$E  # Edge list of the graph
    
    # Convert edge list to adjacency matrix, excluding self-loops
    Adj <- EdgetoAdj_No_loop(E, length(V))

    # --- CALL EXISTING FUNCTION TO FIND CYCLE TRIANGLES ---
    tri_result <- Cycle_Triangles(Adj)

    # Store total number of cycle triangles in this iteration
    trianglecount[i] <- tri_result$NoTriangles

    # Count the number of triangles by cell-type combination
    triangletensor[,,,i] <- Count_Types(tri_result$Triangle_list, V, maxTypes = nCells)
  }
  
  # Compute average triangle counts per cell-type combination across iterations
  av_triag <- apply(triangletensor, c(1,2,3), mean)
  
  # Compute standard deviation of triangle counts per combination across iterations
  std_triag <- apply(triangletensor, c(1,2,3), stats::sd)
  
  # Compute average total number of cycle triangles per iteration
  av_count <- mean(trianglecount)
  
  # Compute standard deviation of total cycle triangles across iterations
  std_count <- stats::sd(trianglecount)
  
  # Return all results as a list
  return(list(
    av_triag = av_triag,      # Average counts by cell-type combination
    std_triag = std_triag,    # Std by cell-type combination
    av_count = av_count,      # Average total cycle triangles
    std_count = std_count     # Std of total cycle triangles
  ))
}

#' Count direct edges across Monte Carlo graph simulations
#'
#' @inheritParams countWedges
#'
#' @return A list of average and standard-deviation direct edge counts.
#' @export
countDirect <- function(Dcell, Dconn, lig, rec, cellnames, N, av, itNo) {
  # Number of different cell types
  nCells <- length(cellnames)
  
  # 3D array to store counts of direct edges by type for each iteration
  # Dimensions: [source cell type, target cell type, Monte Carlo iteration]
  directtensor <- array(0, dim = c(nCells, nCells, itNo))
  
  # Vector to store total number of edges in each iteration
  directCount <- numeric(itNo)
  
  # Loop over Monte Carlo iterations
  for (i in 1:itNo) {
    # Generate random graph for this iteration
    graph <- model1(N, av, lig, rec, Dcell, Dconn, genRandom = FALSE)
    V <- graph$V  # Vector: cell type for each vertex
    E <- graph$E  # Edge list of the graph
    
    # Convert edge list to adjacency matrix, excluding self-loops
    Adj <- EdgetoAdj_No_loop(E, length(V))

    # Total number of edges in this iteration (matches Python: includes self-loops)
    directCount[i] <- nrow(E)
    
    # Count edges between each pair of cell types
    for (t in 1:nCells) {          # Source cell type
      for (s in 1:nCells) {        # Target cell type
        # Extract rows corresponding to source type t
        temp <- Adj[V == t, , drop=FALSE]
        # Count how many edges go to target type s
        directtensor[t, s, i] <- sum(temp[, V == s, drop=FALSE])
      }
    }
  }
  
  # Compute average edges between cell types across iterations
  av_dir <- apply(directtensor, c(1,2), mean)
  
  # Compute standard deviation of edges between cell types across iterations
  std_dir <- apply(directtensor, c(1,2), stats::sd)
  
  # Compute average total number of edges across iterations
  av_count <- mean(directCount)
  
  # Compute standard deviation of total edges across iterations
  std_count <- stats::sd(directCount)
  
  # Return results as a list
  return(list(
    av_dir = av_dir,      # Average edges by cell-type combination
    std_dir = std_dir,    # Std of edges by cell-type combination
    av_count = av_count,  # Average total edges
    std_count = std_count # Std of total edges
  ))
}

#' Count GSCC contributions across Monte Carlo graph simulations
#'
#' @inheritParams countWedges
#'
#' @return A list of average and standard-deviation GSCC contributions.
#' @export
countGSCC <- function(Dcell, Dconn, lig, rec, cellnames, N, av, itNo) {
  # Number of different cell types
  nCells <- length(cellnames)
  
  # 2D array to store fractional contributions of each cell type to GSCC in each iteration
  # Rows = cell types, Columns = Monte Carlo iterations
  GSCCcounttensor <- array(0, dim = c(nCells, itNo))
  
  # Vector to store fractional size of the GSCC for each iteration
  GSCCcount <- numeric(itNo)
  
  # Loop over Monte Carlo iterations
  for (i in 1:itNo) {
    # Generate random graph for this iteration
    graph <- model1(N, av, lig, rec, Dcell, Dconn, genRandom = FALSE)
    V <- graph$V  # Vector: cell type for each vertex
    E <- graph$E  # Edge list of the graph
    
    # Convert edge list to adjacency matrix, excluding self-loops
    Adj <- EdgetoAdj_No_loop(E, length(V))

    # Compute the GSCC (giant strongly connected component)
    # Returns a vector of vertex indices in the GSCC
    gscc <- GSCC(Adj)
    
    # Fractional size of the GSCC (number of nodes in GSCC / total nodes)
    GSCCcount[i] <- length(gscc) / length(V)
    
    # For each cell type, compute the fraction of GSCC nodes of that type
    for (t in 1:nCells) {
      GSCCcounttensor[t, i] <- sum(V[gscc] == t) / length(V)
    }
  }
  
  # Average contribution of each cell type to the GSCC across iterations
  av_GSCC <- apply(GSCCcounttensor, 1, mean)
  
  # Standard deviation of each cell type's contribution across iterations
  std_GSCC <- apply(GSCCcounttensor, 1, stats::sd)
  
  # Average fractional size of the GSCC across iterations
  av_count <- mean(GSCCcount)
  
  # Standard deviation of the GSCC size across iterations
  std_count <- stats::sd(GSCCcount)
  
  # Return all results as a list
  return(list(
    av_GSCC = av_GSCC,    # Average contribution of each cell type
    std_GSCC = std_GSCC,  # Std of contribution of each cell type
    av_count = av_count,  # Average GSCC size
    std_count = std_count # Std of GSCC size
  ))
}

#' Run Monte Carlo simulations for one or more patients
#'
#' @param Lmatrix Cell-by-ligand compatibility matrix.
#' @param Rmatrix Cell-by-receptor compatibility matrix.
#' @param Cmatrix Patient-by-cell-type abundance matrix.
#' @param LRmatrix Ligand-receptor-by-patient tensor.
#' @param cells Character vector of cell-type names.
#' @param communication_type Feature family to simulate (`"D"`, `"W"`, `"TT"`, `"CT"`, or `"GSCC"`).
#' @param pats Number of patients to process, or `"all"`.
#' @param N Number of cells per graph.
#' @param itNo Number of Monte Carlo iterations.
#' @param av Target average degree.
#' @param output_folder Directory used to write the `.out` files.
#' @param file.name Output filename stem.
#' @param norm Logical; if `TRUE`, use a uniformized LR baseline.
#' @param patient_idx Optional single patient index to simulate.
#'
#' @return Invisibly writes the simulation outputs to disk.
#' @export
runSim <- function(Lmatrix, Rmatrix, Cmatrix, LRmatrix, cells, communication_type, pats = "all",
                   N = 10000, itNo = 100, av = 20, output_folder = NULL, file.name = NULL, norm = FALSE,
                   patient_idx = NULL) {
  
  if (!is.null(output_folder) && !dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
  }
  
  cellstring <- paste(cells, collapse = ",")
  
  # Normalization
  if (norm) {
    normvec <- 1 / apply(LRmatrix != 0, 3, sum)
    for (i in 1:dim(LRmatrix)[3]) {
      LRmatrix[,,i][LRmatrix[,,i] != 0] <- normvec[i]
    }
    filename <- paste0(output_folder, "/", file.name, "_norm.out")
  } else {
    filename <- paste0(output_folder, "/", file.name, ".out")
  }
  
  # Determine which patient(s) to process.
  # If patient_idx is provided, use that specific patient (not the same as pats = 1).
  if (!is.null(patient_idx)) {
    if (length(patient_idx) != 1 || patient_idx < 1 || patient_idx > nrow(Cmatrix)) {
      stop("patient_idx must be a single integer between 1 and nrow(Cmatrix)")
    }
    pat_seq <- patient_idx
  } else if (pats == "all") {
    pat_seq <- 1:nrow(Cmatrix)
  } else {
    pat_seq <- seq_len(pats)
  }
  
  con <- file(filename, open = "w")
  
  # Header
  writeLines(communication_type, con)
  writeLines(paste(nrow(Cmatrix), N, itNo, av, sep=","), con)
  
  # ------------------------------------------------------------
  # Loop over patients (pat_seq may be a single index or a sequence)
  # ------------------------------------------------------------
  for (pat in pat_seq) {
    
    writeLines(paste(pat, N, av, sep=","), con)
    writeLines(cellstring, con)
    
    CellD <- Cmatrix[pat, ]
    IntD  <- LRmatrix[,,pat]
    
    # Compute feature
    if (communication_type == "D") {
      res <- countDirect(CellD, IntD, Lmatrix, Rmatrix, cells, N, av, itNo)
      av_mat <- res$av_dir
      std_mat <- res$std_dir
      
    } else if (communication_type == "W") {
      res <- countWedges(CellD, IntD, Lmatrix, Rmatrix, cells, N, av, itNo)
      av_mat <- res$av_triag
      std_mat <- res$std_triag
      
    } else if (communication_type == "TT") {
      res <- countTrustTriangles(CellD, IntD, Lmatrix, Rmatrix, cells, N, av, itNo)
      av_mat <- res$av_triag
      std_mat <- res$std_triag
      
    } else if (communication_type == "GSCC") {
      res <- countGSCC(CellD, IntD, Lmatrix, Rmatrix, cells, N, av, itNo)
      av_vec <- res$av_GSCC
      std_vec <- res$std_GSCC
      
    } else {
      res <- countCycleTriangles(CellD, IntD, Lmatrix, Rmatrix, cells, N, av, itNo)
      av_mat <- res$av_triag
      std_mat <- res$std_triag
    }
    
    # Write counts
    writeLines(paste("Count", res$av_count, res$std_count, sep=","), con)
    
    # ------------------------------------------------------------
    # Write composition
    # ------------------------------------------------------------
    if (communication_type == "D") {
      
      writeLines("Composition - Average:", con)
      for (i in 1:nrow(av_mat)) {
        for (j in 1:ncol(av_mat)) {
          writeLines(paste(i, j, av_mat[i,j], sep=","), con)
        }
      }
      
      writeLines("Composition - Std:", con)
      for (i in 1:nrow(std_mat)) {
        for (j in 1:ncol(std_mat)) {
          writeLines(paste(i, j, std_mat[i,j], sep=","), con)
        }
      }
      
    } else if (communication_type == "GSCC") {
      
      writeLines("Composition - Average:", con)
      for (i in 1:length(av_vec)) {
        writeLines(paste(i, av_vec[i], sep=","), con)
      }
      
      writeLines("Composition - Std:", con)
      for (i in 1:length(std_vec)) {
        writeLines(paste(i, std_vec[i], sep=","), con)
      }
      
    } else {
      
      writeLines("Composition - Average:", con)
      for (i in 1:dim(av_mat)[1]) {
        for (j in 1:dim(av_mat)[2]) {
          for (k in 1:dim(av_mat)[3]) {
            writeLines(paste(i, j, k, av_mat[i,j,k], sep=","), con)
          }
        }
      }
      
      writeLines("Composition - Std:", con)
      for (i in 1:dim(std_mat)[1]) {
        for (j in 1:dim(std_mat)[2]) {
          for (k in 1:dim(std_mat)[3]) {
            writeLines(paste(i, j, k, std_mat[i,j,k], sep=","), con)
          }
        }
      }
    }
  }
  
  close(con)
}

#' Run the full Monte Carlo RaCInG workflow
#'
#' @param counts Gene-by-sample count matrix. Required when `input_data` is not
#'   supplied; ignored otherwise.
#' @param output_folder Directory used to write intermediate and output files.
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
#' @param pt_idx Optional single patient index to simulate.
#' @param file_name File stem used for intermediate files.
#' @param nPatients Number of patients to process, or `"all"`.
#' @param communication_type Feature family to simulate.
#' @param Ncells Number of cells per simulated graph.
#' @param Ngraphs Number of Monte Carlo iterations.
#' @param Ndegree Target average degree.
#' @param remove_direction Logical; if `TRUE`, merge directionally equivalent features.
#' @param norm Logical; if `TRUE`, also run a uniformized baseline simulation for normalization.
#' @param input_data Optional named list of pre-computed input matrices as returned
#'   by [prepare_input_files()].
#'   Must contain `Lmatrix`, `Rmatrix`, `Cmatrix`, `LRmatrix`, `celltypes`,
#'   `ligands`, and `receptors`.
#'   When supplied, the `counts` argument and all preprocessing parameters
#'   (`deconv`, `cc_network`, etc.) are ignored.
#'
#' @return A list with the generated inputs and processed feature matrices.
#' @export
compute_racing_montecarlo = function(counts = NULL, output_folder = "~/Documents/racing/vignettes/", deconv = NULL, cc_network = NULL, fun_LR = min, 
                                     cell_expr_profile = NULL, source = "source_genesymbol", target = "target_genesymbol", signed = FALSE,
                                     deconv_method = "Quantiseq", cbsx.name = NULL, cbsx.token = NULL, pt_idx = NULL, file_name = NULL,
                                     nPatients = "all", communication_type = "W", Ncells = 10000, Ngraphs = 100, Ndegree = 20, remove_direction = TRUE, norm = TRUE,
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

  if (!is.null(pt_idx)) {
    cat("Because patient index is provided, nPatients argument will be ignored.\n")
    cat("Running Monte Carlo simulation for patient index:", pt_idx, "\n")}
  else{
    cat("Running Monte Carlo simulation for ", ifelse(nPatients == "all", "all", nPatients), " patients\n")
  }


  runSim(
    Lmatrix = Lmatrix,
    Rmatrix = Rmatrix,
    Cmatrix = Cmatrix,
    LRmatrix = LRmatrix,
    cells = cellTypes,
    communication_type = communication_type,
    pats = nPatients,
    N = Ncells,
    itNo = Ngraphs,
    av = Ndegree,
    output_folder = output_folder,
    file.name = file_name,
    norm = FALSE,
    patient_idx = pt_idx
  )

  if (norm) {
    runSim(
      Lmatrix = Lmatrix,
      Rmatrix = Rmatrix,
      Cmatrix = Cmatrix,
      LRmatrix = LRmatrix,
      cells = cellTypes,
      communication_type = communication_type,
      pats = nPatients,
      N = Ncells,
      itNo = Ngraphs,
      av = Ndegree,
      output_folder = output_folder,
      file.name = file_name,
      norm = TRUE,
      patient_idx = pt_idx
    )
  }
  
  cat("Processing interaction distribution and generating CSV output...\n")

  res = compute_results_processing(
    celltypes = cellTypes,
    patient_names = rownames(Cmatrix), 
    triangle_type = communication_type,
    sim_raw_file = paste0(output_folder, "/", file_name, ".out"),
    sim_norm_file = if(norm) paste0(output_folder, "/", file_name, "_norm.out") else NULL,
    remove_direction = remove_direction,
    normalized = norm,
    output_folder = output_folder,
    file.name = paste0(file_name, ".csv")
  )
  
  return(list(input = list(LRmatrix = LRmatrix, Lmatrix = Lmatrix, Rmatrix = Rmatrix, Cmatrix = Cmatrix, cellTypes = cellTypes, ligs = ligs, recs = recs), output = res))

}
