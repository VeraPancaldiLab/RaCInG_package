#' Decompose a Direct communication score into per-ligand-receptor-pair contributions
#'
#' `compute_kernel()` aggregates every ligand-receptor pair into
#' `kernel[sender, receiver] = lig_weight[sender, ] %*% LRmatrix %*% rec_weight[receiver, ]`
#' -- a sum over all (ligand, receptor) pairs -- so once the kernel is built,
#' individual ligand-receptor identity is gone from the number alone. Because
#' that formula is linear in `LRmatrix`, the contribution of one specific
#' ligand-receptor pair `(lig, rec)` to `kernel[sender, receiver]` is just the
#' single term `lig_weight[sender, lig] * LRmatrix[lig, rec] * rec_weight[receiver, rec]`,
#' and these terms sum exactly back to the full kernel value. This function
#' computes every such term for one sender-receiver pair, so the (weighted)
#' importance of each ligand-receptor pair to that pair's communication score
#' can be ranked.
#'
#' @param Lmatrix,Rmatrix,Cmatrix,LRmatrix Same inputs as [compute_kernel()].
#' @param celltypes Character vector of cell-type names, in the row order of
#'   `Lmatrix`/`Rmatrix`/columns of `Cmatrix`.
#' @param sender,receiver Cell-type names (must be in `celltypes`) for the
#'   pair of interest.
#' @param Dcell Patient-by-cell-type abundance matrix used to scale raw
#'   contributions to the same units as [calculateDirect()]'s raw score.
#'   Defaults to `Cmatrix`.
#' @param bundle Logical; if `TRUE` (default), pools both directions
#'   (`sender -> receiver` and `receiver -> sender`) under one `Edge` label
#'   and one shared `PctOfTotal` denominator, matching
#'   [calculateDirect()]'s bundled score exactly (a self-pair, `sender ==
#'   receiver`, is doubled the same way `calculateDirect()` doubles it). If
#'   `FALSE`, only the `sender -> receiver` direction is returned.
#' @param patient_names Optional character vector of patient labels, matching
#'   the patient order in `LRmatrix`. Defaults to `Patient_1`, `Patient_2`, ...
#' @param top_n Optional; if given, keep only the top `top_n` ligand-receptor
#'   pairs per patient (by `Contribution`).
#'
#' @return A data frame with one row per (patient, ligand-receptor pair):
#'   `Patient`, `Edge` (a label shared by rows that pool into one
#'   `PctOfTotal` denominator), `Sender`, `Receiver`, `Ligand`, `Receptor`,
#'   `Contribution` (in the same units as [calculateDirect()]'s raw score),
#'   and `PctOfTotal` (share of `Contribution` within that patient and
#'   `Edge`). Contributions that are exactly zero (no possible pathway) are
#'   omitted. `NULL` (with a message) if no ligand-receptor pair has a
#'   nonzero contribution at all.
#' @export
computeLRContributions <- function(Lmatrix, Rmatrix, Cmatrix, LRmatrix, celltypes,
                                   sender, receiver, Dcell = Cmatrix, bundle = TRUE,
                                   patient_names = NULL, top_n = NULL) {
  s <- match(sender, celltypes)
  r <- match(receiver, celltypes)
  if (is.na(s) || is.na(r)) {
    stop("`sender`/`receiver` must both be one of `celltypes`.", call. = FALSE)
  }

  n_patients <- dim(LRmatrix)[3]
  pnames <- if (!is.null(patient_names)) patient_names else paste0("Patient_", seq_len(n_patients))
  ligs <- colnames(Lmatrix)
  recs <- colnames(Rmatrix)
  Dcell <- as.matrix(Dcell)

  edges <- if (bundle && s != r) list(c(s, r), c(r, s)) else list(c(s, r))
  self_pair_mult <- if (bundle && s == r) 2 else 1
  edge_label <- if (bundle) paste0(sender, " <-> ", receiver) else paste0(sender, " -> ", receiver)

  rows <- list()
  for (p in seq_len(n_patients)) {
    Cvec <- Cmatrix[p, ]
    lig_norm <- colSums(sweep(Lmatrix, 1, Cvec, "*")); lig_norm[lig_norm == 0] <- 1
    lig_weight <- sweep(Lmatrix, 2, lig_norm, "/")
    rec_norm <- colSums(sweep(Rmatrix, 1, Cvec, "*")); rec_norm[rec_norm == 0] <- 1
    rec_weight <- sweep(Rmatrix, 2, rec_norm, "/")
    LR_p <- LRmatrix[, , p]
    abund_scale <- Dcell[p, s] * Dcell[p, r] * self_pair_mult

    for (e in edges) {
      contrib <- outer(lig_weight[e[1], ], rec_weight[e[2], ]) * LR_p * abund_scale
      idx <- which(contrib != 0, arr.ind = TRUE)
      if (nrow(idx) == 0) next
      rows[[length(rows) + 1]] <- data.frame(
        Patient = pnames[p], Edge = edge_label,
        Sender = celltypes[e[1]], Receiver = celltypes[e[2]],
        Ligand = ligs[idx[, 1]], Receptor = recs[idx[, 2]],
        Contribution = contrib[idx]
      )
    }
  }

  if (length(rows) == 0) {
    message("computeLRContributions(): no nonzero ligand-receptor contribution found for ", edge_label, ".")
    return(invisible(NULL))
  }

  result <- do.call(rbind, rows)
  result$PctOfTotal <- 100 * result$Contribution / stats::ave(result$Contribution, result$Patient, result$Edge, FUN = sum)
  result <- result[order(result$Patient, -result$Contribution), ]
  if (!is.null(top_n)) {
    result <- do.call(rbind, by(result, result$Patient, function(d) utils::head(d, top_n), simplify = FALSE))
  }
  rownames(result) <- NULL
  result
}

#' Decompose a Wedge's two hub paths into per-edge ligand-receptor contributions
#'
#' A Wedge (see [calculateWedges()]) sums two directed 2-hop paths through hub
#' `j`: `i -> j -> k` and `k -> j -> i`. Since Wedges multiply their two edges
#' together rather than sum them, there is no single well-defined "share of
#' the Wedge score" attributable to one edge's ligand-receptor pair alone --
#' instead, this decomposes each of the 4 edges across both paths separately
#' via [computeLRContributions()] (`bundle = FALSE`, one direction per edge),
#' so `PctOfTotal` in the result is relative to each edge's own kernel value.
#'
#' @inheritParams computeLRContributions
#' @param i,k Cell-type names for the wedge's two leaves.
#' @param j Cell-type name for the wedge's hub.
#'
#' @return A data frame like [computeLRContributions()]'s (pooled across the
#'   4 edges), with an added `Path` column (`"i -> j -> k"` or `"k -> j ->
#'   i"`, using the actual cell-type names).
#' @export
computeWedgeLRContributions <- function(Lmatrix, Rmatrix, Cmatrix, LRmatrix, celltypes,
                                        i, j, k, Dcell = Cmatrix, patient_names = NULL, top_n = NULL) {
  path_ijk <- paste0(i, " -> ", j, " -> ", k)
  path_kji <- paste0(k, " -> ", j, " -> ", i)

  edge_calls <- list(
    list(path = path_ijk, sender = i, receiver = j),
    list(path = path_ijk, sender = j, receiver = k),
    list(path = path_kji, sender = k, receiver = j),
    list(path = path_kji, sender = j, receiver = i)
  )

  result <- do.call(rbind, lapply(edge_calls, function(ec) {
    df <- computeLRContributions(Lmatrix, Rmatrix, Cmatrix, LRmatrix, celltypes,
                                 sender = ec$sender, receiver = ec$receiver, Dcell = Dcell,
                                 bundle = FALSE, patient_names = patient_names, top_n = NULL)
    if (is.null(df)) return(NULL)
    df$Path <- ec$path
    df
  }))

  if (is.null(result)) {
    message("computeWedgeLRContributions(): no nonzero ligand-receptor contribution found for wedge ", i, "-", j, "-", k, ".")
    return(invisible(NULL))
  }

  if (!is.null(top_n)) {
    result <- result[order(result$Patient, result$Edge, -result$Contribution), ]
    result <- do.call(rbind, by(result, list(result$Patient, result$Edge), function(d) utils::head(d, top_n), simplify = FALSE))
  }
  rownames(result) <- NULL
  result
}

#' Decompose a Triangle's edges into per-edge ligand-receptor contributions
#'
#' A bundled Triangle (see [computeTriangles()]) combines all 6 directed
#' edges among 3 cell types (`i`, `j`, `k`) into 8 direction-combination
#' products. As with Wedges, edges are multiplied together rather than summed,
#' so there is no single well-defined "share of the Triangle score"
#' attributable to one edge's ligand-receptor pair alone -- instead, this
#' decomposes each of the 6 directed edges separately via
#' [computeLRContributions()] (`bundle = FALSE`), so `PctOfTotal` in the
#' result is relative to each edge's own kernel value.
#'
#' @inheritParams computeLRContributions
#' @param i,j,k Cell-type names for the triangle's 3 vertices.
#'
#' @return A data frame like [computeLRContributions()]'s, pooled across the
#'   6 directed edges (`Edge` identifies which one, e.g. `"i -> j"`).
#' @export
computeTriangleLRContributions <- function(Lmatrix, Rmatrix, Cmatrix, LRmatrix, celltypes,
                                           i, j, k, Dcell = Cmatrix, patient_names = NULL, top_n = NULL) {
  edge_pairs <- list(c(i, j), c(j, i), c(j, k), c(k, j), c(k, i), c(i, k))

  result <- do.call(rbind, lapply(edge_pairs, function(e) {
    computeLRContributions(Lmatrix, Rmatrix, Cmatrix, LRmatrix, celltypes,
                           sender = e[1], receiver = e[2], Dcell = Dcell,
                           bundle = FALSE, patient_names = patient_names, top_n = NULL)
  }))

  if (is.null(result)) {
    message("computeTriangleLRContributions(): no nonzero ligand-receptor contribution found for triangle ", i, "-", j, "-", k, ".")
    return(invisible(NULL))
  }

  if (!is.null(top_n)) {
    result <- result[order(result$Patient, result$Edge, -result$Contribution), ]
    result <- do.call(rbind, by(result, list(result$Patient, result$Edge), function(d) utils::head(d, top_n), simplify = FALSE))
  }
  rownames(result) <- NULL
  result
}

#' Dot plot of top ligand-receptor pair contributions (LIANA-style)
#'
#' Ranks ligand-receptor pairs by mean contribution (within each `Edge`
#' separately, since different edges of a Wedge/Triangle are typically driven
#' by different pairs) and plots the top ones as a dot per patient x pair:
#' dot size is `PctOfTotal`, dot color is the raw `Contribution`.
#'
#' @param lr_contrib Output of [computeLRContributions()],
#'   [computeWedgeLRContributions()], or [computeTriangleLRContributions()].
#' @param top_n Number of top ligand-receptor pairs to show per `Edge`.
#' @param title Plot title. Defaults to `NULL`, which auto-generates one
#'   naming the cell type(s) involved (e.g. `"Top ligand-receptor pairs: CLL -
#'   NKcells"`), derived from `lr_contrib`'s `Sender`/`Receiver` columns.
#'
#' @return A `ggplot2` object.
#' @export
plotLRContributions <- function(lr_contrib, top_n = 15, title = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package `ggplot2` is required for `plotLRContributions()`. Please install it first.", call. = FALSE)
  }
  required_cols <- c("Patient", "Edge", "Sender", "Receiver", "Ligand", "Receptor", "Contribution", "PctOfTotal")
  if (!all(required_cols %in% colnames(lr_contrib))) {
    stop("`lr_contrib` must be the output of computeLRContributions()/computeWedgeLRContributions()/computeTriangleLRContributions().", call. = FALSE)
  }

  if (is.null(title)) {
    edges <- unique(lr_contrib$Edge)
    pair_label <- if (length(edges) == 1) {
      gsub(" <-> | -> ", " - ", edges)
    } else {
      paste(sort(unique(c(lr_contrib$Sender, lr_contrib$Receiver))), collapse = " - ")
    }
    title <- paste0("Top ", top_n, " ligand-receptor pairs by mean contribution: ", pair_label)
  }

  lr_contrib$LRpair <- paste0(lr_contrib$Ligand, " -> ", lr_contrib$Receptor)

  top_pairs <- lr_contrib %>%
    dplyr::group_by(Edge, LRpair) %>%
    dplyr::summarise(mean_contrib = mean(Contribution), .groups = "drop") %>%
    dplyr::group_by(Edge) %>%
    dplyr::slice_max(mean_contrib, n = top_n, with_ties = FALSE) %>%
    dplyr::ungroup()

  plot_df <- lr_contrib %>% dplyr::inner_join(top_pairs, by = c("Edge", "LRpair"))
  plot_df$LRpair <- stats::reorder(plot_df$LRpair, plot_df$mean_contrib)

  # Many patients/samples on the x axis quickly overlap under a shallow tilt;
  # a vertical (90 deg) label, shrunk once there are enough categories to
  # crowd each other, keeps every label readable without changing plot width.
  n_patients <- length(unique(plot_df$Patient))
  x_text_size <- if (n_patients > 20) 6 else if (n_patients > 10) 8 else 10

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = Patient, y = LRpair, size = PctOfTotal, color = Contribution)) +
    ggplot2::geom_point() +
    ggplot2::scale_color_viridis_c(option = "magma") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = x_text_size)) +
    ggplot2::labs(title = title, x = NULL, y = NULL, size = "% of edge total", color = "Contribution")

  if (length(unique(plot_df$Edge)) > 1) {
    p <- p + ggplot2::facet_wrap(~Edge, scales = "free_y")
  }
  p
}
