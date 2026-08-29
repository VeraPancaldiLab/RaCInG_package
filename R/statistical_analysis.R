#' Run Wilcoxon tests across network features
#'
#' Runs a Wilcoxon rank-sum test per feature (column) between every pair of
#' groups. Also reports the fold change (ratio of group means) so results can
#' be plotted as a standard volcano plot with [volcano_plot()].
#'
#' @param data_matrix Numeric matrix or data frame with patients in rows and features in columns.
#' @param groups Vector of group labels with length matching `nrow(data_matrix)`.
#'   When more than two groups are present, every pairwise comparison is run.
#' @param p_adjust_method Multiple-testing correction method passed to [stats::p.adjust()].
#'
#' @return A data frame with one row per feature per group comparison, including
#'   `Comparison`, `Wilcox_statistic`, `Fold_change`, `Log2FC`, `P_value`, and
#'   `Adjusted_P_value` (adjusted jointly across all comparisons).
#' @export
wilcox_group_test <- function(data_matrix, groups, p_adjust_method = "fdr") {
  if (length(groups) != nrow(data_matrix)) {
    stop("Length of groups must match the number of rows in data_matrix")
  }

  groups <- as.factor(groups)
  levs <- levels(groups)
  if (length(levs) < 2) {
    stop("groups must contain at least two distinct values")
  }
  pairs <- utils::combn(levs, 2, simplify = FALSE)
  features <- colnames(data_matrix)

  comparisons <- lapply(pairs, function(pair) {
    g1 <- pair[1]
    g2 <- pair[2]
    t1_idx <- groups == g1
    t2_idx <- groups == g2

    p_values <- numeric(length(features))
    stat_values <- numeric(length(features))
    fold_change <- numeric(length(features))

    for (i in seq_along(features)) {
      t1 <- data_matrix[t1_idx, i]
      t2 <- data_matrix[t2_idx, i]

      test_res <- stats::wilcox.test(t1, t2)
      p_values[i] <- test_res$p.value
      stat_values[i] <- test_res$statistic
      fold_change[i] <- mean(t1) / mean(t2)
    }

    data.frame(
      Comparison = paste0(g1, "_vs_", g2),
      Feature = features,
      Wilcox_statistic = stat_values,
      Fold_change = fold_change,
      Log2FC = log2(fold_change),
      P_value = p_values,
      stringsAsFactors = FALSE
    )
  })

  result_df <- do.call(rbind, comparisons)
  result_df$Adjusted_P_value <- stats::p.adjust(result_df$P_value, method = p_adjust_method)
  result_df <- result_df[order(result_df$Adjusted_P_value), ]
  rownames(result_df) <- NULL

  return(result_df)
}

#' Create a volcano plot from Wilcoxon results
#'
#' @param wilcox_results Output of [wilcox_group_test()]. If it contains more
#'   than one `Comparison`, only the first is plotted; subset beforehand to
#'   plot a specific comparison.
#' @param top_labels Number of top significant features to label.
#' @param p_threshold Adjusted p-value threshold used to mark significance.
#' @param title Plot title.
#'
#' @return A `ggplot2` object.
#' @export
volcano_plot <- function(wilcox_results, top_labels = 10,
                                     p_threshold = 0.05, title = "Volcano Plot") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package `ggplot2` is required for `volcano_plot()`. Please install it first.", call. = FALSE)
  }

  if (!all(c("Feature", "Adjusted_P_value", "Log2FC") %in% colnames(wilcox_results))) {
    stop("wilcox_results must contain columns: Feature, Log2FC, Adjusted_P_value (the output of wilcox_group_test())")
  }

  if ("Comparison" %in% colnames(wilcox_results) && length(unique(wilcox_results$Comparison)) > 1) {
    first_comparison <- wilcox_results$Comparison[1]
    wilcox_results <- wilcox_results[wilcox_results$Comparison == first_comparison, ]
    title <- paste0(title, " (", first_comparison, ")")
  }

  # Build plot dataframe, dropping non-finite fold changes (e.g. a feature
  # that is all-zero in one of the two groups)
  plot_df <- data.frame(
    Feature = wilcox_results$Feature,
    Effect = wilcox_results$Log2FC,
    Adj_P_value = wilcox_results$Adjusted_P_value
  )
  plot_df <- plot_df[is.finite(plot_df$Effect), ]

  # Significant features
  plot_df$Significant <- ifelse(plot_df$Adj_P_value < p_threshold, "Yes", "No")

  # Identify top features by significance
  top_features <- utils::head(plot_df[plot_df$Significant == "Yes", ][order(plot_df$Adj_P_value), "Feature"], top_labels)

  # Plot
  ggplot2::ggplot(plot_df, ggplot2::aes(x = Effect, y = -log10(Adj_P_value))) +
    ggplot2::geom_point(ggplot2::aes(color = Significant), shape = 4, alpha = 0.7) +
    ggplot2::geom_text(
      data = subset(plot_df, Feature %in% top_features),
      ggplot2::aes(label = Feature),
      size = 3,
      vjust = -0.5,
      check_overlap = TRUE
    ) +
    ggplot2::geom_hline(yintercept = -log10(p_threshold), color = "red", linetype = "dashed") +
    ggplot2::scale_color_manual(values = c("No" = "gainsboro", "Yes" = "blue")) +
    ggplot2::xlim(min(plot_df$Effect) - 1, max(plot_df$Effect) + 1) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = "log2(fold change)",
      y = "-log10(Adjusted P-value)",
      title = title,
      color = "Significant"
    ) +
    ggplot2::theme(legend.position = "top")
}

#' Correlate network features with an external score
#'
#' Computes the Spearman correlation between an external per-patient score
#' (e.g. an immune response score) and every feature column, pooled across
#' all patients ("All") and, optionally, within each group separately.
#'
#' @param feature_matrix Numeric matrix or data frame with patients in rows
#'   (identified by `rownames()`) and features in columns.
#' @param score Named numeric vector of per-patient scores. Names must match
#'   `rownames(feature_matrix)`; patients missing a score are dropped.
#' @param group Optional grouping vector (e.g. cancer type or cohort), in the
#'   same row order as `feature_matrix`. When supplied, correlations are also
#'   computed separately within each group.
#' @param p_adjust_method Multiple-testing correction method passed to [stats::p.adjust()].
#'
#' @return A data frame with `Group`, `Feature`, `Rho`, `P_value`, and `Adjusted_P_value`.
#' @export
correlate_features_with_score <- function(feature_matrix, score, group = NULL, p_adjust_method = "fdr") {
  original_rownames <- rownames(feature_matrix)
  patients <- intersect(original_rownames, names(score))
  if (length(patients) == 0) {
    stop("No patients in common between rownames(feature_matrix) and names(score)")
  }

  if (!is.null(group)) group <- group[match(patients, original_rownames)]
  feature_matrix <- feature_matrix[patients, , drop = FALSE]
  score <- score[patients]

  correlate_one_group <- function(mat, sc, label) {
    features <- colnames(mat)
    rho <- numeric(length(features))
    p_values <- numeric(length(features))
    for (i in seq_along(features)) {
      test_res <- stats::cor.test(mat[, i], sc, method = "spearman")
      rho[i] <- unname(test_res$estimate)
      p_values[i] <- test_res$p.value
    }
    data.frame(Group = label, Feature = features, Rho = rho, P_value = p_values, stringsAsFactors = FALSE)
  }

  result_df <- correlate_one_group(feature_matrix, score, "All")

  if (!is.null(group)) {
    for (g in unique(group)) {
      result_df <- rbind(result_df, correlate_one_group(feature_matrix[group == g, , drop = FALSE], score[group == g], g))
    }
  }

  result_df$Adjusted_P_value <- stats::p.adjust(result_df$P_value, method = p_adjust_method)
  result_df <- result_df[order(result_df$Group, result_df$Adjusted_P_value), ]
  rownames(result_df) <- NULL

  return(result_df)
}

#' Rainfall plot of feature-score correlations
#'
#' Shows the top positive and top negative correlations from
#' [correlate_features_with_score()] as a diverging bar chart.
#'
#' @param correlation_results Output of [correlate_features_with_score()].
#' @param group Which `Group` value to plot. Defaults to the first group present.
#' @param top_n Number of top positive and top negative correlations to show.
#' @param title Plot title.
#'
#' @return A `ggplot2` object.
#' @export
correlation_plot <- function(correlation_results, group = NULL, top_n = 10, title = "Feature correlation") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package `ggplot2` is required for `correlation_plot()`. Please install it first.", call. = FALSE)
  }

  if (is.null(group)) group <- correlation_results$Group[1]
  df <- correlation_results[correlation_results$Group == group, ]

  top_pos <- utils::head(df[order(-df$Rho), ], top_n)
  top_neg <- utils::head(df[order(df$Rho), ], top_n)
  plot_df <- unique(rbind(top_neg, top_pos))
  plot_df <- plot_df[order(plot_df$Rho), ]
  plot_df$Feature <- factor(plot_df$Feature, levels = plot_df$Feature)
  plot_df$Direction <- ifelse(plot_df$Rho >= 0, "Positive", "Negative")

  ggplot2::ggplot(plot_df, ggplot2::aes(x = Feature, y = Rho, fill = Direction)) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(values = c("Positive" = "forestgreen", "Negative" = "firebrick")) +
    ggplot2::coord_flip() +
    ggplot2::ylim(-1, 1) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = NULL, y = "Spearman rho", title = paste0(title, " (", group, ")"), fill = NULL) +
    ggplot2::theme(legend.position = "top")
}
