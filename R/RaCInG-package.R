#' RaCInG: Random Cell-cell Interaction Generator
#'
#' `RaCInG` provides functions to prepare ligand-receptor input matrices,
#' generate patient-specific communication graphs, derive kernel or Monte Carlo
#' features, and perform downstream statistical analysis.
#'
#' @keywords internal
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom foreach %dopar%
"_PACKAGE"

utils::globalVariables(
  c(
    "Adj_P_value", "Celltype", "Datasets", "Direction", "Effect",
    "Expr_Ligand", "Expr_Receptor", "Feature", "Ligand", "Receiver",
    "Receptor", "Rho", "Sender", "Significant", "category_intercell_source",
    "category_intercell_target"
  )
)

#' Check that optional packages are installed
#'
#' Internal helper used by high-level workflows that depend on optional
#' preprocessing packages.
#'
#' @param packages Character vector of package names.
#'
#' @return Invisibly returns `TRUE` when all packages are available.
#' @keywords internal
.check_installed_packages <- function(packages) {
  missing_pkgs <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]

  if (length(missing_pkgs) > 0) {
    stop(
      "The following optional packages are required for this step but are not installed: ",
      paste(missing_pkgs, collapse = ", "),
      ". Install them first or supply precomputed inputs.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

