#' SKCM melanoma example input for RaCInG
#'
#' Pre-processed input matrices derived from TCGA skin cutaneous melanoma (SKCM)
#' bulk RNA-seq data. The object is a named list that can be passed directly to
#' `compute_racing_kernel()` or `compute_racing_montecarlo()` via the
#' `input_data` parameter.
#'
#' @format A named list with the following elements:
#' \describe{
#'   \item{Lmatrix}{Numeric matrix (9 cell types x 276 ligands). Binary
#'     cell-to-ligand compatibility.}
#'   \item{Rmatrix}{Numeric matrix (9 cell types x 298 receptors). Binary
#'     cell-to-receptor compatibility.}
#'   \item{Cmatrix}{Numeric matrix (10 patients x 9 cell types). Row-normalised
#'     cell-type abundance estimates from deconvolution, with M1 and M2
#'     macrophages merged into a single M category. Columns are named as in
#'     `celltypes`; rows are unnamed, so the workflows label patients
#'     `Patient_1` ... `Patient_10`.}
#'   \item{LRmatrix}{3-D numeric array (276 ligands x 298 receptors x 10
#'     patients). Normalised ligand–receptor interaction probability tensor.}
#'   \item{celltypes}{Character vector of 9 cell-type names (alphabetically
#'     sorted).}
#'   \item{ligands}{Character vector of 276 ligand names.}
#'   \item{receptors}{Character vector of 298 receptor names.}
#'   \item{Sign_matrix}{Numeric matrix (276 x 298) of zeros (unknown
#'     interaction signs).}
#' }
#'
#' @source Derived from TCGA SKCM data processed with TMEmod deconvolution
#'   and OmniPath ligand–receptor annotations.
#'
#' @examples
#' data(skcm_example)
#' str(skcm_example, max.level = 1)
#'
#' # Use directly with the kernel workflow
#' # result <- compute_racing_kernel(input_data = skcm_example)
"skcm_example"
