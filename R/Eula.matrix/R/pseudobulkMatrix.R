#' Construct a Pseudobulk Expression Matrix
#'
#' An S3 generic function to generates a pseudobulk matrix by aggregating or
#' averaging columns (e.g. single cells) of a matrix-like object according to
#' grouping labels.
#'
#' @param object A numeric matrix-like object where rows typically represent
#' features (e.g. genes) and columns represent cells or samples.
#' @param ... `r .dot_param`
#'
#' @export pseudobulkMatrix
pseudobulkMatrix <- function(object, ...) {
  UseMethod("pseudobulkMatrix", object)
}

#' @param cell.groups A vector of grouping labels with length equal to
#' `ncol(object)`. Each element specifies the group assignment of a column.
#' @param method The method used for calculating pseudobulk expression; one of:
#' "average" or "aggregate"
#' @param reverse Logical. If `TRUE`, each column of the output matrix will
#' represent summaries without corresponding category.
#' @param bulk Logical. If `TRUE`, an additional `"bulk"` column is added
#' containing sums (`aggregate`) or means (`average`) across all columns.
#' @param sparse Logical. Whether to return a `CsparseMatrix` or `matrix`.
#' @param drop Logical. Whether to drop empty groups (groups with zero columns).
#'

#' @return
#' A dense or sparse matrix with:
#' - Rows corresponding to `rownames(object)`
#' - Columns corresponding to the levels of `cell.groups`
#' - If `bulk = TRUE`, the first column is `"bulk"`.
#'
#' @examples
#' set.seed(1)
#' mat <- matrix(rpois(20, lambda = 5), nrow = 5)
#' groups <- rep(c("A", "B"), each = 2)
#'
#' # Aggregate (sum) by group
#' pseudobulkMatrix(mat, groups, method = "aggregate")
#'
#' # Average by group
#' pseudobulkMatrix(mat, groups, method = "average")
#'
#' # Include bulk column
#' pseudobulkMatrix(mat, groups, bulk = TRUE)
#'
#' @seealso [categoryMatrix()]
#'
#' @importFrom MatrixGenerics rowMeans rowSums
#' @importFrom methods as
#' @importClassesFrom Matrix CsparseMatrix
#'
#' @rdname pseudobulkMatrix
#' @export
#' @method pseudobulkMatrix CsparseMatrix
pseudobulkMatrix.CsparseMatrix <- function(
  object,
  cell.groups,
  method = c("aggregate", "average"),
  reverse = FALSE,
  bulk = FALSE,
  sparse = FALSE,
  drop = TRUE,
  ...
) {
  method <- match.arg(method)
  if (length(cell.groups) != ncol(object)) {
    stop("'cell.groups' must have the same length as 'ncol(object)'")
  }
  cell.groups <- as.factor(cell.groups)
  cat.mat <- categoryMatrix(
    labels = cell.groups,
    method = method,
    reverse = reverse,
    drop = drop
  )
  mat <- object %*% cat.mat
  rownames(mat) <- rownames(object)
  colnames(mat) <- levels(cell.groups)
  if (bulk) {
    func <- switch(method,
      "aggregate" = rowSums,
      "average" = rowMeans
    )
    bulk.mat <- func(object)
    mat <- cbind("bulk" = bulk.mat, mat)
  }
  if (sparse) {
    if (inherits(mat, "CsparseMatrix")) {
      return(mat)
    }
    return(as(mat, "CsparseMatrix"))
  }
  mat <- as.matrix(mat)
  mat
}

#' @export
#' @method pseudobulkMatrix matrix
pseudobulkMatrix.matrix <- function(
  object,
  cell.groups,
  method = c("aggregate", "average"),
  reverse = FALSE,
  bulk = FALSE,
  sparse = FALSE,
  drop = TRUE,
  ...
) {
  pseudobulkMatrix.CsparseMatrix(
    object = object,
    cell.groups = cell.groups,
    method = method,
    reverse = reverse,
    bulk = bulk,
    sparse = sparse,
    drop = drop,
    ...
  )
}

#' Construct a Category (Grouping) Design Matrix
#'
#' Creates a design matrix representing group membership of samples (or cells).
#' This matrix is typically used for pseudobulk aggregation:
#'
#' @param labels A vector of grouping labels.
#' @param method Character string specifying scaling:
#'   - `"aggregate"`: return a binary indicator matrix (0/1)
#'   - `"average"`: scale each column by `1 / group_size`
#' @param reverse Logical. If `TRUE`, invert the indicator matrix (1 becomes 0
#' and 0 becomes 1). When `labels` only contain 1 level, `reverse = TRUE` is
#' invalid and will raise an error.
#' @param drop Logical. If `TRUE`, remove groups with zero samples.
#'
#' @return
#' A numeric matrix (typically sparse) with:
#' - Rows corresponding to elements of `labels`
#' - Columns corresponding to unique label levels
#'
#' @examples
#' labels <- c("A", "A", "B", "B", "B")
#'
#' # Indicator matrix
#' categoryMatrix(labels, method = "aggregate")
#'
#' # Averaging matrix
#' categoryMatrix(labels, method = "average")
#'
#' @seealso [pseudobulkMatrix()]
#'
#' @importFrom MatrixGenerics colMeans colSums
#' @importFrom Matrix Diagonal sparse.model.matrix
#' @importFrom Eula.utils fastWarn
#' @export
categoryMatrix <- function(
  labels,
  method = c("aggregate", "average"),
  reverse = FALSE,
  drop = TRUE
) {
  method <- match.arg(method)
  if (length(unique(labels)) == 1) {
    fastWarn("All grouping variables have 1 value only.")
    if (reverse) {
      stop("Cannot use `reverse = TRUE`.")
    }
    cat.mat <- matrix(1, length(labels), 1)
    colnames(cat.mat) <- unique(as.character(labels))
    if (method == "aggregate") {
      return(cat.mat)
    }
    return(cat.mat / sum(cat.mat))
  }
  labels <- as.factor(labels)
  cat.mat <- sparse.model.matrix(~ 0 + labels)
  colnames(cat.mat) <- levels(labels)
  if (reverse) {
    cat.mat <- 1 - cat.mat
  }
  colsums <- colSums(cat.mat)
  if (drop) {
    cat.mat <- cat.mat[, colsums > 0, drop = FALSE]
    colsums <- colsums[colsums > 0]
  }
  if (method == "aggregate") {
    return(cat.mat)
  }
  cat.mat <- cat.mat %*% Diagonal(x = 1 / colsums, names = colnames(cat.mat))
  cat.mat
}
