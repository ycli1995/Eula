
#' Compute Fold Change Between Two Cell Groups
#'
#' Computes per-feature fold changes between two sets of columns (e.g. two cell
#' groups) in a matrix-like object.
#'
#' @param object A numeric matrix-like object where rows represent features
#' (e.g. genes) and columns represent cells.
#' @param ... `r .dot_param`
#'
#' @export foldChange
foldChange <- function(object, ...) {
  UseMethod("foldChange", object)
}

#' @param cells.1 A vector of column indices or names defining group 1.
#' @param cells.2 A vector of column indices or names defining group 2.
#' @param mean.fxn A function used to calculate the means of features. Default
#' is [rowExpMean()] in non-logged space.
#' @param min.exp Minimum threshold used when computing the percentage of
#' expressing cells. Default is 0.
#' @param pseudocount.use Pseudocount added to group means before computing
#' fold change. Default is 1.
#' @param base Logarithm base used for fold change calculation. Default is 2.
#' - If `base > 0`, log fold change is computed:
#'   ```r
#'   log(mean.1) - log(mean.2)
#'   ```
#' - If `base < 0`, simple mean difference is returned:
#'   ```r
#'   mean.1 - mean.2
#'   ```
#' @param verbose `r .vb_param`
#'
#' @return
#' A `data.frame` with one row per feature and columns:
#'
#' - `avg_log{base}FC` (or `avg_diff` if `base < 0`)
#' - `mean.1`, `mean.2`: average expression per group
#' - `ncells.1`, `ncells.2`: number of cells in each group
#' - `pct.1`, `pct.2`: percentage of expressing cells in each group
#'
#' @importFrom Eula.utils verboseMsg
#' @importFrom MatrixGenerics rowSums
#'
#' @rdname foldChange
#' @export
#' @method foldChange CsparseMatrix
foldChange.CsparseMatrix <- function(
    object,
    cells.1,
    cells.2,
    mean.fxn = rowExpMean,
    min.exp = 0,
    pseudocount.use = 1,
    base = 2,
    verbose = TRUE,
    ...
) {
  verboseMsg("Calculating fold changes for ", nrow(object), " features.")
  out <- rowMeanPct(
    object = object,
    cell.groups = list(cells.1 = cells.1, cells.2 = cells.2),
    mean.fxn = mean.fxn,
    min.exp = min.exp
  )
  colnames(out$n.cells) <- c("ncells.1", "ncells.2")
  colnames(out$avg.pct) <- c("pct.1", "pct.2")
  colnames(out$avg.exp) <- c("mean.1", "mean.2")
  p1 <- pseudocount.use / length(cells.1)
  p2 <- pseudocount.use / length(cells.2)
  if (base < 0) {
    fc <- (out$avg.exp[, 1] + p1) - (out$avg.exp[, 2] + p2)
    fc.name <- "avg_diff"
  } else {
    fc <- log(out$avg.exp[, 1] + p1, base) - log(out$avg.exp[, 2] + p2, base)
    if (base == exp(1)) {
      base <- ""
    }
    fc.name <- paste0("avg_log", base, "FC")
  }
  fc.results <- cbind(fc, out$avg.exp, out$n.cells, out$avg.pct)
  colnames(fc.results)[1] <- fc.name
  fc.results <- as.data.frame(fc.results)
  fc.results
}

#' @rdname foldChange
#' @export
#' @method foldChange matrix
foldChange.matrix <- function(
    object,
    cells.1,
    cells.2,
    mean.fxn,
    min.exp = 0,
    pseudocount.use = 1,
    base = 2,
    verbose = TRUE,
    ...
) {
  foldChange.CsparseMatrix(
    object = object,
    cells.1 = cells.1,
    cells.2 = cells.2,
    mean.fxn = mean.fxn,
    min.exp = min.exp,
    pseudocount.use = pseudocount.use,
    base = base,
    verbose = verbose,
    ...
  )
}

#' Select genes from differential expression analysis
#'
#' Filter fold-change results (typically produced by [foldChange()]) according
#' to multiple expression and detection thresholds.
#'
#' @param fc.results A [data.frame()] from [foldChange()] containing:
#'   - a first column representing fold change (`avg_log2FC`, `avg_diff`, etc.)
#'   - `mean.1`, `mean.2`
#'   - `ncells.1`, `ncells.2`
#'   - `pct.1`, `pct.2`
#' @param min.diff Minimum absolute difference or log fold change) required.
#' Must be `>= 0`. Default is 0.1.
#' @param min.mean.exp Minimum mean expression required in at least one group.
#' @param min.pct Select genes that are detected in a minimum fraction of
#' `min.pct` cells in either of the two populations.
#' @param min.diff.pct Minimum difference in detection fraction between groups.
#' @param min.cells.feature Minimum number of cells expressing the feature
#' in at least one group.
#' @param only.pos Logical. If `TRUE`, only retain features with positive
#' difference.
#' @param verbose `r .vb_param`
#'
#' @return
#' A filtered [data.frame()] containing only features that pass all thresholds.
#'
#' @seealso [foldChange()]
#'
#' @importFrom Eula.utils verboseMsg fastWarn
#' @export
selectDE <- function(
    fc.results,
    min.diff = 0.1,
    min.mean.exp = 0,
    min.pct = 0.01,
    min.diff.pct = -Inf,
    min.cells.feature = 3,
    only.pos = FALSE,
    verbose = TRUE
) {
  use.cols <- as.vector(vapply(
    X = c("mean.", "ncells.", "pct."),
    FUN = paste0,
    FUN.VALUE = character(2L),
    1:2
  ))
  miss.cols <- setdiff(use.cols, colnames(fc.results)[-1])
  if (length(miss.cols) > 0) {
    stop(
      "Missing columns in 'fc.results':\n ",
      paste(miss.cols, collapse = ", ")
    )
  }
  if (min.diff < 0) {
    stop("'min.diff' must be >= 0.")
  }
  if (only.pos) {
    selected <- fc.results[, 1] >= min.diff
  } else {
    selected <- abs(fc.results[, 1]) >= min.diff
  }
  pct.max <- pmax(fc.results$pct.1, fc.results$pct.2)
  pct.diff <- pct.max - pmin(fc.results$pct.1, fc.results$pct.2)
  selected <- selected &
    (pmax(fc.results$ncells.1, fc.results$ncells.2) >= min.cells.feature) &
    (pmax(fc.results$mean.1, fc.results$mean.2) >= min.mean.exp) &
    (pct.max >= min.pct) &
    (pct.diff >= min.diff.pct)
  features <- rownames(fc.results)[selected]
  if (length(features) == 0) {
    fastWarn(
      "No feature pass thresholds:\n ",
      "min.cells.feature = ", min.cells.feature, "\n ",
      "min.pct = ", min.pct, "\n ",
      "min.diff.pct = ", min.diff.pct, "\n ",
      "min.diff = ", min.diff, "\n ",
      "only.pos = ", only.pos
    )
  }
  verboseMsg("Select ", length(features), " features")
  fc.results <- fc.results[features, , drop = FALSE]
  fc.results
}
