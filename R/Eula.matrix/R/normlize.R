#' @export relativeCounts
relativeCounts <- function(object, ...) {
  UseMethod("relativeCounts", object)
}

#' @export
#' @method relativeCounts matrix
relativeCounts.matrix <- function(object, scale.factor = 10000, ...) {
  scale.factor <- .rep_scale_factor(scale.factor, ncol(object))
  old.dimnames <- dimnames(object)
  object <- matrix_rc_norm(object, scale.factor)
  dimnames(object) <- old.dimnames
  object
}

#' @importClassesFrom Matrix dgCMatrix
#' @export
#' @method relativeCounts dgCMatrix
relativeCounts.dgCMatrix <- function(object, scale.factor = 10000, ...) {
  scale.factor <- .rep_scale_factor(scale.factor, ncol(object))
  object@x <- object@x /
    rep.int(Matrix::colSums(object) / scale.factor, diff(object@p))
  object
}

.rep_scale_factor <- function(scale.factor, ncol) {
  if (length(scale.factor) == 1) {
    scale.factor <- rep.int(scale.factor, ncol)
  }
  if (!length(scale.factor) == ncol) {
    stop(
      "scale.factor must have the same length as the number of columns",
      " in the matrix"
    )
  }
  scale.factor
}
