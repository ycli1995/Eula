
#' Normalize and log-transform raw data
#'
#' Normalize each cell by relative counts, then natural log-transform the data.
#'
#' @param mat Matrix with the raw data.
#'
#' @export logNorm
logNorm <- function(object, ...) {
  UseMethod("logNorm", object)
}

#' @param scale.factor The total count after normalization. Default is 10000.
#'
#' @method logNorm matrix
#' @export
#' @rdname logNorm
logNorm.matrix <- function(object, scale.factor = 10000, ...) {
  scale.factor <- .rep_scale_factor(scale.factor, ncol(object))
  old.dimnames <- dimnames(object)
  object <- matrix_log_norm(object, scale.factor)
  dimnames(object) <- old.dimnames
  object
}

#' @importFrom Matrix colSums
#' @importClassesFrom Matrix dgCMatrix
#' @method logNorm dgCMatrix
#' @export
#' @rdname logNorm
logNorm.dgCMatrix <- function(object, scale.factor = 10000, ...) {
  scale.factor <- .rep_scale_factor(scale.factor, ncol(object))
  object@x <- log1p(
    object@x / rep.int(Matrix::colSums(object) / scale.factor, diff(object@p))
  )
  object
}

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

#' @importFrom Matrix colSums
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
