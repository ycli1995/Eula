
#' Convert to dense matrix
#'
#' Convert a matrix-like object into dense matrix in row chunks.
#'
#' @param object An object
#' @param ... `r .dot_param`
#'
#' @returns
#' A dense `matrix` converted from `object`.
#'
#' @name asMatrix
#' @export asMatrix
asMatrix <- function(object, ...) {
  UseMethod("asMatrix", object)
}

#' @param chunk.size An integer indicating number of rows to convert in each
#' chunk.
#' @param verbose `r .vb_param`
#'
#' @importFrom Eula.utils verboseMsg
#' 
#' @rdname asMatrix
#' @export
#' @method asMatrix default
asMatrix.default <- function(object, chunk.size = 1000, verbose = TRUE, ...) {
  breaks <- 0:ceiling(nrow(object) / chunk.size) * chunk.size
  k <- cut(seq_len(nrow(object)), breaks = breaks)
  mat <- matrix(0, nrow(object), ncol(object), dimnames = dimnames(object))
  for (i in levels(k)) {
    verboseMsg("as.matrix : ", paste0(range(which(k == i)), collapse = "-"))
    mat[which(k == i), ] <- as.matrix(object[which(k == i), ])
  }
  mat
}

#' @rdname asMatrix
#' @export
#' @method asMatrix matrix
asMatrix.matrix <- function(object, ...) {
  object
}

#' @rdname asMatrix
#' @export
#' @method asMatrix data.frame
asMatrix.data.frame <- function(object, ...) {
  as.matrix(object)
}
