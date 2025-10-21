
#' @export getFeaturesID
getFeaturesID <- function(object, features, ...) {
  UseMethod("getFeaturesID", object)
}

#' @export
#' @method getFeaturesID data.frame
getFeaturesID.data.frame <- function(object, features, uniq = FALSE, ...) {
  cols <- c("id", "name", "unique_name", "merge_name")
  checkColumns(object, cols = cols)
  for (i in cols) {
    object[[i]] <- toupper(object[[i]])
  }

  features <- sapply(features, function(x) {
    x <- toupper(x)
    g1 <- toupper(object$id) %in% x
    if (any(g1)) {
      return(rownames(object)[g1])
    }
    g2 <- toupper(object$name) %in% x
    if (any(g2)) {
      return(rownames(object)[g2])
    }
    g3 <- toupper(object$unique_name) %in% x
    if (any(g3)) {
      return(rownames(object)[g3])
    }
    g4 <- toupper(object$merge_name) %in% x
    if (any(g4)) {
      return(rownames(object)[g4])
    }
    g5 <- toupper(rownames(object)) %in% x
    if (any(g5)) {
      return(rownames(object)[g5])
    }
    message("[WARNING] '", x, "' not found gene id.")
    return(NULL)
  }, simplify = FALSE)
  features <- unlist(features)
  if (uniq) {
    features <- features[!duplicated(names(features))]
  }
  return(features)
}

#' @export getFeaturesName
getFeaturesName <- function(object, features, ...) {
  UseMethod("getFeaturesName", object)
}

#' @export
#' @method getFeaturesName data.frame
getFeaturesName.data.frame <- function(
    object,
    features,
    col = "unique_name",
    uniq = FALSE,
    ...
) {
  cols <- c("id", "name", "unique_name", "merge_name")
  checkColumns(object, cols = cols)
  col <- match.arg(col, choices = cols)

  features.id <- getFeaturesID(object, features = features)
  out <- object[features.id, col]
  names(out) <- features.id

  if (uniq) {
    out <- setNames(make.unique(out), names(out))
  }
  out
}
