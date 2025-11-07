
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
    g5 <- rownames(object) %in% x
    if (any(g5)) {
      return(rownames(object)[g5])
    }
    x2 <- toupper(x)
    g1 <- object$id %in% x2
    if (any(g1)) {
      return(unique(rownames(object)[g1]))
    }
    g2 <- object$name %in% x2
    if (any(g2)) {
      return(unique(rownames(object)[g2]))
    }
    g3 <- object$unique_name %in% x2
    if (any(g3)) {
      return(unique(rownames(object)[g3]))
    }
    g4 <- object$merge_name %in% x2
    if (any(g4)) {
      return(unique(rownames(object)[g4]))
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
