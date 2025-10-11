
#source("utils0.R")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Helpers ######################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## Get default #################################################################

getDefault <- function(key, map) {
  if (is.null(key)) {
    stop("Cannot found default value for an empty key")
  }
  key <- key[1]
  if (!key %in% names(map)) {
    warning(
      "Cannot found default value for key '", key, "'",
      call. = FALSE, immediate. = TRUE
    )
  }
  return(map[[key]])
}

findRC <- function(n) {
  row <- ceiling(sqrt(n))
  col <- ceiling(n/row)
  c(row, col)
}

FindWHNum <- function(data, ncol = NULL, nrow = NULL, by.col = TRUE) {
  if (length(data) == 1) {
    stopifnot(is.numeric(data))
    data <- round(data)
  } else {
    data <- length(data)
  }
  if (!is.null(ncol) && !is.null(nrow)) {
    stopifnot(ncol * nrow >= data)
    return(c(nrow, ncol))
  }
  if (!is.null(ncol)) {
    nrow <- ceiling(data / ncol)
    return(c(nrow, ncol))
  }
  if (!is.null(nrow)) {
    ncol <- ceiling(data / nrow)
    return(c(nrow, ncol))
  }
  if (by.col) {
    ncol <- ceiling(sqrt(data))
    nrow <- ceiling(data / ncol)
    return(c(nrow, ncol))
  }
  nrow <- ceiling(sqrt(data))
  ncol <- ceiling(data / nrow)
  return(c(nrow, ncol))
}

IfNull <- function(var, default = NULL){
  if (is.null(var)) {
    return(default)
  }
  return(var)
}

## Load data ###################################################################
Load <- function(file) {
  file <- tools::file_path_as_absolute(file)
  Message("Loading: ", file)
  object <- readRDX(file)
  if (!"version" %in% slotNames(object)) {
    return(object)
  }
  if (grepl('^2', object@version)) {
    return(Seurat::UpdateSeuratObject(object))
  }
  return(object)
}

parse_func <- function(expr, name, ...) {
  str <- deparse(expr, ...)
  str[1] <- paste(name, "<-", str[1])
  str
}

## Rename factor ###############################################################

ChangeOUTName <- function(features, fdata) {
  features <- as.character(features)
  fdata <- AddUnderscore(fdata)
  if (!is.null(fdata) && all(features %in% fdata$dash)) {
    underscore_id <- fdata$underscore
    names(underscore_id) <- fdata$dash
    features <- underscore_id[features]
  }
  return(features)
}

AddUnderscore <- function(data){
  if (is.null(data)) {
    return(data)
  }
  if (!exists("underscore", data) || !exists("dash", data)) {
    data$underscore <- rownames(data)
    data$dash <- gsub("_", "-", rownames(data))
  }
  return(data)
}

dgC2dgR <- function(x) {
  return(sparseMatrix(
    j = x@i,
    p = x@p,
    x = x@x,
    repr = "R",
    index1 = FALSE,
    dims = rev(dim(x)),
    dimnames = rev(dimnames(x))
  ))
}

dgR2dgC <- function(x) {
  return(sparseMatrix(
    i = x@j,
    p = x@p,
    x = x@x,
    repr = "C",
    index1 = FALSE,
    dims = rev(dim(x)),
    dimnames = rev(dimnames(x))
  ))
}

AsMatrix <- function(object, max.row = 1000) {
  k <- cut(
    1:nrow(object),
    breaks = 0:ceiling(nrow(object) / max.row) * max.row
  )
  mat <- matrix(0, nrow(object), ncol(object), dimnames = dimnames(object))
  for (i in levels(k)) {
    message("as.matrix : ", paste0(range(which(k == i)), collapse = "-") )
    mat[which(k == i)] <- as.matrix(object[which(k == i), ])
  }
  return(mat)
}

# make percentage #
makePct <- function(x , digits = 3) {
  if (!anyNA(as.numeric(x))) {
    x <- paste0(round(as.numeric(x), digits = digits) *100, '%')
  }
  x
}

# make minMax #
minMax <- function(x) (x - min(x)) / (max(x) - min(x))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S4 methods ###################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## recodeFactor ################################################################

recodeFactor.default <- function(x, map, keep.orders = FALSE, ...) {
  map <- unlistMap(map = map)
  if (keep.orders) {
    recode_func <- dplyr::recode
  } else {
    recode_func <- dplyr::recode_factor
  }
  x %>%
    as.factor() %>%
    recode_func(!!!map) %>%
    as(Class = class(x))
}
registerS3method(
  genname = "recodeFactor",
  class = "default",
  method = recodeFactor.default,
  envir = environment(recodeFactor)
)

# setMethod(
#   f = "recodeFactor",
#   signature = c("vector", "list"),
#   definition = function(x, map, keep.orders = FALSE, ...) {
#     map <- unlistMap(map = map)
#     if (keep.orders) {
#       recode_func <- dplyr::recode
#     } else {
#       recode_func <- dplyr::recode_factor
#     }
#     x %>%
#       as.factor() %>%
#       recode_func(!!!map) %>%
#       as(Class = class(x))
#   }
# )
#
# setMethod(
#   f = "recodeFactor",
#   signature = c("data.frame", "list"),
#   definition = function(x, map, column, keep.orders = FALSE, ...) {
#     column <- column[1]
#     if (!column %in% colnames(x)) {
#       return(x)
#     }
#     message("Replace entries in column '", column, "':")
#     x[, paste0(column, "_old")] <- x[, column]
#     x[, column] <- recodeFactor(x[, column], map, keep.orders = keep.orders)
#     print(table(x[, column], x[, paste0(column, "_old")]))
#     return(x)
#   }
# )
#
# ## getFeaturesID ##############################################################
# setMethod(
#   f = "getFeaturesID",
#   signature = c("data.frame", "character"),
#   definition = function(object, features, unlist = TRUE, uniq = FALSE, ...) {
#     fdata <- AddUnderscore(object)
#
#     t1 <- toupper(rownames(fdata))
#     t2 <- toupper(fdata$name)
#     t3 <- toupper(fdata$merge_name)
#     t4 <- toupper(fdata$underscore)
#
#     features <- features %>%
#       sapply(function(x) {
#         g1 <- t1 %in% toupper(x)
#         if (sum(g1) > 0) {
#           return(rownames(fdata)[g1])
#         }
#         g2 <- t2 %in% toupper(x)
#         if (sum(g2) > 0) {
#           return(rownames(fdata)[g2])
#         }
#         g3 <- t3 %in% toupper(x)
#         if (sum(g3) > 0) {
#           return(rownames(fdata)[g3])
#         }
#         g4 <- t4 %in% toupper(x)
#         if (sum(g4) > 0) {
#           return(rownames(fdata)[g4])
#         }
#         message("[WARNING] '", x, "' not found gene id.")
#         return(NULL)
#       }, simplify = FALSE)
#     if (unlist) {
#       features <- unlist(features)
#     }
#     if (uniq) {
#       features <- features[!duplicated(names(features))]
#     }
#     return(features)
#   }
# )
#
# setMethod(
#   f = "getFeaturesName",
#   signature = c("data.frame", "character"),
#   definition = function(
#     object,
#     features,
#     col = "merge_name",
#     is.fast = FALSE,
#     uniq = FALSE,
#     ...
#   ) {
#     fdata <- AddUnderscore(object)
#     new <- gsub("_", "-", features)
#     if (all(new %in% fdata$dash) ) {
#       rownames(fdata) <- fdata$dash
#       features <- new
#     }
#     if (is.fast) {
#       Name <- fdata[features, col]
#       names(Name) <- features
#       return(Name)
#     }
#     out <- sapply(features, function(x) {
#       id <- fdata[x, col]
#       ifelse(is.null(id) || is.na(id), x, id)
#     }, simplify = FALSE)
#     out <- unlist(out)
#     if (uniq) {
#       out <- setNames(make.unique(out), names(out))
#     }
#     out
#   }
# )

