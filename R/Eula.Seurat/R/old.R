
#' @export
CheckMySeuratObj <- function(object, ...) {
  Message('>>>>> Check Seurat extra information...')
  object <- CheckSeuratMetaData(object, ...)
  object <- CheckSeuratFData(object)
  object
}

#' @export
CheckSeuratMetaData <- function(object, column.map = list(), ...) {
  default.cols <- list(
    Samples = "orig.ident",
    Groups = c("group", "Samples"),
    Clusters = "seurat_clusters"
  )

  old.colors <- list(
    seurat_clusters = "color.cluster",
    orig.ident = "color.sample",
    Groups = "color.group"
  )
  for (i in names(old.colors)) {
    if (length(object@misc$colors[[i]]) == 0) {
      object@misc$colors[[i]] <- object@misc[[old.colors[[i]]]]
    }
  }

  check.cols <- names(column.map)
  for (i in check.cols) {
    if (!i %in% colnames(object[[]])) {
      use.col <- column.map[[i]][1]
      if (length(use.col) == 0) {
        use.col <- intersect(default.cols[[i]], colnames(object[[]]))[1]
      }
      if (length(use.col) == 0) {
        fastWarning("No column found for setting default '", i, "'")
        next
      }
      if (!use.col %in% colnames(object[[]])) {
        fastWarning("'", use.col, "' doesn't exist for default '", i, "'.")
        next
      }
      object@meta.data[, i] <- as.factor(object@meta.data[, use.col])

      if (length(object@misc$colors[[i]]) > 0) {
        object@misc$colors[[i]] <- checkColorMap(
          x = object[[]],
          column = i,
          colors = object@misc$colors[[i]]
        )
      } else {
        object@misc$colors[[i]] <- checkColorMap(
          x = object[[]],
          column = use.col,
          colors = object@misc$colors[[use.col]]
        )
      }
    }
    object@misc$colors[[i]] <- checkColorMap(
      x = object[[]],
      coloumn = i,
      colors = object@misc$colors[[i]]
    )
  }
  return(object)
}

ROW_DATA_NAMES <- c("id", "name", "unique_name", "merge_name")

#' @export
CheckSeuratFData <- function(object) {
  fdata <- AddUnderscore(object@misc$fdata)
  if (is.null(fdata)) {
    return(object)
  }
  if (is.data.frame(object@misc$rowData)) {
    checkColumns(object@misc$rowData, ROW_DATA_NAMES)
    return(object)
  }
  object@misc$rowData <- data.frame(
    id = rownames(fdata),
    name = fdata$name,
    unique_name = make.unique(fdata$name),
    merge_name = fdata$merge_name,
    row.names = fdata$dash
  )
  object
}

FindFeaturesID <- function(object, features, unlist = TRUE) {
  object@misc$fdata <- AddUnderscore(object@misc$fdata)
  if ( all(rownames(object) %in% object@misc$fdata$dash) ) {
    rownames(object@misc$fdata) <- object@misc$fdata$dash
  }
  if (exists("fdata", object@misc)) {
    f1 <- toupper(rownames(object@misc$fdata))
    f2 <- toupper(object@misc$fdata$name)
    f3 <- toupper(object@misc$fdata$merge_name)
    f4 <- toupper(object@misc$fdata$underscore)
  }
  features <- sapply(features, function(x) {
    if ( !exists("fdata", object@misc) ) return(NULL)
    g1 <- f1 %in% toupper(x)
    if ( sum(g1) > 0 ) return(rownames(object@misc$fdata)[g1])
    g2 <- f2 %in% toupper(x)
    if ( sum(g2) > 0 ) return(rownames(object@misc$fdata)[g2])
    g3 <- f3 %in% toupper(x)
    if ( sum(g3) > 0 ) return(rownames(object@misc$fdata)[g3])
    g4 <- f4 %in% toupper(x)
    if ( sum(g4) > 0 ) return(rownames(object@misc$fdata)[g4])
    message("[WARNING] '", x, "' not found gene id.")
    return(NULL)
  }, simplify = FALSE)
  if ( unlist ) features <- unlist(features)
  return(features)
}

FindFeaturesName <- function(object, features, col = "merge_name", is.fast = FALSE) {
  if ( ! exists("fdata", object@misc) )
    return(features)
  object@misc$fdata <- AddUnderscore(object@misc$fdata)
  new <- gsub("_", "-", features)
  if ( all(new %in% object@misc$fdata$dash) ) {
    rownames(object@misc$fdata) <- object@misc$fdata$dash
    features <- new
  }
  if ( is.fast ) {
    Name <- object@misc$fdata[features, col]
    names(Name) <- features
    return(Name)
  }
  Name <- sapply(features, function(x) {
    id <- object@misc$fdata[x, col]
    ifelse(is.null(id)||is.na(id), x, id)
  })
  return(Name)
}

AddUnderscore <- function(data){
  if (!is.null(data)) {
    if ( ! exists("underscore", data) || ! exists("dash", data) ) {
      data$underscore <- rownames(data)
      data$dash <- gsub("_", "-", rownames(data))
    }
  }
  return(data)
}
