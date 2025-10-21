
#' @export
Load <- function(file) {
  object <- Eula.utils::Load(file)
  if (!is(object, "Seurat")) {
    return(object)
  }
  CheckMySeuratObj(object)
}

#' @importFrom SeuratObject FetchData
#' @importFrom Eula.utils validCharacters
#' @export
FetchSeuratData <- function(object, vars, ...) {
  new_names <- names(vars)
  df <- FetchData(object, vars = vars, ...)
  if (any(validCharacters(new_names))) {
    vars <- vars[validCharacters(new_names)]
    vars <- vars[vars %in% colnames(df)]
    df <- df %>%
      dplyr::rename(all_of(vars))
  }
  df <- df %>%
    dplyr::mutate(Cells = rownames(.), .before = 1) %>%
    as.data.frame()
  df
}

#' @importFrom dplyr select
#' @export
MergeFData <- function(object, assay = NULL) {
  fdata <- object@misc$rowData
  if (is.null(fdata)) {
    return(object)
  }
  assay <- assay %||% DefaultAssay(object)
  meta.features <- object[[assay]][[]]
  fdata <- fdata[rownames(meta.features), , drop = FALSE] %>%
    dplyr::select(setdiff(colnames(.), colnames(meta.features))) %>%
    cbind(meta.features)
  if (inherits(object[[assay]], "StdAssay")) {
    object[[assay]]@meta.data <- fdata
  } else {
    object[[assay]]@meta.features <- fdata
  }
  object
}

#' @export
RetrieveCounts <- function(mat, scale.factor = 10000) {
  nCounts <- apply(mat, 2, function(x) {
    exp1 <- sort(unique(x))
    exp1 <- exp1[exp1 != 0]

    ## First suppose the minimum non-zero is 1.
    N1 <- 1 / (expm1(exp1[1]) / scale.factor)

    ## Reverse it.
    c1 <- expm1(exp1) / scale.factor * N1

    ## Get differs between reversed values and the nearest integer.
    f1 <- round(abs(c1 - round(c1)), digits = 7)
    f1 <- f1[f1 != 0]

    x_min <- 1
    if (length(f1) > 0) {
      ## min(fl) == 1 / x_min
      x_min <- round(max(1 / f1))

      ## When f1 very closed to 0?
      x_min[is.infinite(x_min)] <- 1
    }
    return(round(x_min/ (expm1(exp1[1]) / scale.factor)))
  })
  mat@x <- round(expm1(mat@x) / scale.factor * rep.int(nCounts, diff(mat@p)))
  mat
}

#' @export
ShrinkSeuratObject <- function(
    object,
    assays = NULL,
    scale.data = NULL,
    misc.counts = NULL,
    ...
) {
  assays <- assays %||% SeuratObject::Assays(object)
  scale.data <- scale.data %||% FALSE
  misc.counts <- misc.counts %||% TRUE

  if (!all(SeuratObject::Assays(object) %in% assays)) {
    message("keep assays: ", paste(assays, collapse = ", "))
    if (!DefaultAssay(object) %in% assays) {
      DefaultAssay(object) <- assays[1]
    }
    object@assays <- object@assays[assays]
  }

  if (!scale.data) {
    message("remove 'scale.data' ")
    if (packageVersion("Seurat") >= 5) {
      SeuratObject::LayerData(object, layer = "scale.data") <- NULL
    } else {
      object <- SeuratObject::SetAssayData(
        object,
        slot = "scale.data",
        new.data = matrix(0, 0, 0)
      )
    }
    gc(verbose = FALSE)
  }

  if (!misc.counts) {
    message("remove 'obj@misc$counts' ")
    object@misc$counts <- NULL
  }
  object
}

#' @export
#' @method getFeaturesID Seurat
getFeaturesID.Seurat <- function(object, features, uniq = FALSE, ...) {
  object <- CheckSeuratFData(object)
  rowData <- object@misc$rowData
  if (!is.data.frame(rowData)) {
    fastWarning("No 'rowData' in this Seurat object. Cannot get row names")
    return(features)
  }
  getFeaturesID(object = rowData, features = features, uniq = uniq, ...)
}

#' @importFrom Eula.utils getFeaturesName
#' @export
#' @method getFeaturesName Seurat
getFeaturesName.Seurat <- function(
    object,
    features,
    col = "unique_name",
    uniq = FALSE,
    ...
) {
  object <- CheckSeuratFData(object)
  rowData <- object@misc$rowData
  if (!is.data.frame(rowData)) {
    fastWarning("No 'rowData' in this Seurat object. Cannot get feature names")
    return(features)
  }
  getFeaturesID(
    object = rowData,
    features = features,
    col = col,
    uniq = uniq,
    ...
  )
}
