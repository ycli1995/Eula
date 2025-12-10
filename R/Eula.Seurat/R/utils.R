#' @include old.R
NULL

#' @importFrom Eula.utils readRDX
#' @importFrom SeuratObject UpdateSeuratObject
#' @export
Load <- function(file) {
  object <- Eula.utils::readRDX(file)
  if (!is(object, "Seurat")) {
    return(object)
  }
  if (!"version" %in% slotNames(object)) {
    return(object)
  }
  if (grepl('^2', object@version)) {
    return(SeuratObject::UpdateSeuratObject(object))
  }
  CheckMySeuratObj(object)
}

#' @export ModuleScoring
ModuleScoring <- function(object, ...) {
  UseMethod("ModuleScoring", object)
}

#' @importFrom ggplot2 cut_number
#' @export
#' @method ModuleScoring default
ModuleScoring.default <- function(
    object,
    features,
    pool = NULL,
    nbin = 24,
    ctrl = 100,
    seed = 42,
    name = 'Cluster',
    ...
) {
  if (length(features) == 0) {
    stop("Missing input feature list")
  }
  if (is.character(features)) {
    features <- list(features)
    names(features) <- paste0(name, "1")
  }
  old.features <- features
  upper.rownames <- toupper(rownames(object))
  for (i in seq_along(features)) {
    features[[i]] <- intersect(toupper(features[[i]]), upper.rownames)
    which.missing <- which(!toupper(old.features[[i]]) %in% upper.rownames)
    if (length(which.missing) > 0) {
      which.missing <- paste(old.features[[i]][which.missing], collapse = ", ")
      fastWarning("The following features are not found:\n ", which.missing)
    }
    if (length(features[[i]]) == 0) {
      fastWarning("'features[[", i, "]]' contains no valid features.")
    }
  }
  features <- features[lengths(features) > 0]
  if (length(features) == 0) {
    stop("No valid features found.")
  }

  rownames(object) <- upper.rownames

  data.avg <- Matrix::rowMeans(object)
  data.avg <- data.avg[order(data.avg)]
  data.cut <- cut_number(
    x = data.avg + rnorm(n = length(data.avg)) / 1e30,
    n = nbin,
    labels = FALSE,
    right = FALSE
  )
  names(data.cut) <- names(data.avg)
  ctrl.use <- vector(mode = "list", length = length(features))
  for (i in seq_along(features)) {
    features.use <- features[[i]]
    for (j in seq_along(features.use)) {
      bin <- data.cut[features.use[j]]
      set.seed(seed)
      features.ctrl <- names(sample(
        x = data.cut[data.cut == bin],
        size = ctrl,
        replace = FALSE
      ))
      ctrl.use[[i]] <- c(ctrl.use[[i]], features.ctrl)
    }
    ctrl.use[[i]] <- unique(ctrl.use[[i]])
  }
  ctrl.scores <- matrix(0, nrow = length(ctrl.use), ncol = ncol(object))
  for (i in seq_along(ctrl.use)) {
    ctrl.scores[i, ] <- Matrix::colMeans(object[ctrl.use[[i]], , drop = FALSE])
  }
  features.scores <- matrix(0, nrow = length(features), ncol = ncol(object))
  for (i in seq_along(features)) {
    features.scores[i, ] <- Matrix::colMeans(object[features[[i]], ])
  }
  features.scores <- features.scores - ctrl.scores

  if (length(names(features)) == 0) {
    names(features) <- paste0(name, seq_along(features))
  }
  rownames(features.scores) <- names(features)
  features.scores <- as.data.frame(t(features.scores))
  rownames(features.scores) <- colnames(object)
  features.scores
}

#' @importClassesFrom SeuratObject Assay
#' @export
#' @method ModuleScoring Assay
ModuleScoring.Assay <- function(
    object,
    features,
    slot = "data",
    pool = NULL,
    nbin = 24,
    ctrl = 100,
    seed = 42,
    name = 'Cluster',
    ...
) {
  ModuleScoring(
    object = GetAssayData(object, slot),
    features = features,
    pool = pool,
    nbin = nbin,
    ctrl = ctrl,
    seed = seed,
    name = name,
    ...
  )
}

#' @importClassesFrom SeuratObject StdAssay
#' @export
#' @method ModuleScoring StdAssay
ModuleScoring.StdAssay <- function(
    object,
    features,
    slot = "data",
    pool = NULL,
    nbin = 24,
    ctrl = 100,
    seed = 42,
    name = 'Cluster',
    ...
) {
  layer_names <- SeuratObject::Layers(object, search = slot)
  output_list <- lapply(layer_names, function(x) {
    ModuleScoring(
      object = LayerData(object, layer = x),
      features = features,
      pool = pool,
      nbin = nbin,
      ctrl = ctrl,
      seed = seed,
      name = name,
      ...
    )
  })
  output_list <- unname(output_list)
  features.scores.use <- do.call(rbind, output_list)
  common.cells <- intersect(Cells(object), rownames(features.scores.use))
  features.scores.use[common.cells, , drop = FALSE]
}

#' @importClassesFrom SeuratObject Seurat
#' @export
#' @method ModuleScoring Seurat
ModuleScoring.Seurat <- function(
    object,
    features,
    assay = NULL,
    slot = "data",
    pool = NULL,
    nbin = 24,
    ctrl = 100,
    seed = 42,
    name = 'Cluster',
    ...
) {
  assay.old <- DefaultAssay(object = object)
  assay <- assay %||% assay.old
  DefaultAssay(object) <- assay
  features.scores.use <- ModuleScoring(
    object = object[[assay]],
    features = features,
    slot = slot,
    pool = pool,
    nbin = nbin,
    ctrl = ctrl,
    seed = seed,
    name = name,
    ...
  )
  object[[colnames(features.scores.use)]] <- features.scores.use
  DefaultAssay(object) <- assay.old
  object
}

#' @export CellCycle
CellCycle <- function(object, ...) {
  UseMethod("CellCycle", object)
}

#' @export
#' @method CellCycle default
CellCycle.default <- function(
    object,
    s.features,
    g2m.features,
    ctrl = NULL,
    ...
) {
  features <- list('S.Score' = s.features, 'G2M.Score' = g2m.features)
  ctrl <- ctrl %||% min(vapply(features, FUN = length, FUN.VALUE = numeric(1)))
  cc <- ModuleScoring(object = object, features = features, ctrl = ctrl, ...)
  cc[['CC.Difference']] <- cc[['S.Score']] - cc[['G2M.Score']]

  assignments <- apply(
    X = cc,
    MARGIN = 1,
    FUN = function(scores, first = 'S', second = 'G2M', null = 'G1') {
      if (all(scores < 0)) {
        return(null)
      }
      if (length(which(scores == max(scores))) > 1) {
        return('Undecided')
      }
      return(c(first, second)[which.max(scores)])
    }
  )
  cc[["Phase"]] <- factor(assignments, c("G1", "S", "G2M"))
  cc
}

#' @export
#' @method CellCycle Assay
CellCycle.Assay <- function(
    object,
    s.features,
    g2m.features,
    slot = "data",
    ctrl = NULL,
    ...
) {
  CellCycle(
    object = GetAssayData(object, slot),
    s.features = s.features,
    g2m.features = g2m.features,
    slot = slot,
    ctrl = ctrl,
    ...
  )
}

#' @export
#' @method CellCycle StdAssay
CellCycle.StdAssay <- function(
    object,
    s.features,
    g2m.features,
    slot = "data",
    ctrl = NULL,
    ...
) {
  layer_names <- SeuratObject::Layers(object, search = slot)
  output_list <- lapply(layer_names, function(x) {
    CellCycle(
      object = SeuratObject::LayerData(object, layer = x),
      s.features = s.features,
      g2m.features = g2m.features,
      ctrl = ctrl,
      ...
    )
  })
  output_list <- unname(output_list)
  features.scores.use <- do.call(rbind, output_list)
  common.cells <- intersect(Cells(object), rownames(features.scores.use))
  features.scores.use[common.cells, , drop = FALSE]
}

#' @export
#' @importFrom SeuratObject DefaultAssay<-
#' @method CellCycle Seurat
CellCycle.Seurat <- function(
    object,
    s.features,
    g2m.features,
    assay = NULL,
    slot = "data",
    ctrl = NULL,
    set.ident = FALSE,
    ...
) {
  assay.old <- DefaultAssay(object = object)
  assay <- assay %||% assay.old
  DefaultAssay(object) <- assay
  features.scores.use <- CellCycle(
    object = object[[assay]],
    s.features = s.features,
    g2m.features = g2m.features,
    ctrl = ctrl,
    slot = slot,
    ...
  )
  object <- AddMetaData(object, metadata = features.scores.use)
  DefaultAssay(object) <- assay.old
  if (set.ident) {
    Idents(object) <- "Phase"
  }
  object
}

#' @importFrom Eula.utils validCharacters
#' @importFrom SeuratObject FetchData
#' @importFrom dplyr mutate rename
#' @importFrom tidyselect all_of
#' @export
FetchSeuratData <- function(object, vars, ...) {
  new.names <- names(vars)
  df <- FetchData(object, vars = vars, ...)
  if (any(validCharacters(new.names))) {
    vars <- vars[validCharacters(new.names)]
    vars <- vars[vars %in% colnames(df)]
    df <- df %>%
      dplyr::rename(tidyselect::all_of(vars))
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

#' @importFrom Matrix Diagonal
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
  dd <- Matrix::Diagonal(x = 1 / scale.factor * nCounts, names = colnames(mat))
  mat <- round(expm1(mat) %*% dd)
  mat
}

#' @importFrom SeuratObject Assays DefaultAssay DefaultAssay<-
#' @export
ShrinkSeuratObject <- function(
    object,
    assays = NULL,
    scale.data = NULL,
    misc.counts = NULL,
    verbose = TRUE,
    ...
) {
  assays <- assays %||% SeuratObject::Assays(object)
  scale.data <- scale.data %||% FALSE
  misc.counts <- misc.counts %||% TRUE

  if (!all(SeuratObject::Assays(object) %in% assays)) {
    verboseMsg("keep assays: ", paste(assays, collapse = ", "))
    if (!DefaultAssay(object) %in% assays) {
      DefaultAssay(object) <- assays[1]
    }
    object@assays <- object@assays[assays]
  }

  if (!scale.data) {
    verboseMsg("remove 'scale.data' ")
    if (packageVersion("Seurat") >= package_version("5.0.0")) {
      if (any(dim(SeuratObject::LayerData(object, layer = "scale.data")) > 0)) {
        SeuratObject::LayerData(object, layer = "scale.data") <- NULL
      }
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
  getFeaturesName(
    object = rowData,
    features = features,
    col = col,
    uniq = uniq,
    ...
  )
}
