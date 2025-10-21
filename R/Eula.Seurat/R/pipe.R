
#' @export
pipe_RenameObject <- function(obj, parameter = list()) {

  print(str(parameter))

  obj <- RenameSeuratColumns(obj, parameter[['rename_colnames']])

  obj <- RenameSeuratWrapper(obj, parameter[['rename']])

  obj <- pipe_UpdateSeurat(obj, parameter[['UpdateSeurat']])

  obj <- CheckMySeuratObj(obj, parameter$default_colnames)

  obj <- pipe_SubsetObject(obj, parameter[['subset']])

  obj <- pipe_ShrinkSeuratObject(obj, parameter)

  Message('>>>>> Finally check...')
  print(obj)
  print(str(obj@meta.data))
  print(obj@misc$colors)

  obj
}

#' @export
pipe_ShrinkSeuratObject <- function(obj, parameter = list(), ...) {
  parameter$scale.data <- parameter$scale.data %||% FALSE
  parameter$misc.counts <- parameter$misc.counts %||% TRUE

  parameter$assays <- parameter$assays %||% Assays(obj)
  parameter$assays <- intersect(parameter$assays, Assays(obj))

  Message('>>>>> Shrinking Seurat object...')
  obj <- ShrinkSeuratObject(
    object = obj,
    assays = parameter$assays,
    scale.data = parameter$scale.data,
    misc.counts = parameter$misc.counts
  )
  obj
}

#' @export
pipe_SubsetObject <- function(obj, parameter = list(), ...) {
  parameter$cells <- parameter$cells %||% parameter$cells_use
  parameter$features <- parameter$features %||% parameter$features

  Message('>>>>> Subset Seurat object...')
  cells_features <- get_cells_and_features(
    object = obj,
    cells = parameter$cells,
    features = parameter$features,
    cells_exclude = parameter$cells_exclude,
    features_exclude = parameter$features_exclude
  )

  parameter <- parameter[intersect(names(parameter), colnames(obj[[]]))]
  obj <- SubsetObject(
    object = obj,
    column_map = parameter,
    cells = cells_features[['cells']],
    features = cells_features[['features']]
  )
  obj
}

#' @export
pipe_UpdateSeurat <- function(obj, parameter = list(), ...) {
  parameter$Assay5 <- parameter$Assay5 %||% FALSE

  Message('>>>>> Check Seurat version...')
  obj <- UpdateSeuratAll(object = obj, Assay5 = parameter$Assay5)
  obj
}

#' @export
pipe_StatFeaturePercentage <- function(obj, parameter = list(), ...) {
  stat.pct <- parameter[['stat.pct']] %||% TRUE
  parameter[['stat.pct']] <- NULL

  Message('>>>>> Stat features')
  for (i in names(parameter)) {
    i2 <- gsub("\\/| ", "_", i)
    obj <- StatFeatures(
      object = obj,
      features = parameter[[i]],
      col.name = paste0("percent.", i2),
      stat.pct = stat.pct
    )
  }
  obj
}

#' @export
pipe_MakeSeuratObj <- function(parameter) {
  if (length(parameter) == 0) {
    stop("'MakeSeuratObj' parameter is empty.")
  }
  data.paths <- parameter[["paths"]]
  if (length(data.paths) == 0) {
    stop("No 'paths' for matrices.")
  }
  data.names <- parameter[["names"]]
  if (length(data.names) != length(data.paths)) {
    stop("'names' must has the same length as 'paths'")
  }
  name.list <- parameter[['name_list']]
  if (length(name.list) == 0) {
    stop("No 'name_list' for MakeSeuratObj")
  }
  assay <- parameter[['assay']] %||% "RNA"

  Message('>>>>> Making Seurat object')
  obj <- MakeSeuratObj(
    data.paths = data.paths,
    data.names = data.names,
    assay = assay,
    name.list = name.list
  )
  obj
}

#' @export
pipe_MakeSeuratObj_R <- function(parameter) {
  obj <- pipe_MakeSeuratObj(parameter[["MakeSeuratObj"]])

  obj <- pipe_StatFeaturePercentage(obj, parameter[["StatFeaturePercentage"]])

  obj <- RenameSeuratColumns(obj, parameter[['rename_colnames']])

  obj <- RenameSeuratWrapper(obj, parameter[['rename']])

  obj <- pipe_UpdateSeurat(obj, parameter[['UpdataSeurat']])

  obj <- CheckMySeuratObj(obj, parameters = parameter[['default_colnames']])

  obj <- pipe_SubsetObject(obj, parameter[['subset']])
}
