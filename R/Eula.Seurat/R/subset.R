#' @importFrom SeuratObject Assays Cells Graphs Images
NULL

#' @importFrom Eula.utils getArgList
get_cells_and_features <- function(
    object,
    cells = NULL,
    cells_exclude = NULL,
    features = NULL,
    features_exclude = NULL
) {
  all_cells <- rownames(object[[]])

  cells <- getArgList(cells)
  features <- getArgList(features)
  cells_exclude <- getArgList(cells_exclude)
  features_exclude <- getArgList(features_exclude)

  if (length(cells) > 0) {
    all_cells <- intersect(all_cells, cells)
  }
  all_cells <- setdiff(all_cells, cells_exclude)

  all_features <- rownames(object)
  if (length(features) > 0) {
    all_features <- intersect(all_features, features)
  }
  all_features <- setdiff(all_features, features_exclude)

  if (setequal(all_cells, Cells(object))) {
    all_cells <- NULL
  }
  if (setequal(all_features, rownames(object))) {
    all_features <- NULL
  }
  list(cells = all_cells, features = all_features)
}
#
# get_subset_cells <- function(object, column, include = NULL, exclude = NULL) {
#   cells <- rownames(object[[]])
#   if (!column %in% colnames(object@meta.data)) {
#     fastWarning("Subset: '", column, "' not found in object@meta.data")
#     return(cells)
#   }
#   if (length(include) > 0) {
#     include <- as.character(include)
#     message("Search '", column, "' includes: ", paste(include, collapse = ", "))
#     use_cells <- rownames(object[[]])[object@meta.data[, column] %in% include]
#     cells <- intersect(cells, use_cells)
#   }
#   if (length(exclude) > 0) {
#     exclude <- as.character(exclude)
#     message("Search '", column, "' excludes: ", paste(exclude, collapse = ", "))
#     use_cells <- rownames(object[[]])[!object@meta.data[, column] %in% exclude]
#     cells <- intersect(cells, use_cells)
#   }
#   return(cells)
# }

check_image_cells <- function(object) {
  if (!"images" %in% slotNames(object)) {
    return(object)
  }
  for (image in Images(object)) {
    image.cells <- intersect(Cells(object[[image]]), rownames(object[[]]))
    if (length(image.cells) > 0) {
      object@images[[image]] <- object[[image]][image.cells]
    } else {
      object@images[[image]] <- NULL
    }
  }
  return(object)
}

#' @importFrom Eula.utils filterData
#' @export
SubsetObject <- function(
    object,
    column_map = list(),
    cells = NULL,
    features = NULL,
    ...
) {
  cells_kp <- rownames(object[[]])
  for (i in names(column_map)) {
    if (!is.null(column_map[[i]])) {
      invert <- column_map[[i]][['invert']] %||% FALSE
      include <- column_map[[i]][['include']]
      if (length(column_map[[i]][['exclude']]) > 0) {
        include <- column_map[[i]][['exclude']]
        invert <- TRUE
      }
      filter.out <- filterData(
        x = object[[]],
        column = i,
        include = include,
        invert = invert
      )
      use_cells <- rownames(object[[]])[filter.out]
      cells_kp <- intersect(cells_kp, use_cells)
    }
  }
  if (length(cells) > 0) {
    cells_kp <- intersect(cells_kp, cells)
  }

  ## Keep old images
  if ('images' %in% slotNames(object)) {
    old_images <- Images(object)
  }

  ## Avoid subset graphs with only 1 cell.
  for (i in SeuratObject::Graphs(object)) {
    if (length(intersect(cells_kp, Cells(object[[i]]))) < 2) {
      object[[i]] <- NULL
    }
  }

  if (length(cells_kp) > 0 & !setequal(cells_kp, rownames(object[[]]))) {
    object <- subset(object, cells = cells_kp)
  }
  if (length(features) > 0) {
    object <- subset(object, features = features)
  }
  object@meta.data <- droplevels(object@meta.data)

  ## For assays whose cells are not matched with object[[]]
  for (i in SeuratObject::Assays(object)) {
    if (inherits(object[[i]], "StdAssay")) {
      assay_cells <- intersect(Cells(object@assays[[i]]), rownames(object[[]]))
      object@assays[[i]]@cells <- object@assays[[i]]@cells[assay_cells, ]
      if (nrow(object[[i]]@cells) == 0) {
        object[[i]] <- NULL
      }
    }
  }

  if ('images' %in% slotNames(object)) {
    ## when object[, cells], if image's name is like 'A-B', it will be
    ## duplicated with name 'A.B', here substracting right name images
    ## with 'sample.name' in object@meta.data
    object@images <- object@images[Images(object) %in% old_images]
    object <- check_image_cells(object)
  }

  ## Colors in old version
  old.colors <- list(
    seurat_clusters = "color.cluster",
    orig.ident = "color.sample",
    Groups = "color.group"
  )
  for (i in names(old.colors)) {
    misc.name <- old.colors[[i]]
    object <- checkColorMap(
      object,
      column = i,
      colors = object@misc[[misc.name]],
      misc.name = misc.name
    )
  }

  for (i in names(object@misc$colors)) {
    object@misc$colors[[i]] <- checkColorMap(
      object[[]], i,
      colors = object@misc$colors[[i]]
    )
  }
  return(object)
}

.SUBSET_FILES <- c(
  "cells_use",
  "cells_exclude",
  "features_use",
  "features_exclude"
)

#' @importFrom Eula.utils readTable
read_subset_files <- function(parameter, ...) {
  for (i in .SUBSET_FILES) {
    if (!is.null(parameter[[i]])) {
      parameter[[i]] <- readTable(parameter[[i]], header = FALSE, ...)[, 1]
    }
  }
  return(parameter)
}

#' @export
SubsetObjectWrapper <- function(object, parameter = list(), ...) {
  parameter <- read_subset_files(parameter)
  pset1 <- .SUBSET_FILES
  param1 <- parameter[pset1]

  parameter[['check_color']] <- NULL

  param2 <- parameter[setdiff(names(parameter), pset1)]
  cells_features <- get_cells_and_features(
    object = object,
    cells = param1[["cells_use"]],
    cells_exclude = param1[["cells_exclude"]],
    features = param1[["features_use"]],
    features_exclude = param1[["features_exclude"]]
  )
  object <- SubsetObject(
    object,
    column_map = param2,
    cells = cells_features[['cells']],
    features = cells_features[['features']]
  )
  object
}

