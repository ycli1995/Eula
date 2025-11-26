
#' @importFrom Eula.utils dim_plot
#' @export
#' @method dim_plot Seurat
dim_plot.Seurat <- function(
    object,
    reduction = NULL,
    dims = c(1, 2),
    cells = NULL,

    group.by = NULL,
    split.by = NULL,
    shape.by = NULL,
    ncol = NULL,
    key = NULL,

    order = NULL,
    shuffle = FALSE,
    seed = 42,
    combine = TRUE,

    label = NULL,
    label.repel = TRUE,
    label.size = 4,
    label.color = "black",
    label.box = FALSE,

    colors = NULL,
    pt.size = NULL,
    pt.alpha = NULL,
    raster = NULL,
    raster.dpi = c(512, 512),
    facet.args = list(),
    corner.axis = FALSE,
    theme = NULL,
    ...
) {
  reduction <- reduction %||% DefaultDimReduc(object)
  reduc.obj <- object[[reduction]]
  cells <- cells %||% Cells(object, assay = DefaultAssay(reduc.obj))
  data <- FetchData(
    object,
    vars = paste0(Key(reduc.obj), dims),
    cells = cells,
    clean = "project"
  )
  key <- key %||% Key(reduc.obj)

  if (!is.null(shape.by)) {
    shape.by <- shape.by[1]
    data[, shape.by] <- object[[shape.by, drop = TRUE]]
  }
  if (!is.null(split.by)) {
    split.by <- split.by[1]
    split <- FetchData(object, vars = split.by, clean = TRUE)[split.by]
    data <- data[rownames(split), ]
    data[, split.by] <- split
  }

  orig.groups <- group.by
  group.by <- group.by %||% "ident"
  group.data <- FetchData(
    object,
    vars = group.by,
    cells = cells,
    clean = "project"
  )
  group.by <- colnames(group.data)
  group.data <- group.data[rownames(data), , drop = FALSE]
  data <- cbind(data, group.data)

  dim_plot(
    data,
    group.by = group.by,
    split.by = split.by,
    shape.by = shape.by,

    dims = c(1, 2),
    ncol = ncol,
    key = key,

    order = order,
    shuffle = shuffle,
    seed = seed,
    combine = combine,

    label = label,
    label.repel = label.repel,
    label.size = label.size,
    label.color = label.color,
    label.box = label.box,

    colors = colors,
    pt.size = pt.size,
    pt.alpha = pt.alpha,
    raster = raster,
    raster.dpi = raster.dpi,
    facet.args = facet.args,
    corner.axis = corner.axis,
    theme = theme,
    ...
  )
}

#' @importFrom Eula.utils feature_dim_plot
#' @export
#' @method feature_dim_plot Seurat
feature_dim_plot.Seurat <- function(
    object,
    features,
    reduction = NULL,
    dims = c(1, 2),
    cells = NULL,

    assay = NULL,
    slot = "data",

    split.by = NULL,
    ncol = NULL,
    key = NULL,

    order = NULL,
    shuffle = FALSE,
    seed = 42,
    combine = TRUE,

    colors = NULL,
    pt.size = NULL,
    pt.alpha = NULL,
    raster = NULL,
    raster.dpi = c(512, 512),
    facet.args = list(),
    corner.axis = FALSE,
    theme = NULL,
    ...
) {
  reduction <- reduction %||% DefaultDimReduc(object)
  reduc.obj <- object[[reduction]]
  assay <- assay %||% DefaultAssay(object)
  DefaultAssay(object) <- assay

  cells <- cells %||% Cells(object, assay = DefaultAssay(reduc.obj))
  data <- FetchData(
    object,
    vars = c(paste0(Key(reduc.obj), dims), features),
    cells = cells,
    clean = "project"
  )
  key <- key %||% Key(reduc.obj)

  if (!is.null(split.by)) {
    split.by <- split.by[1]
    split <- FetchData(object, vars = split.by, clean = TRUE)[split.by]
    data <- data[rownames(split), ]
    data[, split.by] <- split
  }
  feature_dim_plot(
    object = data,
    features = features,

    split.by = split.by,

    dims = dims,
    ncol = ncol,
    key = key,

    order = order,
    shuffle = shuffle,
    seed = seed,
    combine = combine,

    colors = colors,
    pt.size = pt.size,
    pt.alpha = pt.alpha,
    raster = raster,
    raster.dpi = raster.dpi,
    facet.args = facet.args,
    corner.axis = corner.axis,
    theme = theme,
    ...
  )
}

#' @importFrom Eula.utils dot_plot
#' @export
#' @method dot_plot Seurat
dot_plot.Seurat <- function(
    object,
    features,
    assay = NULL,
    cells = NULL,
    group.by = NULL,
    split.by = NULL,
    mean.fxn = NULL,
    min.exp = 0,
    scale = TRUE,
    split.features = NULL,
    colors = NULL,
    color.limits = c(-2.5, 2.5),
    size.limits = c(0, 100),
    coord.flip = FALSE,
    theme = NULL,
    ...
) {
  if (length(features) == 0) {
    stop("No 'features' passed to `dot_plot`.")
  }
  assay <- assay %||% DefaultAssay(object)
  DefaultAssay(object) <- assay

  orig.groups <- group.by
  group.by <- group.by %||% "ident"
  group.data <- FetchData(
    object,
    vars = group.by,
    cells = cells,
    clean = "project"
  )
  group.data$group.by <- as.factor(group.data[, 1])

  if (length(split.by) > 0) {
    split.data <- FetchData(
      object,
      vars = split.by,
      cells = cells,
      clean = "project"
    )
    group.data$split.by <- as.factor(split.data[, 1])
  }

  data <- GetAssayData(object)[features, , drop = FALSE]
  if (!setequal(rownames(group.data), colnames(data))) {
    data <- data[, rownames(group.data), drop = FALSE]
  }
  if (!all(rownames(data) == features)) {
    data <- data[features, , drop = FALSE]
  }
  group.data <- group.data[colnames(data), , drop = FALSE]

  dot_plot(
    object = data,
    group.by = group.data$group.by,
    split.by = group.data$split.by,
    split.features = split.features,
    mean.fxn = mean.fxn,
    min.exp = min.exp,
    scale = scale,
    colors = colors,
    color.limits = color.limits,
    size.limits = size.limits,
    coord.flip = coord.flip,
    theme = theme,
    ...
  )
}
