
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
  key <- key %||% Key(reduc.obj)
  data <- .fetch_seurat_data(
    object = object,
    reduction = reduction,
    dims = dims,
    cells = cells,
    features = NULL,
    group.by = group.by,
    split.by = split.by,
    shape.by = shape.by
  )
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
    key = NULL,

    order = NULL,
    shuffle = FALSE,
    seed = 42,
    combine = TRUE,
    ncol = NULL,

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
  if (length(features) == 0) {
    stop("No 'features' passed to `feature_dim_plot`.")
  }
  reduction <- reduction %||% DefaultDimReduc(object)
  reduc.obj <- object[[reduction]]
  key <- key %||% Key(reduc.obj)
  data <- .fetch_seurat_data(
    object = object,
    reduction = reduction,
    dims = dims,
    cells = cells,
    features = features,
    group.by = NULL,
    split.by = split.by,
    shape.by = NULL
  )
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
  group.by <- group.by[1]
  data <- .fetch_seurat_data(
    object = object,
    reduction = NULL,
    cells = cells,
    assay = assay,
    features = features,
    group.by = group.by,
    split.by = split.by,
    shape.by = NULL
  )
  group.data <- data[, setdiff(colnames(data), features), drop = FALSE]
  if (length(split.by) > 0) {
    group.data[["split.by"]] <- group.data[[split.by]]
  }
  data <- t(as.matrix(data[, features, drop = FALSE]))
  if (!is.numeric(data)) {
    stop("The input features must be numeric.")
  }
  dot_plot(
    object = data,
    group.by = group.data[[1]],
    split.by = group.data[["split.by"]],
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

#' @importFrom Eula.utils violin_plot
#' @export
#' @method violin_plot Seurat
violin_plot.Seurat <- function(
    object,
    features,
    group.by = NULL,
    split.by = NULL,
    cells = NULL,

    assay = NULL,
    slot = "data",

    combine = TRUE,

    colors = NULL,
    noise.scale = 1e+5,
    scale = "width",
    adjust = 1,

    box.width = 0,
    box.color = "black",
    box.dodge.width = 0.9,
    pt.size = 0,
    pt.alpha = 1,
    jitter.width = 0.3,
    jitter.dodge.width = 0.8,
    raster = NULL,
    raster.dpi = c(512, 512),
    coord.flip = FALSE,
    facet.args = list(),
    theme = NULL,
    ...
) {
  if (length(features) == 0) {
    stop("No 'features' passed to `violin_plot`.")
  }
  data <- .fetch_seurat_data(
    object = object,
    reduction = NULL,
    cells = cells,
    features = features,
    group.by = NULL,
    split.by = split.by,
    shape.by = NULL
  )
  group.by <- colnames(data)[1]
  violin_plot(
    object = data,
    features = features,
    group.by = group.by,
    split.by = split.by,

    combine = combine,

    colors = colors,
    noise.scale = noise.scale,
    scale = scale,
    adjust = adjust,

    box.width = box.width,
    box.color = box.color,
    box.dodge.width = box.dodge.width,
    pt.size = pt.size,
    pt.alpha = pt.alpha,
    jitter.width = jitter.width,
    jitter.dodge.width = jitter.dodge.width,
    raster = raster,
    raster.dpi = raster.dpi,
    coord.flip = coord.flip,
    facet.args = facet.args,
    theme = theme,
    ...
  )
}

.fetch_seurat_data <- function(
    object,
    reduction = NULL,
    dims = c(1, 2),
    assay = NULL,
    cells = NULL,
    features = NULL,
    group.by = NULL,
    split.by = NULL,
    shape.by = NULL
) {
  assay <- assay %||% DefaultAssay(object)
  group.by <- group.by %||% "ident"
  split.by <- split.by[1]
  shape.by <- shape.by[1]
  vars <- c(group.by, split.by, shape.by, features)
  if (!is.null(reduction)) {
    reduc.obj <- object[[reduction]]
    cells <- cells %||% Cells(object, assay = DefaultAssay(reduc.obj))
    vars <- c(paste0(Key(reduc.obj), dims), vars)
  } else {
    cells <- cells %||% Cells(object, assay = assay)
  }
  FetchData(object, vars = vars, cells = cells)
}
