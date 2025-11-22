
#' @export dim_plot
dim_plot <- function(object, ...) {
  UseMethod("dim_plot", object)
}

#' @importFrom patchwork wrap_plots
#' @importFrom Eula.utils theme_base_default theme_facet theme_no_legend
#'
#' @export
#' @method dim_plot data.frame
dim_plot.data.frame <- function(
    object,
    group.by,
    dims = c(1, 2),
    split.by = NULL,
    shape.by = NULL,
    ncol = NULL,
    key = "Dim",

    order = NULL,
    shuffle = FALSE,
    seed = 42,
    combine = TRUE,

    label = NULL,
    label.repel = TRUE,
    label.size = 4,
    label.color = "black",
    label.box = FALSE,

    legend = NULL,

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
  orig.groups <- group.by
  group.by <- intersect(orig.groups, colnames(object))
  if (length(group.by) == 0) {
    stop("Invalid 'group.by':\n  ", paste(orig.groups, collapse = ", "))
  }

  legend <- legend %||% group.by

  colnames(object)[dims] <- c("x", "y")

  split.by <- split.by[1]
  shape.by <- shape.by[1]

  if (isTRUE(shuffle)) {
    set.seed(seed)
    object <- object[sample(seq_len(nrow(object))), , drop = FALSE]
  }

  if (is.character(colors)) {
    colors <- sapply(group.by, function(x) colors, simplify = FALSE)
  }

  labs.args <- setNames(paste0(key, dims), c("x", "y"))
  labs.args <- as.list(labs.args)

  theme <- theme %||% theme_dim_default()

  ncol <- ncol %||% 4
  if (!is.null(split.by)) {
    theme <- theme + theme_facet()
    if (length(group.by) > 1) {
      facet.args[['nrow']] <- 1
      facet.args[['ncol']] <- NULL
    }
  }

  pt.alpha <- pt.alpha %||% 1

  plist <- list()
  for (i in group.by) {
    data <- object[, c("x", "y")]

    data[, 'colour'] <- as.factor(object[, i])
    cols <- colors[[i]]

    label.data <- NULL
    if (i %in% label) {
      label.data <- get_label_data(data)
    }
    labs.args[["colour"]] <- i

    if (!is.null(split.by)) {
      data[, 'split'] <- object[, split.by]
    }
    if (!is.null(shape.by)) {
      data[, 'shape'] <- object[, shape.by]
    }

    if (isTRUE(order)) {
      data <- data[order(data[['colour']]), , drop = FALSE]
    }
    plist[[i]] <- single_dim_plot(
      data = data,
      colors = cols,
      label.data = label.data,
      label.repel = label.repel,
      label.size = label.size,
      label.color = label.color,
      label.box = label.box,
      pt.size = pt.size,
      raster = raster,
      raster.dpi = raster.dpi,
      facet.args = facet.args,
      labs.args = labs.args,
      theme = theme,
      corner.axis = corner.axis,
      alpha = pt.alpha,
      ...
    )
    if (!i %in% legend) {
      plist[[i]] <- plist[[i]] + theme_no_legend()
    }
  }
  if (!is.null(split.by)) {
    ncol <- 1
  }
  if (combine) {
    ncol <- min(ncol, length(plist))
    plist <- wrap_plots(plist, ncol = ncol)
  }
  plist
}

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

#' @export feature_dim_plot
feature_dim_plot <- function(object, ...) {
  UseMethod("feature_dim_plot", object)
}

#' @importFrom Eula.utils theme_dim_default theme_facet
#' @export
#' @method feature_dim_plot data.frame
feature_dim_plot.data.frame <- function(
    object,
    feature.data,
    split.by = NULL,

    dims = c(1, 2),
    ncol = NULL,
    key = "Dim",

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
  if (!is.matrix(feature.data)) {
    stop("'feature.data' must be a matrix.")
  }
  if (nrow(feature.data) == 0) {
    stop("No feature in `feature.data`.")
  }
  if (nrow(object) != ncol(feature.data)) {
    stop("`nrow(object)` must be the same as `ncol(feature.data)`.")
  }
  if (!setequal(rownames(object), colnames(feature.data))) {
    stop("`rownames(object)` must be the same as `colnames(feature.data)`")
  }

  if (isTRUE(shuffle)) {
    set.seed(seed)
    object <- object[sample(seq_len(nrow(object))), , drop = FALSE]
  }

  colnames(object)[dims] <- c("x", "y")
  labs.args <- setNames(paste0(key, dims), c("x", "y"))
  labs.args <- as.list(labs.args)
  labs.args[['colour']] <- "Exp"

  colors <- colors %||% c(
    "#FFF7EC", "#FEE8C8", "#FDD49E", "#FDBB84", "#FC8D59",
    "#EF6548", "#D7301F", "#B30000", "#7F0000"
  )
  theme <- theme %||% theme_dim_default()

  ncol <- ncol %||% 4
  split.by <- split.by[1]
  if (!is.null(split.by)) {
    theme <- theme + theme_facet()
    facet.args[['nrow']] <- 1
    facet.args[['ncol']] <- NULL
    ncol <- 1
  }

  pt.alpha <- pt.alpha %||% 1

  features <- rownames(feature.data)
  feature.data <- t(feature.data[, rownames(object), drop = FALSE])
  plist <- list()
  for (i in features) {
    data <- object[, c("x", "y")]
    data[, 'colour'] <- feature.data[, i]
    labs.args[['title']] <- i
    if (!is.null(split.by)) {
      data[, 'split'] <- object[, split.by]
    }
    if (isTRUE(order)) {
      data <- data[order(data[['colour']]), , drop = FALSE]
    }
    plist[[i]] <- single_dim_plot(
      data = data,
      colors = colors,
      pt.size = pt.size,
      raster = raster,
      raster.dpi = raster.dpi,
      facet.args = facet.args,
      labs.args = labs.args,
      theme = theme,
      corner.axis = corner.axis,
      alpha = pt.alpha,
      ...
    )
  }
  if (combine) {
    ncol <- min(ncol, length(plist))
    plist <- wrap_plots(plist, ncol = ncol)
  }
  plist
}

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
    vars = paste0(Key(reduc.obj), dims),
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
  feature.data <- as.matrix(FetchData(object, vars = features, cells = cells))
  feature_dim_plot(
    object = data,
    feature.data = t(feature.data),

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

#' @export dot_plot
dot_plot <- function(object, ...) {
  UseMethod("dot_plot", object)
}

#' @importFrom ggplot2 facet_grid
#' @importFrom Eula.utils theme_dot_default single_dot_plot
#' @export
#' @method dot_plot CsparseMatrix
dot_plot.CsparseMatrix <- function(
    object,
    group.by,
    split.by = NULL,
    split.features = NULL,
    mean.fxn = NULL,
    min.exp = 0,
    scale = TRUE,
    colors = NULL,
    color.limits = c(-2.5, 2.5),
    size.limits = c(0, 100),
    coord.flip = FALSE,
    theme = NULL,
    ...
) {
  mean.fxn <- mean.fxn %||% function(x) rowExpMean(x, log = TRUE)
  group.by <- as.factor(group.by)

  id <- group.by
  if (length(split.by) > 0) {
    if (length(split.by) != ncol(object)) {
      stop("'split.by' must have the same length as 'ncol(object)'.")
    }
    id <- droplevels(pasteFactors(group.by, split.by))
  }
  cell.groups <- split(colnames(object), f = id)
  out <- rowMeanPct(
    object = object,
    cell.groups = cell.groups,
    mean.fxn = mean.fxn,
    min.exp = min.exp
  )

  group.data <- data.frame(id = id, group.by = group.by)
  if (length(split.by) > 0) {
    group.data$split.by <- split.by
  }
  group.data <- group.data %>%
    group_by(id) %>%
    summarize(across(everything(), unique)) %>%
    as.data.frame()

  df <- .get_dot_plot_data(
    out$avg.exp,
    out$avg.pct,
    group.by = group.data$group.by,
    split.by = group.data$split.by,
    split.features = split.features
  )
  df$x <- df$group.by
  df$y <- df$features
  if (coord.flip) {
    df$x <- df$features
    df$y <- df$group.by
  }
  if (scale) {
    df$colour <- df$avg.scale.exp
  } else {
    df$colour <- df$avg.exp
  }
  if (length(color.limits) == 2) {
    df$colour <- min_max_cut(df$colour, limits = color.limits)
  }
  df$size <- df$avg.pct * 100
  labs.args <- list(
    x = "",
    y = "",
    size = "Expressed percentage",
    color = "Average Expression"
  )
  theme <- theme %||% theme_dot_default()
  p <- single_dot_plot(
    data = df,
    colors = colors,
    color.limits = color.limits,
    size.limits = size.limits,
    theme = theme,
    facet.args = list(),
    labs.args = labs.args,
    ...
  )
  if (coord.flip) {
    p <- p + theme(legend.box = "horizontal")
  }
  if (length(split.features) > 0 & length(split.by) > 0) {
    if (coord.flip) {
      return(p + facet_grid(split.by ~ split.features, scales = "free"))
    }
    return(p + facet_grid(split.features ~ split.by, scales = "free"))
  }
  if (length(split.features) > 0) {
    facet.layer <- if (coord.flip) {
      facet_wrap(
        "~split.features",
        nrow = 1,
        scales = "free_x",
        strip.position = "top"
      )
    } else {
      facet_wrap(
        "~split.features",
        ncol = 1,
        scales = "free_y",
        strip.position = "right"
      )
    }
    return(p + facet.layer)
  }
  if (length(split.by) > 0) {
    facet.layer <- if (coord.flip) {
      facet_wrap(
        "~split.by",
        ncol = 1,
        scales = "free_y",
        strip.position = "right"
      )
    } else {
      facet_wrap(
        "~split.by",
        nrow = 1,
        scales = "free_x",
        strip.position = "top"
      )
    }
    return(p + facet.layer)
  }
  p
}

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

#' @importFrom reshape2 melt
#' @importFrom dplyr left_join
.get_dot_plot_data <- function(
    avg.exp,
    avg.pct,
    group.by = NULL,
    split.by = NULL,
    split.features = NULL,
    ...
) {
  features <- rownames(avg.exp)
  clusters <- colnames(avg.pct)

  group.by <- group.by %||% clusters
  if (length(group.by) > 0) {
    if (length(group.by) != ncol(avg.exp)) {
      stop("'group.by' must have the same length as 'ncol(avg.exp)'.")
    }
    group.by <- as.factor(group.by)
    names(group.by) <- clusters
  }
  if (length(split.by) > 0) {
    if (length(split.by) != ncol(avg.exp)) {
      stop("'split.by' must have the same length as 'ncol(avg.exp)'.")
    }
    split.by <- as.factor(split.by)
    names(split.by) <- clusters
  }
  if (length(split.features) > 0) {
    if (length(split.features) != nrow(avg.exp)) {
      stop("'split.features' must have the same length as 'nrow(avg.exp)'.")
    }
    split.features <- as.factor(split.features)
    names(split.features) <- features
  }

  avg.scale.exp <- t(scale(t(avg.exp)))
  avg.exp <- cbind(
    features = features,
    as.data.frame(avg.exp, check.names = FALSE)
  )
  avg.pct <- cbind(
    features = features,
    as.data.frame(avg.pct, check.names = FALSE)
  )
  avg.scale.exp <- cbind(
    features = features,
    as.data.frame(avg.scale.exp, check.names = FALSE)
  )
  id.vars <- "features"
  variable.name <- "id"
  avg.exp <- melt(
    avg.exp,
    id.vars = id.vars,
    variable.name = variable.name,
    value.name = "avg.exp"
  )
  avg.pct <- melt(
    avg.pct,
    id.vars = id.vars,
    variable.name = variable.name,
    value.name = "avg.pct"
  )
  avg.scale.exp <- melt(
    avg.scale.exp,
    id.vars = id.vars,
    variable.name = variable.name,
    value.name = "avg.scale.exp"
  )
  df <- Reduce(dplyr::left_join, list(avg.exp, avg.scale.exp, avg.pct))
  df$features <- factor(df$features, features)

  if (length(group.by) > 0) {
    df$group.by <- group.by[df$id]
  }
  if (length(split.by) > 0) {
    df$split.by <- split.by[df$id]
  }
  if (length(split.features) > 0) {
    df$split.features <- split.features[df$features]
  }
  df
}




