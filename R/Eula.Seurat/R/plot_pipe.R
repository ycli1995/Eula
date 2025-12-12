#' @importFrom Eula.utils captureMsg getArgList getDefaultArgs ggsave2 pipeMsg
NULL

#' @importFrom SeuratObject Reductions
#' @export
pipe_dim_plot <- function(obj, params = list(), pipe.name = NULL, ...) {
  pipe.name <- c(pipe.name, "dim_plot")
  pipeMsg("Start")

  defaults <- list(
    outdir = getwd(),
    basic.size = 5.5,
    reductions = c("umap", "tsne"),
    corner.axis = TRUE,
    group.by = NULL,
    split.by = NULL,
    label = NULL,
    legend = NULL
  )
  params <- getDefaultArgs(params, defaults)
  print(str(params))
  list2env(params, envir = environment())

  reductions <- intersect(reductions, Reductions(obj))
  if (length(reductions) == 0) {
    stop("No valid reduction in the Seurat object.")
  }
  save_multi_dim_plot(
    obj,
    outdir = outdir,
    basic.size = basic.size,
    reductions = reductions,
    group.by = group.by,
    split.by = split.by,
    label = label,
    legend = legend,
    corner.axis = corner.axis,
    ...
  )
}

#' @export
pipe_dot_plot <- function(obj, params = list(), ...) {
  defaults <- list(
    outdir = getwd(),
    group.by = NULL,
    split.by = NULL,
    assay = NULL,
    theme = NULL,
    min.exp = 0,
    scale = TRUE,
    color.limits = c(-2.5, 2.5),
    size.limits = c(0, 100),
    coord.flip = FALSE
  )
  params <- getDefaultArgs(params, defaults)
  capture.msg(str(params))
  list2env(params, envir = environment())

  features <- getArgList(params[['features']])

  save_multi_dot_plot(
    obj = obj,
    features = features,
    outdir = outdir,
    assay = assay,
    group.by = group.by,
    split.by = split.by,
    split.features = split.features,
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

#' @export
pipe_violin_plot <- function(obj, params = list(), ...) {
  defaults <- list(
    outdir = getwd(),
    basic.size = 4,
    group.by = NULL,
    combine = TRUE,
    split.by = NULL,

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
    theme = NULL
  )
  params <- getDefaultArgs(params, defaults)
  capture.msg(str(params))
  list2env(params, envir = environment())

  features <- getArgList(params[['features']])

  save_violin_plots(
    obj = obj,
    features = features,
    outdir = outdir,
    basic.size = basic.size,
    group.by = group.by,
    combine = combine,
    split.by = split.by,

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

#' @importFrom Eula.utils normalizeFileName
#' @importFrom Seurat SplitObject
#' @export
save_dim_plots <- function(
    obj,
    outdir = getwd(),
    basic.size = 5.5,
    reductions = "umap",
    group.by = NULL,
    split.by = NULL,
    ...
) {
  split.by <- split.by[1]
  for (reduc in reductions) {
    outfile <- paste(reduc, paste(group.by, collapse = "."), "pdf", sep = ".")
    save_dim_plot(
      obj = obj,
      outfile = file.path(outdir, normalizeFileName(outfile)),
      basic.size = basic.size,
      reduction = reduc,
      group.by = group.by,
      ...
    )
  }
  if (length(split.by) == 0) {
    return(invisible(NULL))
  }
  obj.list <- SplitObject(obj, split.by = split.by)
  for (j in names(obj.list)) {
    for (reduc in reductions) {
      outfile <- paste(c(reduc, group.by, j, "pdf"), collapse = ".")
      save_dim_plot(
        obj = obj.list[[j]],
        outfile = file.path(outdir, normalizeFileName(outfile)),
        basic.size = basic.size,
        reduction = reduc,
        group.by = group.by,
        ...
      )
    }
  }
  invisible(NULL)
}

#' @importFrom Eula.utils normalizeFileName
#' @importFrom Seurat SplitObject
#' @export
save_feature_dim_plots <- function(
    obj,
    features,
    outdir = getwd(),
    basic.size = 5.5,
    reductions = "umap",
    combine = TRUE,
    split.by = NULL,
    order = NULL,
    ...
) {
  if (combine) {
    features <- list(combine = features)
  } else {
    features <- getFeaturesID(
      obj,
      features = features,
      keep.not.found = TRUE
    )
    features <- getFeaturesName(
      obj,
      features = features,
      col = "unique_name",
      keep.not.found = TRUE
    )
    names(features) <- features
  }
  split.by <- split.by[1]
  obj.list <- list()
  if (length(split.by) > 0) {
    obj.list <- SplitObject(obj, split.by = split.by)
  }
  for (reduc in reductions) {
    for (i in seq_along(features)) {
      fname <- names(features)[i]
      outfile <- paste(reduc, fname, "pdf", sep = ".")
      save_feature_dim_plot(
        obj = obj,
        features = features[[i]],
        outfile = file.path(outdir, normalizeFileName(outfile)),
        basic.size = basic.size,
        reduction = reduc,
        order = order,
        ...
      )
      for (j in seq_along(obj.list)) {
        nm <- gsub("\\/| ", "_", names(obj.list)[j])
        outfile <- paste(reduc, fname, nm, "pdf", sep = ".")
        save_feature_dim_plot(
          obj = obj.list[[j]],
          features = features[[i]],
          outfile = file.path(outdir, normalizeFileName(outfile)),
          basic.size = basic.size,
          reduction = reduc,
          order = order,
          ...
        )
      }
    }
  }
  invisible(NULL)
}

#' @importFrom Eula.utils normalizeFileName
#' @importFrom Seurat SplitObject
#' @export
save_violin_plots <- function(
    obj,
    features,
    outdir = getwd(),
    basic.size = 4,
    group.by = NULL,
    combine = TRUE,
    split.by = NULL,
    ...
) {
  group.by <- group.by %||% "Idents"
  if (combine) {
    features <- list(combine = features)
  } else {
    features <- getFeaturesID(
      obj,
      features = features,
      keep.not.found = TRUE
    )
    features <- getFeaturesName(
      obj,
      features = features,
      col = "unique_name",
      keep.not.found = TRUE
    )
    names(features) <- features
  }
  split.by <- split.by[1]
  obj.list <- list()
  if (length(split.by) > 0) {
    obj.list <- SplitObject(obj, split.by = split.by)
  }
  for (g in group.by) {
    for (i in seq_along(features)) {
      fname <- names(features)[i]
      outfile <- paste(fname, g, "pdf", sep = ".")
      save_violin_plot(
        obj = obj,
        features = features[[i]],
        outfile = file.path(outdir, normalizeFileName(outfile)),
        basic.size = basic.size,
        group.by = g,
        ...
      )
      for (j in seq_along(obj.list)) {
        nm <- gsub("\\/| ", "_", names(obj.list)[j])
        outfile <- paste(fname, g, nm, "pdf", sep = ".")
        save_violin_plot(
          obj = obj.list[[j]],
          features = features[[i]],
          outfile = file.path(outdir, normalizeFileName(outfile)),
          basic.size = basic.size,
          group.by = g,
          ...
        )
      }
    }
  }
  invisible(NULL)
}

#' @importFrom Eula.utils normalizeFileName
#' @export
save_dot_plots <- function(
    obj,
    features,
    outdir = getwd(),
    group.by = NULL,
    split.by = NULL,
    ...
) {
  group.by <- group.by %0% "Idents"
  for (g in group.by) {
    outfile <- paste("dot_plot", g, "pdf", sep = ".")
    save_dot_plot(
      obj = obj,
      features = features,
      outfile = file.path(outdir, normalizeFileName(outfile)),
      group.by = g,
      split.by = split.by,
      ...
    )
  }
  invisible(NULL)
}

#' @export
save_dim_plot <- function(
    obj,
    outfile = NULL,
    basic.size = 5.5,
    reduction = "umap",
    group.by = NULL,
    colors = NULL,
    label = NULL,
    legend = NULL,
    theme = NULL,
    corner.axis = FALSE,
    ...
) {
  group.by <- group.by %0% "Idents"
  if (identical(group.by, "Idents")) {
    obj$Idents <- Idents(obj)
  }
  obj@meta.data <- droplevels(obj@meta.data)

  for (i in group.by) {
    ## colors for each `group.by`
    obj <- checkColorMap(obj, i, obj@misc$colors[[i]])
  }

  legend <- legend %||% group.by
  legend <- intersect(legend, group.by)
  theme <- theme %||% theme_base_default()

  colors <- colors %||% obj@misc$colors[group.by]
  p <- dim_plot(
    obj,
    group.by = group.by,
    reduction = reduction,
    colors = colors,
    label = label,
    split.by = NULL,
    shape.by = NULL,
    combine = FALSE,
    label.repel = TRUE,
    legend = legend,
    corner.axis = corner.axis,
    theme = theme,
    ...
  )
  if (is.null(outfile)) {
    return(p)
  }
  w <- rep(basic.size, length(p))
  names(w) <- names(p)
  for (i in names(p)) {
    if (i %in% legend) {
      lgd.ncol <- ceiling(length(unique(p[[i]]$data[["colour"]])) / 20)
      lgd.size <- maxNChar(p[[i]]$data[["colour"]]) * 0.08 * lgd.ncol
      w[i] <- w[i] + lgd.size + basic.size / 10
    }
  }

  ## save plot
  p <- wrap_plots(p, nrow = 1)
  ggsave2(outfile, p, width = sum(w), height = basic.size, limitsize = FALSE)
  invisible(NULL)
}

#' @importFrom Eula.utils theme_no_legend
#' @importFrom ggplot2 labs
#' @importFrom patchwork wrap_plots
#' @export
save_feature_dim_plot <- function(
    obj,
    features,
    outfile = NULL,
    basic.size = 5.5,
    reduction = NULL,
    dims = c(1, 2),
    assay = NULL,
    slot = "data",
    ncol = NULL,
    order = NULL,
    shuffle = FALSE,
    seed = 42,
    colors = NULL,
    pt.size = NULL,
    pt.alpha = NULL,
    raster = NULL,
    corner.axis = FALSE,
    theme = NULL,
    ...
) {
  features <- getFeaturesID(obj, features, keep.not.found = TRUE)
  p <- feature_dim_plot(
    object = obj,
    features = features,
    reduction = reduction,
    dims = dims,
    assay = assay,
    slot = slot,
    split.by = NULL,
    ncol = NULL,
    order = order,
    shuffle = shuffle,
    seed = seed,
    combine = FALSE,
    colors = colors,
    pt.size = pt.size,
    pt.alpha = pt.alpha,
    raster = raster,
    corner.axis = corner.axis,
    theme = theme,
    ...
  )
  names(p) <- getFeaturesName(obj, features, "unique_name")
  for (i in seq_along(p)) {
    p[[i]] <- p[[i]] + labs(title = names(p)[i])
  }
  if (is.null(outfile)) {
    return(p)
  }
  ncol <- ncol %||% min(ceiling(sqrt(length(features))), 5)
  if (length(p) == 1) {
    ncol <- 1
    p <- p[[1]]
  } else {
    p <- wrap_plots(p, ncol = ncol)
  }
  nrow <- ceiling(length(features) / ncol)
  w <- basic.size * (6/5) * ncol
  h <- basic.size * nrow
  ggsave2(outfile, p, width = w, height = h, limitsize = FALSE)
  invisible(NULL)
}

#' @importFrom ggplot2 labs
#' @importFrom patchwork wrap_plots
#' @export
save_violin_plot <- function(
    obj,
    features,
    outfile = NULL,
    basic.size = 4,
    assay = NULL,
    slot = "data",
    ncol = NULL,
    group.by = NULL,
    cells = NULL,
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
  group.by <- group.by[1] %0% "Idents"
  if (identical(group.by, "Idents")) {
    obj$Idents <- Idents(obj)
  }
  obj@meta.data <- droplevels(obj@meta.data)
  obj <- checkColorMap(obj, group.by, obj@misc$colors[[group.by]])
  colors <- colors %||% obj@misc$colors[[group.by]]

  features <- getFeaturesID(obj, features, keep.not.found = TRUE)
  p <- violin_plot(
    object = obj,
    features = features,
    assay = assay,
    slot = slot,
    group.by = group.by,
    split.by = NULL,
    ncol = NULL,
    combine = FALSE,
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
  names(p) <- getFeaturesName(
    obj,
    features,
    col = "unique_name",
    keep.not.found = TRUE
  )
  for (i in seq_along(p)) {
    p[[i]] <- p[[i]] + labs(title = names(p)[i]) + theme_no_legend()
  }
  if (is.null(outfile)) {
    return(p)
  }
  if (coord.flip) {
    lab.size <- maxNChar(p[[1]]$data[["y"]]) * 0.12
    w <- basic.size + lab.size + basic.size / 10
    h <- length(unique(p[[1]]$data[["y"]])) * 0.35 + basic.size / 5
  } else {
    lab.size <- maxNChar(p[[1]]$data[["x"]]) * 0.12 * 2/3
    h <- basic.size + lab.size + basic.size / 10
    w <- length(unique(p[[1]]$data[["x"]])) * 0.35 + basic.size / 5
  }
  w <- max(3.5, w)
  h <- max(4, h)
  ncol <- ncol %||% min(ceiling(sqrt(length(features))), 5)
  nrow <- ceiling(length(features) / ncol)
  if (length(p) == 1) {
    ncol <- 1
    p <- p[[1]]
  } else {
    p <- wrap_plots(p, ncol = ncol)
  }
  message("Basic size: w = ", w, ", h = ", h)
  w <- w * ncol
  h <- h * nrow
  ggsave2(outfile, p, width = w, height = h, limitsize = FALSE)
  invisible(NULL)
}

#' @export
save_dot_plot <- function(
    obj,
    features,
    outfile = NULL,
    assay = NULL,
    group.by = NULL,
    split.by = NULL,
    split.features = NULL,
    min.exp = 0,
    scale = TRUE,
    colors = NULL,
    color.limits = c(-2.5, 2.5),
    size.limits = c(0, 100),
    coord.flip = FALSE,
    theme = NULL,
    ...
) {
  group.by <- group.by %0% "Idents"
  if (identical(group.by, "Idents")) {
    obj$Idents <- Idents(obj)
  }
  obj@meta.data <- droplevels(obj@meta.data)

  features <- getFeaturesID(obj, features)
  p <- dot_plot(
    object = obj,
    features = features,
    assay = assay,
    group.by = group.by,
    split.by = split.by,
    min.exp = min.exp,
    scale = scale,
    split.features = split.features,
    colors = colors,
    color.limits = color.limits,
    size.limits = size.limits,
    coord.flip = coord.flip,
    theme = theme,
    ...
  )
  features.col <- "y"
  if (coord.flip) {
    features.col <- "x"
  }
  levels(p$data[[features.col]]) <- getFeaturesName(
    object = obj,
    features = levels(p$data[[features.col]]),
    col = "unique_name",
    uniq = FALSE
  )
  if (is.null(outfile)) {
    return(p)
  }
  w.str <- maxNChar(p$data$y) * 0.08
  h.str <- maxNChar(p$data$x) * 0.075
  w <- nlevels(p$data$x) * 0.25 + w.str + 1
  h <- nlevels(p$data$y) * 0.25 + h.str + 0.5
  ggsave2(outfile, p, width = w, height = h, limitsize = FALSE)
  invisible(NULL)
}
