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
  params <- getDefaultArgs(defaults, params)
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
  params <- fetch_default_params(defaults, params)
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
    features <- getFeaturesID(obj, features)
    features <- getFeaturesName(obj, features, "unique_name")
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
        combine = combine,
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
          combine = combine,
          order = order,
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

  p <- dim_plot(
    obj,
    group.by = group.by,
    reduction = reduction,
    colors = obj@misc$colors,
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
    combine = TRUE,
    colors = NULL,
    pt.size = NULL,
    pt.alpha = NULL,
    raster = NULL,
    corner.axis = FALSE,
    theme = NULL,
    ...
) {
  features <- getFeaturesID(obj, features)
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
