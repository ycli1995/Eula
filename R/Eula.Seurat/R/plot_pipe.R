#' @importFrom Eula.utils ggsave2
NULL

#' @importFrom SeuratObject Reductions
#' @export
pipe_DimPlot <- function(obj, params = list(), ...) {
  params[['outdir']] <- params[['outdir']] %||% getwd()
  params[['basic.size']] <- params[['basic.size']] %||% 5.5
  params[['reductions']] <- params[['reductions']] %||% c("umap", "tsne")
  params[['corner.axis']] <- params[['corner.axis']] %||% TRUE

  print(str(params))
  list2env(params, envir = environment())

  reductions <- intersect(reductions, Reductions(obj))
  if (length(reductions) == 0) {
    stop("No valid reduction in the Seurat object.")
  }

  group.by <- params[['group.by']]
  split.by <- params[['split.by']]
  label <- params[['label']]
  legend <- params[['legend']]

  save_multi_DimPlot(
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

#' @importFrom Eula.utils norm_fname
#' @importFrom Seurat SplitObject
#' @export
save_multi_DimPlot <- function(
    obj,
    outdir = getwd(),
    basic.size = 6,
    reductions = "umap",
    group.by = NULL,
    split.by = NULL,
    ...
) {
  split.by <- split.by[1]
  for (reduc in reductions) {
    outfile <- paste(reduc, paste(group.by, collapse = "."), "pdf", sep = ".")
    save_DimPlot(
      obj = obj,
      outfile = file.path(outdir, norm_fname(outfile)),
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
      save_DimPlot(
        obj = obj.list[[j]],
        outfile = file.path(outdir, norm_fname(outfile)),
        basic.size = basic.size,
        reduction = reduc,
        group.by = group.by,
        ...
      )
    }
  }
  invisible(NULL)
}

#' @importFrom ggplot2 ggsave
#' @importFrom tidydr theme_dr
#' @export
save_DimPlot <- function(
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

#' @export
pipe_DotPlot <- function(obj, params = list(), ...) {
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

  features <- norm_list_param(params[['features']])

  save_multi_DotPlot(
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

#' @importFrom Eula.utils norm_fname
#' @export
save_multi_DotPlot <- function(
    obj,
    features,
    outdir = getwd(),
    group.by = NULL,
    split.by = NULL,
    ...
) {
  group.by <- group.by %0% "Idents"
  for (g in group.by) {
    outfile <- paste("DotPlot", g, "pdf", sep = ".")
    save_DotPlot(
      obj = obj,
      features = features,
      outfile = file.path(outdir, norm_fname(outfile)),
      group.by = g,
      split.by = split.by,
      ...
    )
  }
  invisible(NULL)
}

#' @export
save_DotPlot <- function(
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
