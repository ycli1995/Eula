#' @importFrom Eula.utils ggsave2
NULL

#' @importFrom SeuratObject Reductions
#' @export
pipe_DimPlot <- function(obj, params = list(), ...) {
  params[['outdir']] <- params[['outdir']] %||% getwd()
  params[['basic.size']] <- params[['basic.size']] %||% 6
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
    basic.size = 6,
    reduction = "umap",
    group.by = NULL,
    label = NULL,
    legend = NULL,
    theme = NULL,
    corner.axis = FALSE,
    ...
) {
  if (length(group.by) == 0) {
    obj$Idents <- Idents(obj)
    group.by <- "Idents"
    obj@misc$colors$Idents <- setColors(obj$Idents)
  }
  obj@meta.data <- droplevels(obj@meta.data)

  for (i in group.by) {
    ## colors for each `group.by`
    obj <- checkColorMap(obj, i, obj@misc$colors[[i]])
  }

  legend <- legend %||% group.by
  legend <- intersect(legend, group.by)
  theme <- theme %||% dot_theme_default()
  if (corner.axis) {
    theme <- theme + theme_dr()
  }

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
  ggsave(outfile, p, width = sum(w), height = basic.size, limitsize = FALSE)
  invisible(NULL)
}
