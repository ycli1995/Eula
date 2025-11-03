
#' @export dim_plot
dim_plot <- function(object, ...) {
  UseMethod("dim_plot", object)
}

#' @importFrom patchwork wrap_plots
#' @importFrom ggplot2 theme_classic
#'
#' @export
#' @method dim_plot data.frame
dim_plot.data.frame <- function(
    object,
    group.by,
    split.by = NULL,
    shape.by = NULL,

    dims = c(1, 2),
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

    legend = NULL,

    colors = NULL,
    pt.size = NULL,
    pt.alpha = NULL,
    raster = NULL,
    raster.dpi = c(512, 512),
    facet.args = list(),
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

  theme <- theme %||% theme_classic()

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
      alpha = pt.alpha,
      ...
    )
    if (i %in% legend) {
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
    dims = c(1, 2),
    reduction = NULL,
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
    theme = theme,
    ...
  )
}

#' @importFrom ggplot2 element_blank element_text theme
#' @export
theme_facet <- function(...) {
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 12),
    validate = TRUE,
    ...
  )
}

#' @export
theme_no_legend <- function(...) {
  theme(legend.position = "none", validate = TRUE, ...)
}
