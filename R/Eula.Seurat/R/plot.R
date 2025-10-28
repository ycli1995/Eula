
#' @export dim_plot
dim_plot <- function(object, ...) {
  UseMethod("dim_plot", object)
}

#' @importFrom patchwork wrap_plots
#' @importFrom ggplot2 theme_classic
#'
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
  colnames(data)[1:2] <- c("x", "y")

  if (!is.null(shape.by)) {
    shape.by <- shape.by[1]
    data[, "shape"] <- object[[shape.by, drop = TRUE]]
  }
  if (!is.null(split.by)) {
    split.by <- split.by[1]
    split <- FetchData(object, vars = split.by, clean = TRUE)[split.by]
    data <- data[rownames(split), ]
    data[, "split"] <- split
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

  if (isTRUE(shuffle)) {
    set.seed(seed)
    data <- data[sample(seq_len(nrow(data))), , drop = FALSE]
  }

  if (is.character(colors)) {
    colors <- sapply(group.by, function(x) colors, simplify = FALSE)
  }

  key <- key %||% Key(reduc.obj)
  labs.args <- setNames(paste0(key, dims), c("x", "y"))
  labs.args <- as.list(labs.args)

  ncol <- ncol %||% 4
  if (length(group.by) > 1 & "split" %in% colnames(data)) {
    facet.args[['nrow']] <- 1
    facet.args[['ncol']] <- NULL
  }

  theme <- theme %||% theme_classic()
  if ("split" %in% colnames(data)) {
    theme <- theme + theme_facet()
  }

  pt.alpha <- pt.alpha %||% 1

  plist <- list()
  for (i in group.by) {
    data$colour <- as.factor(data[[i]])
    cols <- colors[[i]]

    label.data <- NULL
    if (i %in% label) {
      label.data <- get_label_data(data)
    }
    labs.args[["colour"]] <- i
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
  }
  if ("split" %in% colnames(data)) {
    ncol <- 1
  }
  if (combine) {
    ncol <- min(ncol, length(plist))
    plist <- wrap_plots(plist, ncol = ncol)
  }
  plist
}

#' @importFrom ggplot2 element_blank element_text theme
#' @export
theme_facet <- function(...) {
  theme(
    strip.background = element_blank(),
    strip.text = element_text(),
    validate = TRUE,
    ...
  )
}
