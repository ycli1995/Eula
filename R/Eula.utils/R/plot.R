
#' @importFrom ggplot2 geom_jitter geom_violin position_jitterdodge
#' scale_fill_gradientn scale_fill_manual
#' @export
single_violin_plot <- function(
    data,
    colors = NULL,
    fills = NULL,
    fill.limits = NULL,
    noise.scale = 1e+5,
    scale = "width",
    adjust = 1,
    pt.size = 0,
    pt.alpha = 1,
    raster = NULL,
    raster.dpi = c(512, 512),
    facet.args = list(),
    labs.args = list(),
    theme = NULL,
    ...
) {
  cols.required <- c("x", "y")
  cols.optional <- c("alpha")
  data <- normalize_ggplot_data(data, cols = cols.required)
  mapping <- get_mapping_from_data(data, cols = c(cols.required, cols.optional))

  if (noise.scale > 0) {
    noise <- rnorm(nrow(data)) / noise.scale
    data[['y']] <- data[['y']]
  }
  p <- ggplot(data) +
    geom_violin(
      mapping = mapping,
      scale = scale,
      adjust = adjust,
      ...
    )
  if (pt.size > 0) {
    geom.jitter <- geom_jitter(
      position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.9),
      size = pt.size,
      alpha = pt.alpha,
      show.legend = FALSE
    )
    if (isTRUE(raster)) {
      checkPackages("ggrastr")
      geom.jitter <- ggrastr::rasterize(geom.jitter, dpi = raster.dpi)
    }
    p <- p + geom.jitter
  }
  if ("colour" %in% colnames(data) & !is.null(colors)) {
    if (is.numeric(data$colour)) {
      stop("'colour' must be factor or character for `single_violin_plot`.")
    }
    p <- p + scale_color_manual(values = colors)
  }
  if ("fill" %in% colnames(data) & !is.null(fills)) {
    scale.fill <- if (is.numeric(data$fill)) {
      scale_fill_gradientn(colors = colors, limits = fill.limits)
    } else {
      scale_fill_manual(values = colors)
    }
    p <- p + scale.fill
  }
  p <- .add_my_facet_split(p = p, args = facet.args)
  p <- add_my_labs(p, args = labs.args)
  p <- add_my_theme(p, args = theme)
  p
}

#' @importFrom ggplot2 guide_colorbar guide_legend guides scale_size
#' @export
single_dot_plot <- function(
    data,
    colors = NULL,
    color.limits = NULL,
    size.limits = NULL,
    facet.args = list(),
    labs.args = list(),
    theme = NULL,
    ...
) {
  colors <- colors %||% c(
    "#FFF7EC", "#FEE8C8", "#FDD49E", "#FDBB84", "#FC8D59",
    "#EF6548", "#D7301F", "#B30000", "#7F0000"
  )
  cols.required <- c("x", "y", "colour", "size")
  cols.optional <- c("alpha")
  data <- normalize_ggplot_data(data, cols = cols.required)
  mapping <- get_mapping_from_data(data, cols = c(cols.required, cols.optional))
  p <- ggplot(data) +
    geom_point(mapping = mapping, ...) +
    scale_color_gradientn(colors = colors, limits = color.limits) +
    scale_size(limits = size.limits) +
    guides(
      size = guide_legend(order = 1),
      colour = guide_colorbar(order = 2)
    )
  p <- .add_my_facet_split(p = p, args = facet.args)
  p <- add_my_labs(p, args = labs.args)
  p <- add_my_theme(p, args = theme)
  p
}

#' @importFrom ggplot2 facet_wrap geom_point ggplot scale_color_gradientn
#' scale_color_manual
#' @export
single_dim_plot <- function(
    data,
    colors = NULL,
    color.limits = NULL,
    label.data = NULL,
    label.repel = TRUE,
    label.size = 4,
    label.color = "black",
    label.box = FALSE,
    pt.size = NULL,
    raster = NULL,
    raster.dpi = c(512, 512),
    facet.args = list(),
    labs.args = list(),
    corner.axis = FALSE,
    theme = NULL,
    guides = NULL,
    ...
) {
  if (nrow(data) > 1e+5 & is.null(raster)) {
    message("Rasterizing points since number of points > 100,000.")
    message("To disable rasterizing, set `raster = FALSE`")
  }
  raster <- raster %||% nrow(data) > 1e+5
  if (!is.null(raster.dpi)) {
    if (!is.numeric(raster.dpi) || length(raster.dpi) != 2) {
      stop("'raster.dpi' must be a two-length numeric vector.")
    }
  }
  if (!requireNamespace("scattermore", quietly = TRUE)) {
    fastWarning("Package 'scattermore' is not installed. Set `raster = FALSE`")
    raster <- FALSE
  }

  pt.size <- pt.size %||% auto_point_size(data, raster)

  cols.required <- c("x", "y")
  cols.optional <- c("colour", "size", "alpha")
  data <- normalize_ggplot_data(data, cols = cols.required)
  data[['size']] <- NULL
  mapping <- get_mapping_from_data(data, cols = c(cols.required, cols.optional))

  p <- ggplot(data)
  if (isTRUE(raster)) {
    p <- p +
      scattermore::geom_scattermore(
        mapping = mapping,
        pixels = raster.dpi,
        pointsize = pt.size,
        ...
      )
  } else {
    p <- p + geom_point(mapping = mapping, size = pt.size, ...)
  }
  if ("colour" %in% colnames(data)) {
    if (is.numeric(data[['colour']])) {
      if (!is.null(colors)) {
        p <- p + scale_color_gradientn(colors = colors, limits = color.limits)
      }
    } else {
      if (!is.null(colors)) {
        if (length(colors) < length(unique(data$colour))) {
          fastWarning("Insufficient colors for manual. Use default colors")
        } else {
          p <- p + scale_color_manual(values = colors)
        }
      }
      p <- p + dim_plot_guides()
    }
  }
  if (!is.null(label.data)) {
    p <- p + add_labels(
      data = label.data,
      repel = label.repel,
      box = label.box,
      size = label.size,
      color = label.color
    )
  }
  p <- .add_my_facet_split(p = p, args = facet.args)
  p <- add_my_labs(p, args = labs.args)
  p <- add_my_theme(p, args = theme)
  if (corner.axis) {
    p <- add_corner_axis(p)
  }
  p
}

#' @importFrom dplyr across
#' @export
get_median_position <- function(data, group.by) {
  cols <- c("x", "y")
  checkColumns(data, cols)
  group.by <- intersect(group.by, colnames(data))
  if (length(group.by) == 0) {
    stop("No 'group.by' input")
  }
  label.data <- data %>%
    group_by(across(group.by)) %>%
    summarize(x = median(x), y = median(y)) %>%
    as.data.frame()
  label.data
}

AES.USE <- list(
  point = c("x", "y", "colour", "shape", "size", "alpha"),
  dotplot = c("x", "y", "colour", "size", "alpha")
)

#' @export
normalize_ggplot_data <- function(data, cols = NULL, ...) {
  if (!is.data.frame(data)) {
    data <- as.data.frame(data)
  }
  cols <- unique(cols)
  checkColumns(data, cols)
  if ("split" %in% colnames(data)) {
    data[["split"]] <- as.factor(data[["split"]])
  }
  data
}

#' @importFrom ggplot2 aes
#' @export
get_mapping_from_data <- function(data, cols) {
  if ("color" %in% colnames(data) & !"colour" %in% colnames(data)) {
    data[['colour']] <- data[["color"]]
    data[['color']] <- NULL
  }
  mapping <- aes()
  data <- data[, intersect(colnames(data), cols)]
  for (i in colnames(data)) {
    mapping[[i]] <- as.symbol(i)
  }
  mapping
}

#' @importFrom ggplot2 geom_label geom_text
#' @importFrom ggrepel geom_label_repel geom_text_repel
#' @export
add_labels <- function(data, repel = TRUE, box = FALSE, ...) {
  cols <- c("x", "y", "label")
  data <- normalize_ggplot_data(data, cols = cols)
  geom.use <- if (box) {
    ifelse(repel, geom_label_repel, geom_label)
  } else {
    ifelse(repel, geom_text_repel, geom_text)
  }
  mapping <- get_mapping_from_data(data, cols)
  geom.use(data = data, mapping = mapping, show.legend = FALSE, ...)
}

.add_my_facet_split <- function(p, args = list(), ...) {
  if (!"split" %in% colnames(p$data)) {
    return(p)
  }
  add_my_facet(p = p, facets = ~split, args = args, ...)
}

#' @export
get_label_data <- function(data, group.by = NULL, ...) {
  group.by <- group.by %||% "colour"
  group.by <- c(group.by, "split")
  group.by <- intersect(group.by, colnames(data))
  if (length(group.by) == 0) {
    stop("No valid 'group.by' for summarizing the label positions.")
  }
  label.data <- data %>%
    dplyr::group_by(across(all_of(group.by))) %>%
    dplyr::summarize(x = median(x), y = median(y))
  label.data$label <- label.data[[group.by[1]]]
  label.data
}

#' @importFrom ggplot2 facet_wrap
add_my_facet <- function(p, facets, args = list(), ...) {
  fw <- facet_wrap(facets = facets)
  use.args <- intersect(names(args), names(fw$params))
  for (i in use.args) {
    fw$params[[i]] <- args[[i]]
  }
  p + fw
}

#' @importFrom ggplot2 labs
#' @export
add_my_labs <- function(p, args = list(), ...) {
  if (length(args) == 0) {
    return(p)
  }
  my.labs <- labs()
  for (i in names(args)) {
    my.labs[[i]] <- args[[i]]
  }
  p + my.labs
}

#' @importFrom ggplot2 theme
#' @export
add_my_theme <- function(p, args = list(), ...) {
  if (length(args) == 0) {
    return(p)
  }
  if (is(args, "ggplot2::theme")) {
    return(p + args)
  }
  my.theme <- theme()
  for (i in seq_along(args)) {
    my.theme[[i]] <- args[[i]]
  }
  p + my.theme
}

#' @export
auto_point_size <- function(data, raster = NULL) {
  ifelse(isTRUE(raster), 1, min(1, 1583 / nrow(data)))
}

#' @importFrom ggplot2 guides guide_legend
#' @export
dim_plot_guides <- function(...) {
  args <- list(...)
  if ("colour" %in% names(args) | "color" %in% names(args)) {
    return(guides(...))
  }
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)), ...)
}

