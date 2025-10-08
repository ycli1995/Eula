

single_dim_plot <- function(
    data,
    colors = NULL,
    pt.size = NULL,
    raster = NULL, 
    raster.dpi = c(512, 512),
    labs.args = list(),
    theme = NULL,
    ...
) {
  if (nrow(data) > 1e+5 & is.null(raster)) {
    message("Rasterizing points since number of points > 100,000.")
    message("To disable rasterizing, set `raster = FALSE`")
  }
  raster <- raster %||% nrow(data) > 1e+5
  pt.size <- pt.size %||% auto_point_size(data, raster)
  if (!is.null(raster.dpi)) {
    if (!is.numeric(raster.dpi) || length(raster.dpi) != 2) {
      stop("'raster.dpi' must be a two-length numeric vector.")
    }
  }
  if (!is.data.frame(data)) {
    data <- as.data.frame(data)
  }
  if (!all(c("x", "y") %in% colnames(data))) {
    stop("'x' and 'y' must be included in columns of 'data'.")
  }
  data[['size']] <- NULL
  p <- ggplot(data)
  mapping <- get_mapping_from_data(data, type = "point")
  if (isTRUE(raster)) {
    p <- p + 
      scattermore::geom_scattermore(
        mapping = mapping,
        pixels = raster.dpi,
        pointsize = pt.size
      )
  } else {
    p <- p + 
      geom_point(mapping = mapping, size = pt.size)
  }
  if ("color" %in% colnames(data) & !is.null(colors)) {
    if (is.numeric(data$color)) {
      p <- p + scale_color_gradientn(values = colors)
    } else {
      p <- p + scale_color_manual(values = colors)
    }
  }
  p <- add_my_labs(p, args = labs.args)
  p <- add_my_theme(p, args = theme)
  p
}

auto_point_size <- function(data, raster = NULL) {
  ifelse(isTRUE(raster), 1, min(1, 1583 / nrow(data)))
}

aes.use <- list(
  point = c("x", "y", "color", "shape", "size", "alpha")
)

get_mapping_from_data <- function(data, type = "point") {
  mapping <- aes()
  data <- data[, intersect(colnames(data), aes.use[[type]])]
  for (i in colnames(data)) {
    mapping[[i]] <- as.symbol(i)
  }
  mapping
}

add_my_labs <- function(p, args = list(), ...) {
  if (length(args) == 0) {
    return(p)
  }
  my.labs <- labs()
  for (i in seq_along(args)) {
    my.labs[[i]] <- args[[i]]
  }
  p + my.labs
}

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

