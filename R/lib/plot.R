

single_dim_plot <- function(
    data,
    colors = NULL,
    label.data = NULL,
    label.repel = TRUE,
    label.size = 4,
    label.color = "black",
    label.box = FALSE,
    pt.size = NULL,
    raster = NULL, 
    raster.dpi = c(512, 512),
    nrow = NULL,
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
  .check_columns(data, c("x", "y"))
  if ("split" %in% colnames(data)) {
    data[["split"]] <- as.factor(data[["split"]])
  }
  data[['size']] <- NULL
  p <- ggplot(data)
  mapping <- get_mapping_from_data(data, aes.use[["point"]])
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
  if (!is.null(label.data)) {
    p <- p + add_labels(
      data = label.data, 
      repel = label.repel, 
      box = label.box, 
      size = label.size, 
      color = label.color
    )
  }
  if ("split" %in% colnames(data)) {
    p <- p + facet_wrap(~split, nrow = nrow, scales = "fixed")
  }
  p <- add_my_labs(p, args = labs.args)
  p <- add_my_theme(p, args = theme)
  p
}

auto_point_size <- function(data, raster = NULL) {
  ifelse(isTRUE(raster), 1, min(1, 1583 / nrow(data)))
}

aes.use <- list(
  point = c("x", "y", "color", "shape", "size", "alpha", "split")
)

get_mapping_from_data <- function(data, cols) {
  mapping <- aes()
  data <- data[, intersect(colnames(data), cols)]
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

get_median_position <- function(data, group.by) {
  cols <- c("x", "y")
  .check_columns(data, cols)
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

add_labels <- function(data, repel = TRUE, box = FALSE, ...) {
  cols <- c("x", "y", "label")
  .check_columns(data, cols)
  geom.use <- if (box) {
    ifelse(repel, geom_label_repel, geom_label)
  } else {
    ifelse(repel, geom_text_repel, geom_text)
  }
  mapping <- get_mapping_from_data(data, cols)
  geom.use(data = data, mapping = mapping, show.legend = FALSE, ...)
}

.check_columns <- function(data, cols) {
  if (!all(cols %in% colnames(data))) {
    cols <- paste(cols, collapse = ", ")
    stop("Missing columns for the input 'data': \n  ", cols)
  }
  invisible(NULL)
}

facet_theme <- function(...) {
  theme(
    strip.background = element_blank(), 
    strip.text = element_text(face = "plain"), 
    validate = TRUE, 
    ...
  )
}
