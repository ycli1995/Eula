
#' @export checkColorMap
checkColorMap <- function(x, ...) {
  UseMethod("checkColorMap", x)
}

#' @export
#' @method checkColorMap default
checkColorMap.default <- function(x, colors = NULL, ...) {
  if (!is.factor(x)) {
    x <- factor(x)
  }
  all_levels <- levels(x)
  if (length(colors) == 0) {
    fastWarning(
      "'colors' is NULL. Calling 'setColors' for the following levels: \n  ",
      paste(all_levels, collapse = ", ")
    )
    return(SetColor(x, ...))
  }
  if (length(names(colors)) == 0) {
    color_names <- all_levels[1:min(length(all_levels), length(colors))]
    fastWarning(
      "No name in 'colors', set color names with:\n  ",
      paste(color_names, collapse = ", ")
    )
    names(colors) <- color_names
  }
  if (length(colors) < length(all_levels)) {
    miss_levels <- setdiff(all_levels, names(colors))
    fastWarning(
      "Calling 'setColors' since the following levels mapped to no color: \n  ",
      paste(miss_levels, collapse = ", ")
    )
    colors0 <- SetColor(x, ...)
    colors0[names(colors)] <- colors
    return(colors0)
  }
  if (!all(all_levels %in% names(colors))) {
    miss_levels <- setdiff(all_levels, names(colors))
    fastWarning(
      "Calling 'SetColor' since the following levels mapped to no color: \n  ",
      paste(miss_levels, collapse = ", ")
    )
    colors0 <- SetColor(x, ...)
    cmm.nm <- intersect(names(colors0), names(colors))
    colors0[cmm.nm] <- colors[cmm.nm]
    return(colors0)
  }
  colors <- colors[all_levels]
  return(colors)
}

#' @export
#' @method checkColorMap data.frame
checkColorMap.data.frame <- function(x, column, colors = NULL, ...) {
  column <- column[1]
  checkColumns(x, column)
  checkColorMap(x[, column], colors = colors, ...)
}

#' @export
setColors <- function(x, type = "tsne", tag = "set1", ...) {
  if (!is.factor(x)) {
    x <- as.factor(x)
  }
  color <- fetch_color(n = nlevels(x), type = type, tag = tag, ...)
  names(color) <- levels(x)
  color
}

#' @export
load_color_list <- function() {
  dir <- system.file("extdata", package = .packageName)
  readYAML(file.path(dir, "Colors.yml"))
}

#' @export
fetch_color <- function(n = 0, type = NULL, tag = NULL, verbose = FALSE) {
  if (!is.numeric(n)) {
    stop("'n' must be numeric.")
  }
  color.list <- load_color_list()
  all.types <- c(names(color.list), "random")
  type <- type %||% all.types[1]
  type <- match.arg(type, all.types)

  if (type == "random") {
    color.use <- fetch_random_color(n = n, usepalette = TRUE)
    verboseMsg("color : ", type, "->", tag, "->", n)
    return(color.use)
  }

  tag  <- match.arg(tag, names(color.list[[type]]))
  n.available <- as.numeric(names(color.list[[type]][[tag]]))
  n.select <- n.available[n.available >= n]
  if (length(n.select) > 0) {
    if (n == min(n.select)) {
      n.select <- as.character(min(n.select))
      color.use <- color.list[[type]][[tag]][[n.select]]
    } else {
      n.select <- as.character(min(n.select))
      color.use <- head(color.list[[type]][[tag]][[n.select]], n)
    }
    verboseMsg("color : ", type, "->", tag, "->", n.select)
    return(color.use)
  }
  if (length(n.available[n.available == 0]) > 0) {
    n.select <- as.character(0)
  } else {
    n.select <- as.character(max(n.available))
  }
  verboseMsg("color : ", type, "->", tag, "->", n.select)
  color.use <- color.list[[type]][[tag]][[n.select]]
  if (length(color.use) >= n) {
    return(head(color.use, n))
  }
  colorRampPalette(color.use)(n)
}

#' @export
fetch_random_color <- function(
    n = 1,
    usepalette = FALSE,
    hue = " ",
    luminosity = " ",
    ...
) {
  checkPackages("randomcoloR")
  if (usepalette) {
    set.seed(1)
    return(randomcoloR::distinctColorPalette(k = n), ...)
  }
  randomcoloR::randomColor(count = n, hue = hue, luminosity = luminosity)
}

# show_all_color <- function() {
#   library(dplyr, warn.conflicts = F)
#   dt <- reshape2::melt(color.list)
#   dt <- dt %>% group_by(L1,L2,L3) %>% mutate(x = 1:n()) %>% arrange(L1, L2, as.numeric(L3), x)
#   dt$L3 <- factor(dt$L3, levels = as.character(0:max(as.numeric(dt$L3))))
#
#   cc <- unique(as.character(dt$value))
#   names(cc) <- cc
#
#   library(ggplot2, warn.conflicts = F)
#   p <- list()
#   for ( i in unique(dt$L1) ) {
#     p[[i]] <-  ggplot(dt %>% filter(L1 == i), aes(x = x, y = L3)) + geom_tile(aes(fill = value), color = "grey") +
#       facet_grid(L2 ~ ., scales = "free_y", switch = "y") +
#       scale_fill_manual(values = cc) +
#       scale_x_continuous(expand = expand_scale()) +
#       ylab(i) +
#       theme_minimal() +
#       theme(legend.position = "none", strip.placement = "outside", panel.grid = element_blank(),
#             axis.title.x = element_blank(), axis.text.x = element_blank())
#   }
#   cowplot::plot_grid(plotlist = p, ncol = 1, axis = "ltrb", align = "hv")
# }

