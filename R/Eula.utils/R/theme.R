
#' @importFrom ggplot2 element_blank element_line element_text margin theme
#' theme_bw
#' @importFrom grid unit
#' @export
theme_base_default <- function(family = "Arial", ...) {
  theme_bw() +
    theme(
      plot.title = element_text(
        color = "#000000",
        size = 16,
        face = "plain",
        hjust = 0.5
      ),
      plot.margin = unit(c(10, 10, 10, 10), "mm"),

      axis.text = element_text(color = "#000000", size = 11),
      axis.title = element_text(color = "#000000", size = 14, face = "plain"),
      axis.title.x = element_text(margin = margin(2.5, 0, 2.5, 0, "mm")),
      axis.title.y = element_text(margin = margin(0, 2.5, 0, 2.5, "mm")),
      axis.ticks = element_line(color = "#000000", linewidth = 0.5),
      axis.ticks.length = unit(0.1, 'cm'),

      axis.text.y = element_text(hjust = 1, vjust = 0.5),

      legend.title = element_blank(),
      legend.text = element_text(size = 11),
      ...
    ) +
    theme(text = element_text(family = family))
}

#' @importFrom ggplot2 element_blank element_rect theme
#' @export
theme_border_default <- function(...) {
  theme(
    panel.border = element_rect(color = "#000000", linewidth = 0.8, fill = NA),
    axis.line = element_blank(),
    ...
  )
}

#' @importFrom ggplot2 theme
#' @export
theme_no_axis <- function(...) {
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),

    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),

    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    ...
  )
}

#' @importFrom ggplot2 theme
#' @export
theme_no_legend <- function(...) {
  theme(legend.position = "none", ...)
}

#' @importFrom ggplot2 element_blank theme
#' @export
theme_no_grid <- function(...) {
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
}

#' @importFrom ggplot2 element_blank element_text theme
#' @export
theme_facet <- function(...) {
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 12),
    ...
  )
}

#' @importFrom ggplot2 theme
#' @export
theme_rotate_axis <- function(...) {
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    ...
  )
}

#' @export
theme_dim_default <- function(family = "Arial", ...) {
  theme_base_default(family = family, ...)
}

#' @export
theme_dot_default <- function(family = "Arial", ...) {
  theme_base_default(family = family) +
    theme_no_grid() +
    theme_border_default() +
    theme_rotate_axis(...)
}

#' @importFrom ggplot2 get_last_plot ggsave
#' @importFrom grDevices cairo_pdf
#' @importFrom tools file_ext
#' @export
ggsave2 <- function(filename, plot = get_last_plot(), device = NULL, ...) {
  if (requireNamespace("extrafont", quietly = TRUE)) {
    suppressMessages(require("extrafont", quietly = TRUE))
    suppressMessages(require("extrafontdb", quietly = TRUE))
    suppressMessages(require("Rttf2pt1", quietly = TRUE))
  }
  ext <- tools::file_ext(filename)
  if (tolower(ext) == "pdf") {
    device <- device %||% cairo_pdf
  }
  ggsave(filename = filename, plot = plot, device = device, ...)
}
