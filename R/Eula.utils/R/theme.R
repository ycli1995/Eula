
#' @importFrom ggplot2 element_blank element_line element_rect element_text
#' margin theme theme_bw
#' @importFrom grid unit
#' @export
dot_theme_default <- function(family = "Arial", ...) {
  theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(
        color = "#000000",
        linewidth = 0.8,
        fill = NA
      ),

      axis.text = element_text(color = "#000000", size = 11),
      axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
      axis.text.y = element_text(hjust = 0.5, vjust = 0.5),
      axis.title = element_text(color = "#000000", size = 14, face = "plain"),
      axis.title.x = element_text(margin = margin(2.5,0,2.5,0, "mm")),
      axis.title.y = element_text(margin = margin(0,2.5,0,2.5, "mm")),
      axis.ticks = element_line(color = "#000000", linewidth = 0.5),
      axis.ticks.length = unit(0.1, 'cm'),

      legend.title = element_blank(),
      legend.text = element_text(size = 11),
      plot.title = element_text(
        color = "#000000",
        size = 16,
        face = "plain",
        hjust = 0.5
      ),
      plot.margin = unit(c(10, 10, 10, 10), "mm")
    )  +
    theme(text = element_text(family = family))
}

#' @export
bar_theme_default <- function(family = "Arial", ...) {
  options(scipen = -1)
  theme_bw() +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(
        color = "#000000",
        linewidth = 0.8,
        fill = NA
      ),

      axis.text  = element_text(color = "#000000", size = 11),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.text.y = element_text(hjust = 1, vjust = 0.5),
      axis.title = element_text(color = "#000000", size = 14, face = "plain"),
      axis.title.x = element_text(margin = margin(2.5,0,2.5,0, "mm")),
      axis.title.y = element_text(margin = margin(0,2.5,0,2.5, "mm")),
      axis.ticks = element_line(color = "#000000", linewidth = 0.5),
      axis.ticks.length = unit(0.1, 'cm'),

      plot.title = element_text(size = 16, face = "plain", hjust = 0.5),
      plot.margin = unit(c(10,10,10,10), "mm")
    ) +
    theme(text = element_text(family = family))
}

#' @importFrom ggplot2 get_last_plot ggsave
#' @importFrom grDevices cairo_pdf
#' @export
ggsave2 <- function(filename, plot = get_last_plot(), device = NULL, ...) {
  ext <- tools::file_ext(filename)
  if (tolower(ext) == "pdf") {
    device <- device %||% cairo_pdf
  }
  ggsave(filename = filename, plot = plot, device = device, ...)
}
