
#' @importFrom ggplot2 element_blank element_rect element_text theme theme_bw
#' @export
dot_theme_default <- function(family = "", ...) {
  theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "#000000", size = 0.8, fill = NA),

      axis.text = element_text(color = "#000000", size = 11),
      axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
      axis.text.y = element_text(hjust = 0.5, vjust = 0.5),
      axis.title = element_text(color = "#000000", size = 14, face = "plain"),
      axis.title.x = element_text(margin = margin(2.5,0,2.5,0, "mm")),
      axis.title.y = element_text(margin = margin(0,2.5,0,2.5, "mm")),
      axis.ticks = element_line(color = "#000000", size = 0.5),
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

