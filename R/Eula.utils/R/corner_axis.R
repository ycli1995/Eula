
#' @importFrom ggplot2 annotate facet_null geom_text get_labs get_panel_scales
#' lims rel theme_bw
#' @export
add_corner_axis <- function(
    plot,
    length = 0.2,
    lab.size = 3.5,
    fontface = 1,
    family = "Arial",
    clip = "off",
    ...
) {
  tmp.plot <- plot + ggplot2::facet_null()
  xy <- ggplot2::get_panel_scales(plot = tmp.plot)
  x.range <- xy$x$get_limits()
  y.range <- xy$y$get_limits()
  x.diff <- diff(x.range)
  y.diff <- diff(y.range)

  x.lim <- c(x.range[1] - 0.05 * x.diff, x.range[2])
  y.lim <- c(y.range[1] - 0.05 * y.diff, y.range[2])

  plot.labs <- ggplot2::get_labs(plot)

  plot +
    ggplot2::lims(x = x.lim, y = y.lim) +
    Eula_corner_axis(
      length = length,
      x.lab = plot.labs$x,
      y.lab = plot.labs$y,
      x.lim = x.lim,
      y.lim = y.lim,
      clip = clip,
    ) +
    theme_no_axis(plot.margin = unit(c(10, 10, 10, 10), "mm"))
}

#' @importFrom ggplot2 CoordCartesian ggproto element_render
#' @export
Eula_corner_axis <- function(
    length = 0.2,
    x.lab = "x",
    y.lab = "y",
    x.lim = NULL,
    y.lim = NULL,
    clip = "off",
    ...
) {
  ggproto(
    NULL,
    CoordCartesian,
    limits = list(x = x.lim, y = y.lim),
    expand = TRUE,
    default = FALSE,
    clip = clip,

    render_fg = function(panel_params, theme) {
      element_render(
        theme = theme,
        element = "Eula.corner.axis.title",
        label = c(x.lab, y.lab),
        x = unit(c(0.05, 0.025), "npc"),
        y = unit(c(0.01, 0.06), "npc"),
        angle = c(0, 90),
        hjust = c(0, 0),
        vjust = c(0, 0.5)
      )
    },
    render_bg = function(self, panel_params, theme) {
      element_render(
        theme = theme,
        element = "Eula.corner.axis",
        x = unit(c(0.05 + length, 0.05, 0.05), "npc"),
        y = unit(c(0.05, 0.05, 0.05 + length), "npc")
      )
    }
  )
}

#' @importFrom ggplot2 element_line
#' @importFrom grid arrow
.corner.axis.element <- element_line(
  arrow = grid::arrow(
    length = unit(0.15, "inches"),
    type = "closed",
    ends = "both",
    angle = 15
  )
)

#' @importFrom ggplot2 element_text
.corner.axis.title.element <- element_text(
  color = "black"
)



