
#' @importFrom rlang %||%
#' @importFrom vctrs %0%
#' @keywords internal
"_PACKAGE"

#' @importFrom ggplot2 el_def register_theme_elements
.onLoad <- function(libname, pkgname) {
  ggplot2::register_theme_elements(
    Eula.corner.axis.title = .corner.axis.title.element,
    Eula.corner.axis = .corner.axis.element,

    element_tree = list(
      Eula.corner.axis.title = ggplot2::el_def("element_text", "text"),
      Eula.corner.axis = ggplot2::el_def("element_line", "line")
    )
  )
}
