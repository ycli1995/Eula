
.warn <- function(..., call. = FALSE, immediate. = TRUE) {
  warning(..., call. = call., immediate. = immediate.)
}

load_pkg <- function(name, pkg_conf) {
  name <- name[1]
  if (!name %in% names(pkg_conf)) {
    stop("Package '", name, "' is not found.")
  }
  for (i in pkg_conf[[name]]$import) {
    if (paste0("package:", i) %in% search()) {
      pkgload::unload(i)
    }
    load_pkg(i, pkg_conf = pkg_conf)
  }
  pkgload::load_all(name)
  invisible(NULL)
}
