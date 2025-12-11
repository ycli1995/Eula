
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Roxygen2 calls ###############################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.doc_links <- function(nm) {
  pkg <- .pkg_map[nm]
  paste0("\\code{\\link[", pkg, "]{", nm, "}}")
}

.dot_param <- "Arguments passed to other methods."
.val_param <- "An object of a class specified in the S4 method signature."
.vb_param <- "Logical, whether to print the progress information to console."

.pkg_map <- c(
  "validObject" = "methods:validObject",

  "SVT_SparseMatrix" = "SparseArray"
)
