# # This env must be attached first to be base of S4 Methods
# env.name <- "Eula.S4Generics"
# Eula.S4.env <- new.env()
#
# setGeneric(
#   "recodeFactor",
#   function(x, map, ...) standardGeneric("recodeFactor"),
#   where = Eula.S4.env
# )
#
# setGeneric(
#   "checkColorMap",
#   function(x, ...) standardGeneric("checkColorMap"),
#   where = Eula.S4.env
# )
#
# setGeneric(
#   "statByGroups",
#   function(x, group.by, ...) standardGeneric("statByGroups"),
#   where = Eula.S4.env
# )
#
# setGeneric(
#   "statCrossTab",
#   function(x, group.by, stat.what, ...) standardGeneric("statCrossTab"),
#   where = Eula.S4.env
# )
#
# setGeneric(
#   "plotStatCrossTab",
#   function(x, group.by, stat.what, ...) standardGeneric("plotStatCrossTab"),
#   where = Eula.S4.env
# )
#
# setGeneric(
#   "setFData",
#   function(object, fdata, ...) standardGeneric("setFData"),
#   where = Eula.S4.env
# )
#
# setGeneric(
#   "getFeaturesID",
#   function(object, features, ...) standardGeneric("getFeaturesID"),
#   where = Eula.S4.env
# )
#
# setGeneric(
#   "getFeaturesName",
#   function(object, features, ...) standardGeneric("getFeaturesName"),
#   where = Eula.S4.env
# )
#
# setGeneric(
#   "DoAvgExp",
#   function(object, group.by, ...) standardGeneric("DoAvgExp"),
#   where = Eula.S4.env
# )
#
# if (env.name %in% search()) {
#   detach(env.name, character.only = TRUE)
# }
# attach(Eula.S4.env, name = env.name)

recodeFactor <- function(x, map, ...) {
  UseMethod("recodeFactor", x)
}


