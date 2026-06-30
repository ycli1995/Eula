
#' @importFrom xgboost xgb.DMatrix xgb.load
#' @importFrom stats predict
#' @export
predictMalignantCells <- function(
  expr, 
  model.ref = NULL, 
  genes = NULL, 
  thres = 0.5
) {
  genes.path <- system.file(
    "model", "genes-scRNA-tcga.txt",
    package = packageName()
  )
  genes <- genes %||% readLines(genes.path)

  model.path <- system.file(
    "model", "sc_xgboost.model",
    package = packageName()
  )
  model.ref <- model.ref %||% xgb.load(model.path)

  comm.features <- intersect(genes, colnames(expr))
  expr <- as.matrix(expr[, comm.features, drop = FALSE])

  testdata <- xgb.DMatrix(expr)
  predict.label <- predict(model.ref, testdata)

  # store results
  cell.annotation <- data.frame(
    Malign.score = predict.label,
    Malign.type = "Non-malignant",
    row.names = rownames(expr)
  )
  cell.annotation$Malign.type[predict.label > thres] <- "Malignant"
  cell.annotation$Malign.type <- factor(cell.annotation$Malign.type)
  cell.annotation
}
