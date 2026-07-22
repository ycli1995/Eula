test_rowMeanPct <- function(
  object,
  cell.groups,
  mean.fxn = rowExpMean,
  ...
) {
  if (is.list(cell.groups)) {
    all.groups <- names(cell.groups)
  } else {
    cell.groups <- as.factor(cell.groups)
    all.groups <- levels(cell.groups)
  }
  for (min.exp in c(0, 0.05, 0.5)) {
    out <- rowMeanPct(
      object = object,
      cell.groups = cell.groups,
      mean.fxn = mean.fxn,
      min.exp = min.exp
    )
    for (i in seq_along(out)) {
      expect_equal(rownames(out[[i]]), rownames(object))
      expect_equal(colnames(out[[i]]), all.groups)
    }
  }
}

test_that("`rowMeanPct` for >1 levels.", {
  m <- random_matrix(20, 10)
  labels <- c("AA", "BB", "CC")
  g <- sample(labels, ncol(m), replace = TRUE)
  test_rowMeanPct(m, g)
})

test_that("`rowMeanPct` for 1 levels.", {
  m <- random_matrix(20, 10)
  g <- rep.int("AA", ncol(m))
  test_rowMeanPct(m, g)
})
