slow_rowExpMean <- function(object, log1p = FALSE, ...) {
  fun <- function(x) mean(exp(x) - 1)
  if (isTRUE(log1p)) {
    fun <- function(x) log(mean(exp(x) - 1) + 1)
  }
  if (inherits(object, c("matrix", "CsparseMatrix"))) {
    return(apply(X = object, MARGIN = 1, FUN = fun))
  }
  fun(object)
}

test_rowExpMean <- function(object, log1p = FALSE, ...) {
  expect_equal(
    rowExpMean(object, log1p = log1p),
    slow_rowExpMean(object, log1p = log1p),
    ignore_attr = TRUE
  )
}

test_that("`rowExpMean`", {
  m <- random_matrix(20, 10)

  for (log1p in c(TRUE, FALSE)) {
    test_rowExpMean(m, log1p = log1p)
    test_rowExpMean(as.matrix(m), log1p = log1p)

    for (n in seq_len(5)) {
      i <- sample(seq_len(nrow(m)), n)
      test_rowExpMean(m[i, , drop = FALSE], log1p = log1p)
      test_rowExpMean(m[i, , drop = TRUE], log1p = log1p)
    }
  }
})
