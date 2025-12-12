
test_that("asMatrix", {
  mat <- Matrix::rsparsematrix(20, 10, density = 0.5)
  rownames(mat) <- letters[1:20]
  colnames(mat) <- LETTERS[1:10]

  expect_equal(asMatrix(mat, chunk.size = 10, verbose = FALSE), as.matrix(mat))
})

test_that("rowExpMean", {
  mat <- Matrix::rsparsematrix(20, 10, density = 0.5)

  x1 <- rowExpMean(mat)
  x2 <- rowExpMean(as.matrix(mat))
  expect_length(x1, nrow(mat))
  expect_equal(x1, x2)

  x1 <- rowExpMean(mat, log = TRUE)
  expect_equal(expm1(x1), x2)
})

test_that("rowMeanPct", {
  mat <- Matrix::rsparsematrix(10, 20, density = 0.5)
  cell.groups <- sample(LETTERS[1:3], ncol(mat), replace = TRUE)
  out <- rowMeanPct(mat, cell.groups)

  out.names <- c("n.cells", "avg.exp", "avg.pct")
  expect_true(is.list(out))
  expect_length(out, 3)
  expect_named(out, out.names)
  for (i in out.names) {
    expect_true(is.matrix(out[[i]]))
    expect_equal(dim(out[[i]]), c(nrow(mat), length(unique(cell.groups))))
  }
})
