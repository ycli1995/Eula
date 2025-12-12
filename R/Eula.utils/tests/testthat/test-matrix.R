
test_that("asMatrix", {
  mat <- Matrix::rsparsematrix(20, 10, density = 0.5)
  rownames(mat) <- letters[1:20]
  colnames(mat) <- LETTERS[1:10]

  expect_equal(asMatrix(mat, chunk.size = 10, verbose = FALSE), as.matrix(mat))
})
