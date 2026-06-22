test_that("asMatrix", {
  mat <- random_matrix(20, 10)

  expect_equal(asMatrix(mat, chunk.size = 10, verbose = FALSE), as.matrix(mat))
})
