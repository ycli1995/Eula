test_that("relativeCounts works", {
  rc_norm_slow <- function(object, scale.factor = 10000) {
    object <- as.matrix(object)
    Matrix::t(Matrix::t(object) / Matrix::colSums(object) * scale.factor)
  }

  m <- random_matrix(20, 10)

  ans <- rc_norm_slow(as.matrix(m))
  m1 <- as.matrix(relativeCounts(m))
  m2 <- relativeCounts(as.matrix(m))
  expect_equal(m1, ans)
  expect_equal(m2, ans)

  scale.factor <- log10(Matrix::colMeans(m))
  ans <- rc_norm_slow(as.matrix(m), scale.factor = scale.factor)
  m1 <- as.matrix(relativeCounts(m, scale.factor = scale.factor))
  m2 <- relativeCounts(as.matrix(m), scale.factor = scale.factor)
  expect_equal(m1, ans)
  expect_equal(m2, ans)
})

test_that("relativeCounts works with Seurat", {
  skip_if_not_installed("Seurat")

  m <- random_matrix(20, 10)

  ans <- Seurat::RelativeCounts(m, scale.factor = 10000, verbose = FALSE)
  m1 <- as.matrix(relativeCounts(m))
  m2 <- relativeCounts(as.matrix(m))
  expect_equal(m1, as.matrix(ans))
  expect_equal(m2, as.matrix(ans))
})
