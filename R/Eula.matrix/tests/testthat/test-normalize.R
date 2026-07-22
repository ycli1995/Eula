
test_norm_func <- function(m, func, func.slow, ...) {
  ans <- func.slow(as.matrix(m), ...)

  m1 <- as.matrix(func(m, ...))
  m2 <- func(as.matrix(m), ...)
  expect_equal(m1, ans)
  expect_equal(m2, ans)
}

test_that("relativeCounts works", {
  rc_norm_slow <- function(object, scale.factor = 10000) {
    object <- as.matrix(object)
    Matrix::t(Matrix::t(object) / Matrix::colSums(object) * scale.factor)
  }

  m <- random_matrix(20, 10)
  test_norm_func(m, relativeCounts, rc_norm_slow)

  scale.factor <- log10(Matrix::colMeans(m))
  test_norm_func(m, relativeCounts, rc_norm_slow, scale.factor = scale.factor)
})

test_that("relativeCounts works with Seurat", {
  skip_if_not_installed("Seurat")

  m <- random_matrix(20, 10)

  seurat.func <- function(m, ...) {
    as.matrix(Seurat::RelativeCounts(m, verbose = FALSE, ...))
  }

  test_norm_func(m, relativeCounts, seurat.func, scale.factor = 10000)
})

test_that("logNorm works", {
  log_norm_slow <- function(object, scale.factor = 10000) {
    object <- as.matrix(object)
    object <- Matrix::t(
      Matrix::t(object) / Matrix::colSums(object) * scale.factor
    )
    object <- log1p(object)
    object
  }

  m <- random_matrix(20, 10)
  test_norm_func(m, logNorm, log_norm_slow)

  scale.factor <- log10(Matrix::colMeans(m))
  test_norm_func(m, logNorm, log_norm_slow, scale.factor = scale.factor)
})

test_that("logNorm works with Seurat", {
  skip_if_not_installed("Seurat")

  m <- random_matrix(20, 10)

  seurat.func <- function(m, ...) {
    as.matrix(Seurat::LogNormalize(m, verbose = FALSE, ...))
  }

  test_norm_func(m, logNorm, seurat.func, scale.factor = 10000)
})
