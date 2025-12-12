
test_that("filterData character", {
  old <- c("AA", "BB", "CC", "aa", "bb", "cc")
  include <- c("AA", "BB", "bb", "DD")
  exclude <- c("CC", "aa", "dd")
  x0 <- sample(old, 20, replace = TRUE)

  x <- x0[filterData(x0, include = include)]
  included <- intersect(include, x)
  expect_contains(include, x)
  expect_true(all(x %in% include))

  x <- x0[filterData(x0, include = exclude, invert = TRUE)]
  filtered <- setdiff(old, x)
  expect_contains(exclude, filtered)
  expect_false(any(x %in% exclude))
})

test_that("filterData numeric", {
  x0 <- runif(100, 0, 100)
  include <- c(20, 50)
  exclude <- list(c(0, 20), c(70, 90))

  x <- x0[filterData(x0, include = include)]
  expect_true(all(x > 20 & x < 50))

  x <- x0[filterData(x0, include = exclude, invert = TRUE)]
  expect_true(all((x > 20 & x < 70) | x > 90))
})
