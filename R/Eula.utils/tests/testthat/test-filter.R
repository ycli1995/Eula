
test_that("filterData character", {
  old <- c("AA", "BB", "CC", "aa", "bb", "cc")
  include <- c("AA", "BB", "bb", "DD")
  exclude <- c("CC", "aa", "dd")
  x0 <- rep(old, 3)

  # test `include`
  expect_true(all(filterData(x0, include = NULL)))
  expect_false(any(filterData(x0, include = c("XX", "YY"))))

  x <- x0[filterData(x0, include = include)]
  included <- intersect(include, x)
  expect_contains(include, x)
  expect_true(all(x %in% include))
  expect_true(all(filterData(x, include = include)))

  # test `invert`
  expect_false(any(filterData(as.factor(x0), include = NULL, invert = TRUE)))
  expect_true(all(filterData(x0, include = c("XX", "YY"), invert = TRUE)))

  x <- x0[filterData(as.factor(x0), include = exclude, invert = TRUE)]
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

  # Empty `include`, select all
  expect_true(all(filterData(x0, include = NULL)))

  # Invalid `include`
  wrong.include <- list(2, c(0, 50, 70), numeric(), NULL)
  expect_error(filterData(x0, include = wrong.include))
})

test_that("filterData data.frame", {
  df <- data.frame(
    num1 = runif(100, 0, 100),
    num2 = runif(100, 0, 50),
    logi = sample(c(TRUE, FALSE), 100, replace = TRUE),
    char = sample(LETTERS[1:5], 100, replace = TRUE)
  )

  expect_equal(filterData(df, "logi", verbose = FALSE), df[["logi"]])
  expect_warning(
    expect_true(all(filterData(df, "miss column", verbose = FALSE)))
  )

  include <- c("A", "B")
  df2 <- df[filterData(df, "char", include = include, verbose = FALSE), ]
  expect_contains(df2[["char"]], include)
})
