
test_that("checkColumns", {
  data("cars")
  expect_invisible(checkColumns(cars, c("dist")))
  expect_invisible(checkColumns(cars, c("speed", "dist")))
  expect_error(checkColumns(cars, cols = NULL))
  expect_error(checkColumns(cars, c("speed", "dist", "col1")))
  expect_error(checkColumns(cars, c("dist", "col1")))
  expect_error(checkColumns(cars, c("col1", "col2")))
})

test_that("checkKeys", {
  keys <- c("A", "BB", "c")

  expect_error(checkKeys(character(), keys))
  expect_error(checkKeys(1:10, keys))
  expect_error(checkKeys(NULL, keys))

  ll <- list(A = NULL, BB = 1:2, c = character())
  expect_error(checkKeys(ll, keys))

  ll$A <- c("A", "B")
  expect_error(checkKeys(ll, keys))

  ll$c <- NA
  expect_invisible(checkKeys(ll, keys))
})

test_that("validCharacters", {
  expect_false(validCharacters(NA))
  expect_false(validCharacters(NA_character_))

  x <- c('', NA_character_, "A")
  y <- c(FALSE, FALSE, TRUE)
  expect_length(validCharacters(x), length(x))
  expect_equal(validCharacters(x), y)
})

test_that("checkPackages", {
  ok.pkgs <- c("base", "utils")
  miss.pkgs <- c("not_found1", "not_found2")

  expect_invisible(checkPackages(ok.pkgs))
  expect_null(checkPackages(ok.pkgs))

  expect_error(checkPackages(miss.pkgs))
  expect_warning(checkPackages(miss.pkgs, required = FALSE))
})
