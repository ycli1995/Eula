
test_that("recodeFactor", {
  x <- c(rep.int("cluster1", 3), rep.int("cluster2", 2), rep.int("C3", 4))
  x <- as.factor(x)

  # one-to-one
  map1 <- list(c2 = "cluster2", c1 = "cluster1", c3 = "cluster3")
  x1 <- recodeFactor(x, map1)
  expect_equal(levels(x1), c("c2", "c1", "C3"))
  x1 <- recodeFactor(x, map1, keep.order = TRUE)
  expect_equal(levels(x1), c("C3", "c1", "c2"))

  # one-to-many
  map1 <- list(c2 = c("cluster2", "cluster1"), c3 = "C3")
  x1 <- recodeFactor(x, map1)
  expect_equal(levels(x1), c("c2", "c3"))
  x1 <- recodeFactor(x, map1, keep.order = TRUE)
  expect_equal(levels(x1), c("c3", "c2"))

  # Duplicated old levels
  map1 <- list(c2 = c("cluster2", "C3"), c3 = "C3")
  x1 <- recodeFactor(x, map1)
  expect_equal(levels(x1), c("c2", "cluster1"))
  x1 <- recodeFactor(x, map1, keep.order = TRUE)
  expect_equal(levels(x1), c("c2", "cluster1"))
})

test_that("pasteFactors", {
  x <- sample(c("A", "B"), 10, replace = TRUE)
  y <- sample(c("a", "b"), 10, replace = TRUE)

  xy <- pasteFactors(x, y)
  expect_s3_class(xy, "factor")
  expect_equal(levels(xy), c("A_a", "A_b", "B_a", "B_b"))

  xy <- pasteFactors(x, y, rev = TRUE)
  expect_equal(levels(xy), c("A_a", "B_a", "A_b", "B_b"))
})
