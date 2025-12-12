
test_that("readRDX", {
  rds.file <- withr::local_tempfile(fileext = ".Rds")
  rda.file <- withr::local_tempfile(fileext = ".Rda")

  aa <- list(A = 1, B = 2, C = letters[1:3], D = logical())
  bb <- list(a = 3, b = 4, c = LETTERS[1:5], d = NULL)

  saveRDS(aa, rds.file)
  aa2 <- readRDX(rds.file, verbose = FALSE)
  expect_equal(aa2, aa)

  save(aa, file = rda.file)
  aa3 <- readRDX(rda.file, verbose = FALSE)
  expect_equal(aa3, aa)

  # Only fetch the first object in an .Rda file
  save(aa, bb, file = rda.file)
  aa3 <- readRDX(rda.file, verbose = FALSE)
  expect_equal(aa3, aa)
})

test_that("readTable", {
  df.file <- withr::local_tempfile(fileext = ".Rds")

  df <- data.frame(
    Cells = letters[1:10],
    n1 = runif(10) * 100,
    n2 = 1:10,
    "logi var" = sample(c(TRUE, FALSE), 10, replace = TRUE),
    char = sample(letters, 10)
  )
  writeTable(df, df.file)

  df2 <- readTable(df.file)
  expect_equal(df2, df)

  colnames(df)[4] <- "logi var"
  writeTable(df, df.file)

  df2 <- readTable(df.file)
  expect_equal(df2, df)

  rownames(df) <- df[['Cells']]
  df2 <- readTable(df.file, row.names = "Cells")
  expect_equal(df2, df)

  df2 <- readTable(df.file, row.names = "Cells", keep.row.names = FALSE)
  expect_equal(df2, df[, -1])
})


