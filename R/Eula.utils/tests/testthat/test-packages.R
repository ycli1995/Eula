
test_that("loadPackages", {
  pkgs <- c("Matrix", "dplyr")
  loadPackages(pkgs, verbose = FALSE)
  expect_in(paste0("package:", pkgs), search())

  miss.pkgs <- c("pkg1", "pkg2")
  expect_error(loadPackages(miss.pkgs, verbose = FALSE))
})
