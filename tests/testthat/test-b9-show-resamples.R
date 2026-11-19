# B9: KM_diff() no longer plots the resampled curves unless asked to.

test_that("show_resamples defaults to FALSE", {
  expect_false(formals(KM_diff)$show_resamples)
})

test_that("draws.band > 0 draws nothing by default and the resamples when asked; results do not depend on it", {
  d <- ws_sim(120, 91)
  n_mat <- 0
  local_mocked_bindings(matplot = function(...) n_mat <<- n_mat + 1)
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  a <- KM_diff(d, "time", "status", "treat", draws.band = 50)
  expect_equal(n_mat, 0)
  b <- KM_diff(d, "time", "status", "treat", draws.band = 50, show_resamples = TRUE)
  expect_equal(n_mat, 1)
  expect_identical(a, b)
})
