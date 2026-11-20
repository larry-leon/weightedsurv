# Rc: plotKM.band_subgroups() default draws.band = 1000 (it was 20); nothing else changes.

test_that("Rc: the default draws.band is 1000 and the default call equals the explicit one", {
  expect_identical(formals(plotKM.band_subgroups)$draws.band, 1000)
  d <- ws_sim(150, 21)
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  a <- plotKM.band_subgroups(d, "time", "status", "treat")
  b <- plotKM.band_subgroups(d, "time", "status", "treat", draws.band = 1000)
  expect_identical(a, b)
  expect_identical(ncol(a$fit_itt$dhat_star), 1000L)
  k <- KM_diff(d, "time", "status", "treat", at_points = a$fit_itt$at_points, draws.band = 1000)
  expect_identical(unname(a$fit_itt$c_alpha_band), unname(k$c_alpha_band))
})

test_that("Rc: a call that supplies draws.band is as before", {
  d <- ws_sim(150, 21)
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  a <- plotKM.band_subgroups(d, "time", "status", "treat", draws.band = 20)
  k <- KM_diff(d, "time", "status", "treat", at_points = a$fit_itt$at_points, draws.band = 20)
  expect_identical(a$fit_itt$dhat_star, k$dhat_star)
  expect_identical(ncol(a$fit_itt$dhat_star), 20L)
})
