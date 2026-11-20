# B12: plotKM.band_subgroups() default draws = 0 (the pointwise variance is analytic unless draws is supplied).

pkb <- function(d, ...) { grDevices::pdf(NULL); on.exit(grDevices::dev.off()); plotKM.band_subgroups(d, "time", "status", "treat", ...) }

test_that("B12: the default of draws is 0 and nothing else in the signature changed", {
  expect_identical(formals(plotKM.band_subgroups)$draws, 0)
  expect_identical(formals(plotKM.band_subgroups)$draws.band, 20)
})

test_that("B12: with weights and draws not supplied, the SE and the band standardisation are analytic", {
  d <- ws_sim_ipw(200, 81)
  r <- pkb(d, weight.name = "w", draws.band = 100)$fit_itt
  k <- KM_diff(d, "time", "status", "treat", weight.name = "w", at_points = r$at_points, draws = 0, draws.band = 100)
  expect_identical(r$sig2_dhat, k$sig2_dhat)
  expect_identical(r$lower, k$lower); expect_identical(r$upper, k$upper)
  expect_identical(unname(r$c_alpha_band), unname(k$c_alpha_band))
  expect_identical(r$sb_lower, k$sb_lower); expect_identical(r$sb_upper, k$sb_upper)
  # draws supplied: the resampling variance, as before
  r20 <- pkb(d, weight.name = "w", draws = 20, draws.band = 100)$fit_itt
  k20 <- KM_diff(d, "time", "status", "treat", weight.name = "w", at_points = r$at_points, draws = 20, draws.band = 100)
  expect_identical(r20$sig2_dhat, k20$sig2_dhat)
  expect_false(isTRUE(all.equal(r20$sig2_dhat, r$sig2_dhat)))
})

test_that("B12: without weights (or with equal weights) draws makes no difference", {
  d <- ws_sim_ipw(200, 82, ties = TRUE)
  a <- pkb(d, draws.band = 100); b <- pkb(d, draws = 20, draws.band = 100)
  expect_identical(a, b)
  expect_identical(pkb(d, weight.name = "one", draws.band = 100), pkb(d, weight.name = "one", draws = 20, draws.band = 100))
})
