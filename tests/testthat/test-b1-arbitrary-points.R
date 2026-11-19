# B1: KM_diff() returns the Kaplan-Meier value at any valid at_points.

sf_at <- function(d, g, times, w = FALSE) {
  dd <- d[d$treat == g, ]
  fit <- if (w) survival::survfit(survival::Surv(time, status) ~ 1, data = dd, weights = w)
         else survival::survfit(survival::Surv(time, status) ~ 1, data = dd)
  summary(fit, times = times, extend = TRUE)
}

test_that("surv0/surv1 equal survfit and unweighted variances equal Greenwood at arbitrary points", {
  d <- ws_sim(150, 31)
  grids <- list(seq(0, 30, 5), c(2.5, 2.5, 7, 7, 7, 40), seq(3.3, 12, by = 0.7), 12)
  for (ap in grids) {
    k <- KM_diff(d, "time", "status", "treat", at_points = ap)
    expect_identical(k$at_points, ap)
    for (g in 0:1) {
      ref <- sf_at(d, g, ap)
      expect_equal(if (g == 0) k$surv0 else k$surv1, ref$surv, tolerance = 1e-10)
      expect_equal(sqrt(if (g == 0) k$sig2_surv0 else k$sig2_surv1), ref$std.err, tolerance = 1e-10)
    }
  }
})

test_that("weighted: surv equals survfit(weights =) at arbitrary points", {
  d <- ws_sim(150, 32, ties = TRUE, weighted = TRUE)
  ap <- c(1, 4.2, 9.9, 15)
  k <- KM_diff(d, "time", "status", "treat", weight.name = "w", at_points = ap)
  for (g in 0:1) expect_equal(if (g == 0) k$surv0 else k$surv1, sf_at(d, g, ap, w = TRUE)$surv, tolerance = 1e-10)
})

test_that("modify_tau = TRUE: curves equal survfit at the grid", {
  d <- ws_sim(120, 33)
  k <- KM_diff(d, "time", "status", "treat", draws.band = 50, show_resamples = FALSE, modify_tau = TRUE)
  for (g in 0:1) expect_equal(if (g == 0) k$surv0 else k$surv1, sf_at(d, g, k$at_points)$surv, tolerance = 1e-10)
})

test_that("unsorted, NA or non-finite at_points are errors", {
  d <- ws_sim(60, 34)
  expect_error(KM_diff(d, "time", "status", "treat", at_points = c(5, 2, 9)), "non-decreasing")
  expect_error(KM_diff(d, "time", "status", "treat", at_points = c(1, NA, 9)), "finite")
  expect_error(KM_diff(d, "time", "status", "treat", at_points = c(1, Inf)), "finite")
})

test_that("default grid: identical to the direct calculation on at_points", {
  d <- ws_sim(100, 35, weighted = TRUE)
  k <- KM_diff(d, "time", "status", "treat", weight.name = "w")
  ds <- d[order(d$time), ]; g <- ds[ds$treat == 1, ]   # KM_diff works on time-sorted data
  re <- calculate_risk_event_counts(g$time, g$status, g$w, sort(d$time))
  km <- KM_estimates(re$ybar, re$nbar, re$sig2w_multiplier)
  expect_identical(k$surv1, km$S_KM)
  expect_identical(k$sig2_surv1, km$sig2_KM)
})

test_that("resampled processes on a sub-grid are the rows of those on the full grid", {
  d <- ws_sim(100, 36)
  full <- KM_diff(d, "time", "status", "treat", draws.band = 40, show_resamples = FALSE)
  sub_pts <- full$at_points[seq(1, length(full$at_points), by = 7)]
  sub <- KM_diff(d, "time", "status", "treat", at_points = sub_pts, draws.band = 40, show_resamples = FALSE)
  rows <- match(sub_pts, full$at_points)
  expect_equal(sub$surv0_star, full$surv0_star[rows, , drop = FALSE], tolerance = 1e-12)
  expect_equal(sub$dhat_star, full$dhat_star[rows, , drop = FALSE], tolerance = 1e-12)
})
