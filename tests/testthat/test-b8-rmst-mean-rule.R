# B8: cumulative_rmst_bands(mean_rule = c("trapezoid", "rectangle")).

b8_data <- function() {
  g <- survival::gbsg
  list(veteran = ws_veteran(),
       gbsg = data.frame(time = g$rfstime / 30.4375, status = g$status, treat = g$hormon),
       sim = ws_sim(200, 81, ties = TRUE))
}
sf_rmean_diff <- function(d, tau) {
  fit <- survival::survfit(survival::Surv(time, status) ~ treat, data = d)
  rm <- summary(fit, rmean = tau)$table[, "rmean"]
  unname(rm[2] - rm[1])
}
check_rect_vs_survfit <- function(d, fit, label) {
  r <- cumulative_rmst_bands(d, fit, "time", "status", "treat", draws_sb = 40, plot = FALSE, mean_rule = "rectangle")
  ap <- r$at_points
  ev <- sort(unique(d$time[d$status == 1 & d$time <= max(ap)]))
  expect_true(all(ev %in% ap), info = label)  # the grid holds every event time
  k <- length(ap)
  j <- unique(c(round(k * c(0.25, 0.5, 0.75)), k))
  for (jj in j) expect_lte(abs(r$rmst_time[jj] - sf_rmean_diff(d, ap[jj])), 1e-10)
  r
}

test_that("mean_rule is the last argument, defaults to the trapezoid, and is returned", {
  f <- formals(cumulative_rmst_bands)
  expect_identical(tail(names(f), 1), "mean_rule")
  expect_identical(eval(f$mean_rule), c("trapezoid", "rectangle"))
  d <- ws_sim(120, 82)
  k <- KM_diff(d, "time", "status", "treat", draws.band = 30)
  a <- cumulative_rmst_bands(d, k, "time", "status", "treat", draws_sb = 30, plot = FALSE)
  b <- cumulative_rmst_bands(d, k, "time", "status", "treat", draws_sb = 30, plot = FALSE, mean_rule = "trapezoid")
  expect_identical(a, b)
  expect_identical(a$mean_rule, "trapezoid")
  # the trapezoid estimate is the average of consecutive values times the step, from 0
  dt <- diff(k$at_points)
  expect_equal(a$rmst_time, c(0, cumsum((head(k$dhat, -1) + tail(k$dhat, -1)) / 2 * dt)), tolerance = 0)
  expect_error(cumulative_rmst_bands(d, k, "time", "status", "treat", plot = FALSE, mean_rule = "simpson"), "should be one of")
})

test_that("rectangle equals survfit's restricted-mean difference on the default KM_diff grid", {
  D <- b8_data()
  for (nm in names(D)) {
    d <- D[[nm]]
    k <- KM_diff(d, "time", "status", "treat", draws.band = 40)
    check_rect_vs_survfit(d, k, nm)
  }
})

test_that("rectangle equals survfit's restricted-mean difference on fit_itt from plotKM.band_subgroups()", {
  D <- b8_data()
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  for (nm in names(D)) {
    d <- D[[nm]]
    p <- plotKM.band_subgroups(d, "time", "status", "treat", draws.band = 40)
    check_rect_vs_survfit(d, p$fit_itt, paste(nm, "fit_itt"))
  }
})

test_that("the rule applies to the estimate and to both resample sets; draws and band machinery unchanged", {
  d <- ws_sim(150, 83)
  k <- KM_diff(d, "time", "status", "treat", draws.band = 50)
  rt <- cumulative_rmst_bands(d, k, "time", "status", "treat", draws_sb = 60, plot = FALSE)
  rr <- cumulative_rmst_bands(d, k, "time", "status", "treat", draws_sb = 60, plot = FALSE, mean_rule = "rectangle")
  dt <- diff(k$at_points)
  rect <- function(x) c(0, cumsum(head(x, -1) * dt))
  expect_equal(rr$rmst_time, rect(k$dhat), tolerance = 0)
  # variance from the first draw set (the fit's resampled curves), integrated with the same rule
  s1 <- apply(apply(k$dhat_star, 2, rect), 1, var, na.rm = TRUE)
  expect_equal(rr$sig2_rmst_time, s1, tolerance = 1e-12)
  expect_equal(rr$rmst_time_lower, rr$rmst_time - 1.96 * sqrt(s1), tolerance = 1e-12)
  # second draw set: same KM_diff call as under the trapezoid (same seed and draws)
  k2 <- KM_diff(d, "time", "status", "treat", at_points = k$at_points, alpha = 0.05, risk.points = k$risk.points,
                modify_tau = FALSE, draws.band = 60, seedstart = 99999, show_resamples = FALSE)
  sups <- apply(abs(apply(k2$dhat_star, 2, rect) / sqrt(s1)), 2, max, na.rm = TRUE)
  expect_equal(rr$c_alpha_band, quantile(sups, 0.95), tolerance = 1e-12)
  expect_equal(rr$rmst_time_sb_upper, rr$rmst_time + rr$c_alpha_band * sqrt(s1), tolerance = 1e-12)
  expect_identical(rr$at_points, rt$at_points)
  expect_identical(rr$mean_rule, "rectangle")
})
