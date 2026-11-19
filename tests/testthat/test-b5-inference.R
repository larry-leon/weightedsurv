# B5: inference details (I1, I2a, I2b, S1, S2, S4).

test_that("I1: the band critical value is the 1 - alpha quantile of the resampled suprema", {
  d <- ws_sim(120, 51)
  a <- KM_diff(d, "time", "status", "treat", draws.band = 200, alpha = 0.05, show_resamples = FALSE)
  b <- KM_diff(d, "time", "status", "treat", draws.band = 200, alpha = 0.20, show_resamples = FALSE)
  sups <- apply(abs(b$Zdhat_star), 2, max, na.rm = TRUE)
  expect_equal(unname(b$c_alpha_band), unname(stats::quantile(sups, 0.80)))
  expect_equal(unname(a$c_alpha_band), unname(stats::quantile(sups, 0.95)))
  expect_lt(b$c_alpha_band, a$c_alpha_band)
})

test_that("I2a: wlr_dhat_estimates() takes the pooled survival at the last grid point <= tzero", {
  dc <- df_counting(ws_sim(120, 52), "time", "status", "treat")
  t0 <- 12
  got <- wlr_dhat_estimates(dc, tzero = t0)
  dN <- diff(c(0, dc$nbar0 + dc$nbar1)); yb <- dc$ybar0 + dc$ybar1
  dJ <- ifelse(yb == 1, 0, (dN - 1) / (yb - 1)); dL <- ifelse(yb == 0, 0, dN / yb)
  h2 <- ifelse(dc$ybar0 * dc$ybar1 > 0, yb / (dc$ybar0 * dc$ybar1), 0) * (dc$at_points <= t0)
  i <- max(which(dc$at_points <= t0))
  expect_equal(got$sig2_dhat, dc$survP[i]^2 * sum(h2 * (1 - dJ) * dL))
  expect_equal(got$dhat, dc$surv1[i] - dc$surv0[i])
})

test_that("I2b: tzero before the first or after the last grid point never gives a length-0 result", {
  dc <- df_counting(ws_sim(120, 53), "time", "status", "treat")
  early <- wlr_dhat_estimates(dc, tzero = -1)
  expect_identical(early$dhat, 0)
  expect_identical(early$sig2_dhat, 0)
  late <- wlr_dhat_estimates(dc, tzero = 1e6)
  n <- length(dc$at_points)
  expect_length(late$dhat, 1)
  expect_equal(late$dhat, dc$surv1[n] - dc$surv0[n])
  expect_equal(late, wlr_dhat_estimates(dc, tzero = max(dc$at_points)))
})

test_that("S1: MB weights use S(t*) carried forward at/beyond the last time, and are always finite", {
  tp <- 1:10; S <- seq(0.95, 0.5, length.out = 10); S_left <- c(1, S[-10])
  w_eq <- wt.rg.S(S, "MB", mb_tstar = 10, tpoints = tp)
  expect_equal(w_eq, 1 / pmax(S_left, S[10]))
  expect_equal(wt.rg.S(S, "MB", mb_tstar = 50, tpoints = tp), w_eq)
  expect_equal(wt.rg.S(S, "MB", mb_tstar = 4.5, tpoints = tp), 1 / pmax(S_left, S[4]))
  w0 <- wt.rg.S(c(0.9, 0.5, 0, 0, 0), "MB", mb_tstar = 99, tpoints = 1:5)
  expect_true(all(is.finite(w0)))
})

test_that("S2: Schemper weights need survG (clear error); no partial matching of survG", {
  dc <- df_counting(ws_sim(120, 54), "time", "status", "treat")
  dc2 <- dc; dc2$survG <- NULL
  expect_error(wlr_dhat_estimates(dc2, scheme = "schemper", tzero = 10), "survG")
  expect_error(wlr_cumulative(dc2[c("at_points", "nbar0", "nbar1", "ybar0", "ybar1", "survP")], scheme = "schemper"), "survG")
  expect_equal(wlr_dhat_estimates(dc2, tzero = 10), wlr_dhat_estimates(dc, tzero = 10))
})

test_that("S4: cumulative_rmst_bands() asks for a fit made with draws.band > 0", {
  d <- ws_sim(80, 55)
  k <- KM_diff(d, "time", "status", "treat")
  expect_error(cumulative_rmst_bands(d, k, "time", "status", "treat", plot = FALSE), "draws.band > 0")
})
