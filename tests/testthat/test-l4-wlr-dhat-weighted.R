# L4: wlr_dhat_estimates() with subject weights -- log-rank variance, variance of the survival
# difference and their covariance from the same per-subject contributions.

dcw <- function(d, w = "w", ...) ws_quiet(df_counting(d, "time", "status", "treat", weight.name = w, ...))
rel <- function(a, b) abs(a / b - 1)

test_that("L4: sig2_dhat equals KM_diff(draws = 0) at the same grid point; z as df_counting(); |cor| <= 1", {
  for (ties in c(FALSE, TRUE)) {
    d <- ws_sim_ipw(300, 71, ties)
    for (sc in ws_schemes(d)[c("fh00", "fh01")]) {
      dc <- dcw(d, scheme = sc[[1]], scheme_params = sc[[2]])
      km <- KM_diff(d, "time", "status", "treat", weight.name = "w", at_points = dc$at_points, draws = 0)
      for (t0 in c(5, 18, stats::median(dc$at_points), 40)) {
        g <- wlr_dhat_estimates(dc, scheme = sc[[1]], scheme_params = sc[[2]], tzero = t0)
        j <- max(which(dc$at_points <= t0))
        expect_lt(rel(g$sig2_dhat, km$sig2_dhat[j]), 1e-8)
        expect_equal(g$dhat, km$dhat[j], tolerance = 1e-12)
        expect_identical(g$score, dc$lr)
        expect_lt(rel(g$sig2.score, dc$sig2_lr), 1e-12)
        expect_lt(rel(g$z.score, dc$z.score), 1e-12)
        expect_lte(abs(g$cor_wlr_dhat), 1)
        expect_equal(g$cor_wlr_dhat, g$cov_wlr_dhat / sqrt(g$sig2.score * g$sig2_dhat))
      }
      # without tzero: the log-rank part only
      g0 <- wlr_dhat_estimates(dc, scheme = sc[[1]], scheme_params = sc[[2]])
      expect_named(g0, c("score", "sig2.score", "z.score"))
      expect_lt(rel(g0$sig2.score, dc$sig2_lr), 1e-12)
    }
  }
})

test_that("L4: the covariance is that of the per-subject contributions (checked by weight derivatives)", {
  # FH(0,0): a_i * d(score)/d(a_i) = -psi_i(0) and a_i * d(dhat)/d(a_i) = d_i, so the covariance is
  # sum_i a_i^2 * d(score)/d(a_i) * d(dhat)/d(a_i); the derivatives are taken numerically.
  d <- ws_sim_ipw(60, 72); t0 <- 15
  stat <- function(w) {
    dd <- d; dd$w <- w
    dc <- ws_quiet(df_counting(dd, "time", "status", "treat", weight.name = "w", get.cox = FALSE, check.KM = FALSE))
    j <- max(which(dc$at_points <= t0))
    c(dc$lr, dc$surv1[j] - dc$surv0[j])
  }
  J <- t(vapply(seq_len(nrow(d)), function(i) {
    h <- 1e-5 * d$w[i]; up <- d$w; up[i] <- up[i] + h; dn <- d$w; dn[i] <- dn[i] - h
    (stat(up) - stat(dn)) / (2 * h) * d$w[i]
  }, numeric(2)))
  g <- wlr_dhat_estimates(dcw(d), tzero = t0)
  expect_equal(g$sig2.score, sum(J[, 1]^2), tolerance = 1e-6)
  expect_equal(g$sig2_dhat, sum(J[, 2]^2), tolerance = 1e-6)
  expect_equal(g$cov_wlr_dhat, sum(J[, 1] * J[, 2]), tolerance = 1e-6)
})

test_that("L4: tzero before the first grid point gives 0 variance and covariance; scale-free correlation", {
  d <- ws_sim_ipw(200, 73); d$w10 <- 10 * d$w
  dc <- dcw(d)
  early <- wlr_dhat_estimates(dc, tzero = -1)
  expect_identical(early$dhat, 0); expect_identical(early$sig2_dhat, 0); expect_identical(early$cov_wlr_dhat, 0)
  expect_true(is.na(early$cor_wlr_dhat))
  a <- wlr_dhat_estimates(dc, tzero = 18); b <- wlr_dhat_estimates(dcw(d, "w10"), tzero = 18)
  expect_lt(rel(b$cor_wlr_dhat, a$cor_wlr_dhat), 1e-9)
  expect_lt(rel(b$sig2_dhat, a$sig2_dhat), 1e-9)
  expect_lt(rel(b$z.score, a$z.score), 1e-9)
})

test_that("L4: no weights and equal weights keep the earlier formulas", {
  d <- ws_sim_ipw(200, 74, ties = TRUE)
  u <- ws_quiet(df_counting(d, "time", "status", "treat")); c1 <- dcw(d, "one")
  expect_identical(wlr_dhat_estimates(u, tzero = 18), wlr_dhat_estimates(c1, tzero = 18))
  got <- wlr_dhat_estimates(u, tzero = 18)
  dN <- diff(c(0, u$nbar0)) + diff(c(0, u$nbar1)); yb <- u$ybar0 + u$ybar1
  dJ <- ifelse(yb == 1, 0, (dN - 1) / (yb - 1)); dL <- ifelse(yb == 0, 0, dN / yb)
  h2 <- ifelse(u$ybar0 * u$ybar1 > 0, yb / (u$ybar0 * u$ybar1), 0) * (u$at_points <= 18)
  i <- max(which(u$at_points <= 18))
  expect_equal(got$sig2_dhat, u$survP[i]^2 * sum(h2 * (1 - dJ) * dL))
  expect_equal(got$cov_wlr_dhat, u$survP[i] * sum((u$at_points <= 18) * (1 - dJ) * dL))
})
