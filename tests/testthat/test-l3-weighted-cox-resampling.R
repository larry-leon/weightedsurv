# L3: cox_rhogamma() resampling with subject weights -- resampled scores from the per-subject contributions.

dcw <- function(d, w = "w", ...) ws_quiet(df_counting(d, "time", "status", "treat", weight.name = w, ...))

test_that("L3: resampled scores are -sum_i G_ib psi_i, with the multipliers generated as before", {
  for (ties in c(FALSE, TRUE)) {
    d <- ws_sim_ipw(200, 61, ties)
    dc <- dcw(d, scheme = "fh", scheme_params = list(rho = 0, gamma = 1))
    B <- 50
    cx <- ws_quiet(cox_rhogamma(dc, scheme = "fh", scheme_params = list(rho = 0, gamma = 1), draws = B))
    f <- cx$fit; n0 <- sum(f$z == 0); n1 <- sum(f$z == 1)
    set.seed(8316951)                                   # seedstart of cox_rhogamma_resample()
    G0 <- matrix(stats::rnorm(B * n0), ncol = B); G1 <- matrix(stats::rnorm(B * n1), ncol = B)
    G <- matrix(0, length(f$z), B); G[f$z == 0, ] <- G0; G[f$z == 1, ] <- G1
    pb <- ws_score_resid_ref(f$time, f$delta, f$z, f$w_hat, K = f$wt_rg, beta = f$bhat)
    p0 <- ws_score_resid_ref(f$time, f$delta, f$z, f$w_hat, K = f$wt_rg, beta = 0)
    rs <- cx$fit_resamples
    expect_equal(rs$score_star, -colSums(G * pb$psi), tolerance = 1e-10)
    expect_equal(rs$score_star_null, -colSums(G * p0$psi), tolerance = 1e-10)
    # everything computed from them is formed as before
    expect_equal(rs$bhat_center_star, rs$score_star / pb$info, tolerance = 1e-8)
    expect_equal(rs$Z_bstar, rs$bhat_center_star / f$sig_bhat_asy)
    expect_equal(rs$sig_bhat_star, stats::sd(rs$bhat_center_star))
    expect_equal(f$sig_bhat_star, rs$sig_bhat_star)
    expect_equal(f$bhat_debiased, f$bhat - mean(rs$bhat_center_star))
    expect_equal(cx$sig2_score_star, stats::var(rs$score_star))
    expect_equal(cx$z.score_debiased, (f$u.zero - mean(rs$score_star_null)) / sqrt(f$sig2_score))
    expect_equal(cx$hr_ci_star$sig_bhat, rs$sig_bhat_star)
  }
})

test_that("L3: with many draws the resampling SE agrees with the analytic SE (weighted)", {
  d <- ws_sim_ipw(300, 5)
  cx <- ws_quiet(cox_rhogamma(dcw(d), scheme = "fh", scheme_params = list(rho = 0, gamma = 0), draws = 4000))
  r <- cx$fit$sig_bhat_star / cx$fit$sig_bhat_asy
  expect_gt(r, 0.95); expect_lt(r, 1.05)
  v <- stats::var(cx$fit_resamples$score_star_null) / cx$fit$sig2_score
  expect_gt(v, 0.90); expect_lt(v, 1.10)
})

test_that("L3: same seed, same output; the caller's RNG state is left as it was", {
  d <- ws_sim_ipw(150, 62)
  dc <- dcw(d)
  set.seed(99); before <- .Random.seed
  a <- ws_quiet(cox_rhogamma(dc, scheme = "fh", scheme_params = list(rho = 0, gamma = 0), draws = 100))
  expect_identical(.Random.seed, before)
  b <- ws_quiet(cox_rhogamma(dc, scheme = "fh", scheme_params = list(rho = 0, gamma = 0), draws = 100))
  expect_identical(a, b)
})

test_that("L3: no weights and equal weights resample as before", {
  d <- ws_sim_ipw(150, 63, ties = TRUE)
  u <- ws_quiet(cox_rhogamma(ws_quiet(df_counting(d, "time", "status", "treat")), draws = 100))
  c1 <- ws_quiet(cox_rhogamma(dcw(d, "one"), draws = 100))
  expect_identical(u, c1)
  # the earlier resampling perturbs the events only: one resampled score from the counts
  f <- u$fit; B <- 100; n0 <- sum(f$z == 0); n1 <- sum(f$z == 1)
  set.seed(8316951)
  G0 <- matrix(stats::rnorm(B * n0), ncol = B); G1 <- matrix(stats::rnorm(B * n1), ncol = B)
  t0 <- f$time[f$z == 0]; t1 <- f$time[f$z == 1]; d0 <- f$delta[f$z == 0]; d1 <- f$delta[f$z == 1]
  y0 <- colSums(outer(t0, f$time, ">=")); y1 <- colSums(outer(t1, f$time, ">=")) * exp(f$bhat)
  K <- ifelse(y0 + y1 > 0, f$wt_rg * y1 * y0 / (y0 + y1), 0)
  dn0 <- diff(c(0, colSums(outer(t0, f$time, "<=") * (d0 * G0[, 1])))); dn1 <- diff(c(0, colSums(outer(t1, f$time, "<=") * (d1 * G1[, 1]))))
  s1 <- sum(K * (ifelse(y0 > 0, dn0 / y0, 0) - ifelse(y1 > 0, dn1 / y1, 0)))
  expect_equal(u$fit_resamples$score_star[1], s1, tolerance = 1e-10)
})
