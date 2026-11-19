# R2 (Stage 2c): for weights that are not all equal, a resampled Kaplan-Meier deviation is
#   -S(t) * sum_i G_ib w_i psi_i(t)
# with psi_i the per-subject influence contribution of R1 and G the same multiplier draws as before.
# No weights / equal weights: the events-only perturbation, unchanged.

psi_mat <- function(U, D, W, at) {  # independent: n x length(at) matrix of psi_i(t)
  s <- sort(unique(U[D == 1]))
  Yw <- vapply(s, function(x) sum(W[U >= x]), 0); dNw <- vapply(s, function(x) sum(W[U == x & D == 1]), 0)
  ok <- Yw > 0 & (Yw - dNw) > 1e-10 * Yw
  sapply(at, function(t) vapply(seq_along(U), function(i) {
    first <- if (D[i] == 1 && U[i] <= t) { j <- match(U[i], s); if (ok[j]) 1 / (Yw[j] - dNw[j]) else 0 } else 0
    sel <- s <= min(t, U[i]) & ok
    first - sum(dNw[sel] / (Yw[sel] * (Yw[sel] - dNw[sel])))
  }, 0))
}
old_resampling <- function(U, W, D, at, surv, G) {  # the events-only perturbation (Stage 2b code)
  event_mat <- outer(U, at, FUN = "<="); risk_w <- colSums(outer(U, at, FUN = ">=") * W)
  dN <- apply(t(event_mat * W) %*% (D * G), 2, function(z) diff(c(0, z)))
  dr <- sweep(dN, 1, risk_w, "/"); dr[is.infinite(dr) | is.nan(dr)] <- 0
  (-1) * surv * apply(dr, 2, cumsum)
}
ipw_data <- function(N = 300, seed = 2) {
  old <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) get(".Random.seed", envir = globalenv()) else NULL
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else assign(".Random.seed", old, envir = globalenv()))
  set.seed(seed); x <- stats::rnorm(N); e <- stats::plogis(-0.5 + x); a <- stats::rbinom(N, 1, e)
  T <- stats::rexp(N, 0.08 * exp(0.7 * x - 0.4 * a)); cc <- stats::runif(N, 0, 40)
  data.frame(time = pmin(T, cc), status = as.numeric(T <= cc), treat = a, w = ifelse(a == 1, 1 / e, 1 / (1 - e)))
}
arm <- function(d, a) d[d$treat == a, ]

test_that("resampling_survival() with unequal weights is -S(t) * sum_i G_ib w_i psi_i(t), with and without ties", {
  for (ties in c(FALSE, TRUE)) {
    d <- ws_sim(120, 301, ties = ties, weighted = TRUE)
    x <- arm(d, 0); at <- sort(unique(c(x$time, 4.1)))
    re <- calculate_risk_event_counts(x$time, x$status, x$w, at)
    S <- KM_estimates(re$ybar, re$nbar, re$sig2w_multiplier)$S_KM
    G <- matrix(stats::rnorm(nrow(x) * 7), ncol = 7)
    got <- resampling_survival(x$time, x$w, x$status, at, 7, S, G)
    ref <- -S * t(psi_mat(x$time, x$status, x$w, at)) %*% (G * x$w)
    expect_equal(unname(got), unname(ref), tolerance = 1e-10)
  }
})

test_that("equal weights and no weights keep the events-only perturbation exactly", {
  d <- ws_sim(150, 302, ties = TRUE)
  x <- arm(d, 1); at <- sort(unique(x$time)); n <- nrow(x)
  G <- matrix(stats::rnorm(n * 5), ncol = 5)
  for (cc in c(1, 3)) {
    W <- rep(cc, n)
    re <- calculate_risk_event_counts(x$time, x$status, W, at)
    S <- KM_estimates(re$ybar, re$nbar, re$sig2w_multiplier)$S_KM
    expect_identical(resampling_survival(x$time, W, x$status, at, 5, S, G), old_resampling(x$time, W, x$status, at, S, G))
  }
  # KM_diff without weights: resampled curves are the old ones, drawn from the same seed
  f <- KM_diff(d, "time", "status", "treat", draws.band = 50, seedstart = 11)
  set.seed(11); G0 <- matrix(stats::rnorm(50 * sum(d$treat == 0)), ncol = 50)
  x0 <- arm(d, 0); x0 <- x0[order(x0$time), ]; g <- sort(unique(c(x0$time, f$at_points)))
  re <- calculate_risk_event_counts(x0$time, x0$status, rep(1, nrow(x0)), g)
  S <- KM_estimates(re$ybar, re$nbar, re$sig2w_multiplier)$S_KM
  expect_identical(f$surv0_star, old_resampling(x0$time, rep(1, nrow(x0)), x0$status, g, S, G0)[match(f$at_points, g), , drop = FALSE])
})

test_that("variance of the resampled weighted curves agrees with the analytic (robust) variance", {
  d <- ipw_data()
  f <- KM_diff(d, "time", "status", "treat", weight.name = "w", draws.band = 4000, show_resamples = FALSE)
  for (g in 0:1) {
    S <- f[[paste0("surv", g)]]; v <- f[[paste0("sig2_surv", g)]]; vs <- apply(f[[paste0("surv", g, "_star")]], 1, stats::var)
    k <- S > 0.1 & S < 0.9 & v > 0; r <- vs[k] / v[k]
    expect_gt(sum(k), 20); expect_gt(min(r), 0.85); expect_lt(max(r), 1.15)
  }
})

test_that("draws > 0 variance path: the SE agrees with the analytic SE, and is the variance over draws of the influence sums", {
  d <- ipw_data()
  f0 <- KM_diff(d, "time", "status", "treat", weight.name = "w")
  f4 <- KM_diff(d, "time", "status", "treat", weight.name = "w", draws = 4000)
  for (g in 0:1) {
    S <- f0[[paste0("surv", g)]]; k <- S > 0.1 & S < 0.9 & f0[[paste0("sig2_surv", g)]] > 0
    r <- sqrt(f4[[paste0("sig2_surv", g)]][k] / f0[[paste0("sig2_surv", g)]][k])
    expect_gt(min(r), 0.93); expect_lt(max(r), 1.07)
  }
  x <- arm(ws_sim(90, 303, weighted = TRUE), 1); at <- sort(unique(x$time))
  re <- calculate_risk_event_counts(x$time, x$status, x$w, at, draws = 30, seedstart = 5)
  set.seed(5); G <- matrix(stats::rnorm(30 * nrow(x)), ncol = 30)
  M <- t(psi_mat(x$time, x$status, x$w, at)) %*% (G * x$w)
  expect_equal(cumsum(re$sig2w_multiplier), apply(M, 1, stats::var), tolerance = 1e-9)
})

test_that("weighted resampling is reproducible and leaves the caller's RNG state alone", {
  d <- ipw_data(200, 4)
  set.seed(123); before <- .Random.seed
  a <- KM_diff(d, "time", "status", "treat", weight.name = "w", draws.band = 100)
  expect_identical(.Random.seed, before)
  b <- KM_diff(d, "time", "status", "treat", weight.name = "w", draws.band = 100)
  expect_identical(a, b)
  r1 <- cumulative_rmst_bands(d, a, "time", "status", "treat", weight.name = "w", draws_sb = 100, plot = FALSE)
  r2 <- cumulative_rmst_bands(d, a, "time", "status", "treat", weight.name = "w", draws_sb = 100, plot = FALSE)
  expect_identical(r1, r2)
  expect_identical(.Random.seed, before)
})
