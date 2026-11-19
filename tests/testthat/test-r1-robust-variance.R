# R1 (Stage 2c): for weights that are not all equal, the analytic variance of a weighted
# Kaplan-Meier curve is the robust (infinitesimal-jackknife) variance,
#   Var[S(t)] = S(t)^2 * sum_i w_i^2 psi_i(t)^2,
#   psi_i(t) = D_i 1(U_i <= t) / (Y_w(U_i) - dN_w(U_i)) - sum_{event s <= min(t, U_i)} dN_w(s) / (Y_w(s) (Y_w(s) - dN_w(s))),
# which is survfit(weights = w, robust = TRUE)'s variance. No weights / equal weights: Greenwood, unchanged.

ij_ref <- function(U, D, W, at) {  # independent: per time point and per subject, no cumulative sums
  s <- sort(unique(U[D == 1]))
  Yw <- vapply(s, function(x) sum(W[U >= x]), 0); dNw <- vapply(s, function(x) sum(W[U == x & D == 1]), 0)
  ok <- Yw > 0 & (Yw - dNw) > 1e-10 * Yw
  vapply(at, function(t) {
    psi <- vapply(seq_along(U), function(i) {
      first <- if (D[i] == 1 && U[i] <= t) { j <- match(U[i], s); if (ok[j]) 1 / (Yw[j] - dNw[j]) else 0 } else 0
      sel <- s <= min(t, U[i]) & ok
      first - sum(dNw[sel] / (Yw[sel] * (Yw[sel] - dNw[sel])))
    }, 0)
    sum(W^2 * psi^2)
  }, 0)
}
robust_se <- function(U, D, W, at, timefix = TRUE) {
  y <- if (timefix) survival::aeqSurv(survival::Surv(U, D)) else survival::Surv(U, D)
  summary(survival::survfit(y ~ 1, weights = W, robust = TRUE, timefix = FALSE), times = at, extend = TRUE)$std.err
}
rel_ok <- function(got, ref, tol) { k <- is.finite(ref) & ref > 0; expect_gt(sum(k), 3); expect_lt(max(abs(got[k] / ref[k] - 1)), tol) }
arm <- function(d, a) d[d$treat == a, ]

test_that("the multiplier cumulates to the per-subject formula computed independently, with and without ties", {
  for (ties in c(FALSE, TRUE)) {
    d <- ws_sim(160, 202, ties = ties, weighted = TRUE)
    for (a in 0:1) {
      x <- arm(d, a); at <- sort(unique(c(0, x$time, 7.25)))
      re <- calculate_risk_event_counts(x$time, x$status, x$w, at)
      expect_equal(cumsum(re$sig2w_multiplier), ij_ref(x$time, x$status, x$w, at), tolerance = 1e-9)
    }
  }
})

test_that("calculate_risk_event_counts() gives survfit robust's SE (timefix = FALSE), with and without ties", {
  for (ties in c(FALSE, TRUE)) {
    d <- ws_sim(300, 203, ties = ties, weighted = TRUE)
    x <- arm(d, 1); at <- sort(unique(x$time))
    re <- calculate_risk_event_counts(x$time, x$status, x$w, at)
    km <- KM_estimates(re$ybar, re$nbar, re$sig2w_multiplier)
    rel_ok(sqrt(pmax(km$sig2_KM, 0)), robust_se(x$time, x$status, x$w, at, timefix = FALSE), 1e-8)
  }
})

test_that("KM_diff() and df_counting() weighted SEs equal survfit robust on the time-fixed data; estimates unchanged", {
  for (ties in c(FALSE, TRUE)) {
    d <- ws_sim(260, 204, ties = ties, weighted = TRUE)
    kd <- KM_diff(d, "time", "status", "treat", weight.name = "w")
    dc <- df_counting(d, "time", "status", "treat", weight.name = "w")
    for (a in 0:1) {
      x <- arm(d, a)
      rel_ok(sqrt(kd[[paste0("sig2_surv", a)]]), robust_se(x$time, x$status, x$w, kd$at_points), 1e-8)
      rel_ok(sqrt(dc[[paste0("sig2_surv", a)]]), robust_se(x$time, x$status, x$w, dc$at_points), 1e-8)
      fa <- survival::survfit(survival::Surv(time, status) ~ 1, data = x, weights = w)
      expect_equal(kd[[paste0("surv", a)]], summary(fa, times = kd$at_points, extend = TRUE)$surv, tolerance = 1e-10)
    }
    rel_ok(sqrt(dc$sig2_survP), robust_se(d$time, d$status, d$w, dc$at_points), 1e-8)
  }
})

test_that("KM_diff() at arbitrary at_points (not event times) still gives survfit robust", {
  d <- ws_sim(240, 205, weighted = TRUE)
  at <- c(0.5, 2, 3.3, 7.25, 12, 20)
  kd <- KM_diff(d, "time", "status", "treat", weight.name = "w", at_points = at)
  for (a in 0:1) { x <- arm(d, a); rel_ok(sqrt(kd[[paste0("sig2_surv", a)]]), robust_se(x$time, x$status, x$w, at), 1e-8) }
})

test_that("the robust SE does not depend on the scale of the weights", {
  d <- ws_sim(220, 206, ties = TRUE, weighted = TRUE)
  se <- function(w) { dd <- d; dd$w <- w; k <- KM_diff(dd, "time", "status", "treat", weight.name = "w"); sqrt(c(k$sig2_surv0, k$sig2_surv1)) }
  ref <- se(d$w); k <- ref > 0
  for (w2 in list(10 * d$w, 0.01 * d$w, d$w / sum(d$w))) expect_lt(max(abs(se(w2)[k] / ref[k] - 1)), 1e-9)
})

test_that("equal weights (and no weights) keep the Greenwood code path exactly", {
  d <- ws_sim(200, 207, ties = TRUE)
  x <- arm(d, 0); at <- sort(unique(x$time)); n <- nrow(x)
  for (cc in c(1, 2.5)) {
    W <- rep(cc, n)
    re <- calculate_risk_event_counts(x$time, x$status, W, at)
    dN <- diff(c(0, re$nbar))
    dN2 <- diff(c(0, colSums(outer(x$time[x$status == 1], at, FUN = "<=") * W[x$status == 1]^2)))
    expect_identical(re$sig2w_multiplier, ifelse(re$ybar > 0 & (re$ybar - dN) > 1e-10 * re$ybar, dN2 / (re$ybar * (re$ybar - dN)), 0.0))
  }
  # at unit weights the robust formula itself equals Greenwood to rounding
  v_ij <- ij_ref(x$time, x$status, rep(1, n), at)
  re <- calculate_risk_event_counts(x$time, x$status, rep(1, n), at)
  expect_equal(cumsum(re$sig2w_multiplier), v_ij, tolerance = 1e-9)
})

test_that("unsorted at_points give the same variance at each point; zero weights contribute nothing", {
  d <- ws_sim(150, 208, weighted = TRUE)
  x <- arm(d, 1); at <- sort(unique(x$time))
  v <- function(at) { re <- calculate_risk_event_counts(x$time, x$status, x$w, at); cumsum(re$sig2w_multiplier) }
  p <- rev(seq_along(at))
  expect_equal(v(at[p]), v(at)[p], tolerance = 1e-12)
  W0 <- x$w; W0[1:5] <- 0
  re0 <- calculate_risk_event_counts(x$time, x$status, W0, at)
  re1 <- calculate_risk_event_counts(x$time[-(1:5)], x$status[-(1:5)], x$w[-(1:5)], at)
  expect_equal(cumsum(re0$sig2w_multiplier), cumsum(re1$sig2w_multiplier), tolerance = 1e-12)
})
