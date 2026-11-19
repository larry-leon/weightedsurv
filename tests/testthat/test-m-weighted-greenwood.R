# M: the analytic weighted Kaplan-Meier variance is weighted Greenwood,
#   Var[S(t)] = S(t)^2 * sum_{s <= t} [sum_i w_i^2 dN_i(s)] / [Y_w(s) (Y_w(s) - dN_w(s))].
# (No test compares a weighted SE with survfit.)

wg_ref <- function(U, D, W, at) {  # independent: per time point, no cumulative sums
  vapply(at, function(t) {
    s <- sort(unique(U[D == 1 & U <= t]))
    if (!length(s)) return(0)
    sum(vapply(s, function(u) {
      Yw <- sum(W[U >= u]); dNw <- sum(W[U == u & D == 1]); d2 <- sum(W[U == u & D == 1]^2)
      if (Yw > 0 && Yw - dNw > 1e-10 * Yw) d2 / (Yw * (Yw - dNw)) else 0
    }, 0))
  }, 0)
}
mult <- function(U, D, W, at, ...) calculate_risk_event_counts(U, D, W, at, ...)$sig2w_multiplier
arm <- function(d, a) d[d$treat == a, ]

test_that("unit weights: Greenwood, identical to survfit's standard error", {
  for (ties in c(FALSE, TRUE)) {
    d <- arm(ws_sim(200, 101, ties = ties), 1)
    at <- sort(unique(d$time))
    re <- calculate_risk_event_counts(d$time, d$status, rep(1, nrow(d)), at)
    km <- KM_estimates(re$ybar, re$nbar, re$sig2w_multiplier)
    dN <- diff(c(0, re$nbar))
    expect_equal(re$sig2w_multiplier, ifelse(re$ybar > dN, dN / (re$ybar * (re$ybar - dN)), 0), tolerance = 1e-12)
    fa <- survival::survfit(survival::Surv(time, status) ~ 1, data = d)
    i <- match(fa$time, at); ok <- fa$surv > 0
    expect_equal(sqrt(km$sig2_KM[i][ok]), (fa$surv * fa$std.err)[ok], tolerance = 1e-12)
  }
})

test_that("the multiplier equals the formula computed independently, with and without ties", {
  for (ties in c(FALSE, TRUE)) {
    d <- ws_sim(250, 102, ties = ties, weighted = TRUE)
    for (a in 0:1) {
      x <- arm(d, a); at <- sort(unique(c(0, x$time, 7.25)))
      expect_equal(cumsum(mult(x$time, x$status, x$w, at)), wg_ref(x$time, x$status, x$w, at), tolerance = 1e-9)
    }
  }
})

test_that("the SE does not depend on the scale of the weights", {
  for (ties in c(FALSE, TRUE)) {
    d <- ws_sim(200, 103, ties = ties, weighted = TRUE)
    x <- arm(d, 0); at <- sort(unique(x$time)); w <- x$w
    f <- function(w) mult(x$time, x$status, w, at)
    for (w2 in list(10 * w, 0.01 * w, w / sum(w))) expect_equal(f(w2), f(w), tolerance = 1e-9)
    se <- function(w) { dd <- d; dd$w[dd$treat == 0] <- w; k <- KM_diff(dd, "time", "status", "treat", weight.name = "w"); sqrt(k$sig2_surv0) }
    expect_equal(se(10 * w), se(w), tolerance = 1e-9)
    expect_equal(se(w / sum(w)), se(w), tolerance = 1e-9)
  }
})

test_that("constant weights give the unit-weight SE", {
  d <- ws_sim(180, 104, ties = TRUE)
  x <- arm(d, 1); at <- sort(unique(x$time)); n <- nrow(x)
  for (cc in c(2, 0.37, 50)) expect_equal(mult(x$time, x$status, rep(cc, n), at), mult(x$time, x$status, rep(1, n), at), tolerance = 1e-12)
  d$wc <- 3
  a <- df_counting(d, "time", "status", "treat", weight.name = "wc")
  b <- df_counting(d, "time", "status", "treat")
  expect_equal(a$sig2_surv0, b$sig2_surv0, tolerance = 1e-12)
  expect_equal(a$sig2_surv1, b$sig2_surv1, tolerance = 1e-12)
})

test_that("df_counting() and KM_diff() report the formula; the survival estimate is unchanged", {
  d <- ws_sim(220, 105, ties = TRUE, weighted = TRUE)
  dc <- df_counting(d, "time", "status", "treat", weight.name = "w")
  kd <- KM_diff(d, "time", "status", "treat", weight.name = "w")
  for (a in 0:1) {
    x <- arm(d, a)
    S <- dc[[paste0("surv", a)]]
    expect_equal(dc[[paste0("sig2_surv", a)]], S^2 * wg_ref(x$time, x$status, x$w, dc$at_points), tolerance = 1e-9)
    Sk <- kd[[paste0("surv", a)]]
    expect_equal(kd[[paste0("sig2_surv", a)]], Sk^2 * wg_ref(x$time, x$status, x$w, kd$at_points), tolerance = 1e-9)
    fa <- survival::survfit(survival::Surv(time, status) ~ 1, data = x, weights = w)
    expect_equal(Sk, summary(fa, times = kd$at_points, extend = TRUE)$surv, tolerance = 1e-10)
  }
})

test_that("resampling path (draws > 0, unequal weights) is the one used before and is unchanged", {
  d <- ws_sim(150, 106, weighted = TRUE)
  x <- arm(d, 1); at <- sort(unique(x$time)); U <- x$time; D <- x$status; W <- x$w
  got <- mult(U, D, W, at, draws = 100, seedstart = 77)
  old <- local({
    event_mat <- outer(U, at, FUN = "<="); risk_w <- colSums(outer(U, at, FUN = ">=") * W)
    set.seed(77); G <- matrix(stats::rnorm(100 * length(U)), ncol = 100)
    dN <- apply(t(event_mat * W) %*% (D * G), 2, function(z) diff(c(0, z)))
    dr <- sweep(dN, 1, risk_w, "/"); dr[is.infinite(dr) | is.nan(dr)] <- 0
    apply(dr, 1, stats::var)
  })
  expect_identical(got, old)
  # equal weights with draws > 0 take the formula, as before
  expect_identical(mult(U, D, rep(1, length(U)), at, draws = 100), mult(U, D, rep(1, length(U)), at))
})

test_that("an empty group gives zeros", {
  expect_identical(mult(numeric(0), numeric(0), numeric(0), c(0, 1, 2)), c(0, 0, 0))
})
