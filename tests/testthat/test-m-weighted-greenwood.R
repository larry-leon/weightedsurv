# M (Stage 2b): unit and constant weights give Greenwood; the analytic weighted SE is scale-free.
# Stage 2c (R1) replaced the two blocks of this file that encoded the weighted Greenwood formula
# ("the multiplier equals the formula computed independently ..." and "df_counting() and KM_diff()
# report the formula ..."); the weighted variance is now the robust variance, tested in
# test-r1-robust-variance.R. Stage 2c (R2) replaced the block that required the weighted
# resampling path (draws > 0, unequal weights) to be unchanged; its equal-weights half is kept,
# and the new resampling is tested in test-r2-weighted-resampling.R.
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

test_that("draws > 0 with equal weights takes the formula, as before", {
  d <- ws_sim(150, 106, weighted = TRUE)
  x <- arm(d, 1); at <- sort(unique(x$time)); U <- x$time; D <- x$status
  expect_identical(mult(U, D, rep(1, length(U)), at, draws = 100), mult(U, D, rep(1, length(U)), at))
})

test_that("an empty group gives zeros", {
  expect_identical(mult(numeric(0), numeric(0), numeric(0), c(0, 1, 2)), c(0, 0, 0))
})
