# L1: log-rank with subject weights -- variance of the score = sum of squared per-subject score contributions.

dcw <- function(d, w = "w", ...) ws_quiet(df_counting(d, "time", "status", "treat", weight.name = w, ...))
rel <- function(a, b) abs(a / b - 1)

test_that("L1: FH(0,0) z^2 equals the robust score test of coxph (ties and strata included)", {
  for (ties in c(FALSE, TRUE)) {
    d <- ws_sim_ipw(300, 4, ties)
    cp <- survival::coxph(survival::Surv(time, status) ~ treat, data = d, weights = w, robust = TRUE, ties = "breslow")
    cps <- survival::coxph(survival::Surv(time, status) ~ treat + survival::strata(s), data = d, weights = w, robust = TRUE, ties = "breslow")
    dc <- dcw(d); dcs <- dcw(d, strata.name = "s")
    expect_lt(rel(dc$z.score^2, cp$rscore), 1e-8)
    expect_lt(rel(dcs$z.score_stratified^2, cps$rscore), 1e-8)
    expect_equal(dc$sig2_lr, dc$get_score$sig2.score)
    expect_equal(dc$z.score, dc$lr / sqrt(dc$sig2_lr))
    expect_equal(dcs$z.score_stratified, dcs$lr_stratified / sqrt(dcs$sig2_lr_stratified))
  }
})

test_that("L1: every scheme's sig2.score equals the reference sum of psi^2 with the scheme's time weights", {
  for (ties in c(FALSE, TRUE)) {
    d <- ws_sim_ipw(300, 14, ties)
    for (sc in ws_schemes(d)) {
      dc <- dcw(d, scheme = sc[[1]], scheme_params = sc[[2]])
      K <- ws_quiet(cox_rhogamma(dc, scheme = sc[[1]], scheme_params = sc[[2]], draws = 0))$fit$wt_rg   # time weight per observation
      ref <- ws_score_resid_ref(dc$time, dc$delta, dc$z, dc$w_hat, K = K, beta = 0)
      expect_lt(rel(dc$sig2_lr, sum(ref$psi^2)), 1e-8)
      expect_lt(rel(dc$lr, -ref$U), 1e-8)                 # the score is the one the package already computes
      # a df_counting() result passed directly to wlr_cumulative() gives the same
      g <- wlr_cumulative(dc, sc[[1]], sc[[2]], return_cumulative = FALSE)
      expect_equal(g$sig2.score, dc$sig2_lr, tolerance = 1e-12)
      expect_equal(g$score, dc$lr, tolerance = 1e-12)
    }
  }
})

test_that("L1: cumulative variance equals the direct subject-by-time computation; last value = non-cumulative", {
  for (ties in c(FALSE, TRUE)) {
    d <- ws_sim_ipw(150, 24, ties)
    dc <- dcw(d, scheme = "fh", scheme_params = list(rho = 0, gamma = 1))
    cum <- wlr_cumulative(dc, "fh", list(rho = 0, gamma = 1), return_cumulative = TRUE)
    one <- wlr_cumulative(dc, "fh", list(rho = 0, gamma = 1), return_cumulative = FALSE)
    U <- dc$time; D <- dc$delta; Z <- dc$z; a <- dc$w_hat; g <- dc$at_points
    Kobs <- ws_quiet(cox_rhogamma(dc, scheme = "fh", scheme_params = list(rho = 0, gamma = 1), draws = 0))$fit$wt_rg
    s <- sort(unique(U))
    Ks <- Kobs[match(s, U)]
    S0 <- vapply(s, function(x) sum(a[U >= x]), 0); zb <- vapply(s, function(x) sum(a[U >= x & Z == 1]), 0) / S0
    dL <- vapply(s, function(x) sum(a[U == x & D == 1]), 0) / S0
    V <- vapply(g, function(t) {
      ci <- vapply(seq_along(U), function(i) {
        ev <- if (D[i] == 1 && U[i] <= t) Ks[match(U[i], s)] * (Z[i] - zb[match(U[i], s)]) else 0
        k <- s <= min(t, U[i])
        a[i] * (ev - sum(Ks[k] * (Z[i] - zb[k]) * dL[k]))
      }, 0)
      sum(ci^2)
    }, 0)
    pos <- V > 0
    expect_lt(max(abs(cum$sig2.score[pos] / V[pos] - 1)), 1e-10)
    expect_equal(cum$sig2.score[!pos], V[!pos])
    expect_lt(rel(cum$sig2.score[length(g)], one$sig2.score), 1e-10)
    expect_equal(cum$score[length(g)], one$score, tolerance = 1e-12)
    expect_equal(cum$z.score, cum$score / sqrt(cum$sig2.score))
  }
})

test_that("L1: z.score does not depend on the scale of the weights and is finite for normalised weights", {
  for (ties in c(FALSE, TRUE)) {
    d <- ws_sim_ipw(300, 34, ties); d$w10 <- 10 * d$w; d$wn <- d$w / sum(d$w)
    for (sc in ws_schemes(d)) {
      z <- dcw(d, scheme = sc[[1]], scheme_params = sc[[2]])$z.score
      z10 <- dcw(d, "w10", scheme = sc[[1]], scheme_params = sc[[2]])$z.score
      zn <- dcw(d, "wn", scheme = sc[[1]], scheme_params = sc[[2]])$z.score
      expect_true(all(is.finite(c(z, z10, zn))))
      expect_lt(rel(z10, z), 1e-9); expect_lt(rel(zn, z), 1e-9)
    }
    zs <- dcw(d, strata.name = "s")$z.score_stratified
    expect_lt(rel(dcw(d, "w10", strata.name = "s")$z.score_stratified, zs), 1e-9)
    expect_lt(rel(dcw(d, "wn", strata.name = "s")$z.score_stratified, zs), 1e-9)
  }
})

test_that("L1: no weights, equal weights and bare counts keep the earlier formula", {
  old_sig2 <- function(dc, wt = 1, cumulative = FALSE) {         # the formula of wlr_cumulative() before Stage 2d, FH(0,0)
    dN <- diff(c(0, dc$nbar0)) + diff(c(0, dc$nbar1)); yb <- dc$ybar0 + dc$ybar1
    K <- ifelse(yb > 0, wt * (dc$ybar0 * dc$ybar1) / yb, 0)
    h0 <- ifelse(dc$ybar0 == 0, 0, K^2 / dc$ybar0); h1 <- ifelse(dc$ybar1 == 0, 0, K^2 / dc$ybar1)
    dJ <- ifelse(yb == 1, 0, (dN - 1) / (yb - 1)); dL <- ifelse(yb == 0, 0, dN / yb)
    v <- (h0 + h1) * (1 - dJ) * dL
    if (cumulative) cumsum(v) else sum(v)
  }
  for (ties in c(FALSE, TRUE)) {
    d <- ws_sim_ipw(200, 44, ties)
    u <- ws_quiet(df_counting(d, "time", "status", "treat")); c1 <- dcw(d, "one")
    expect_identical(u$z.score, c1$z.score); expect_identical(u$sig2_lr, c1$sig2_lr)
    expect_identical(u$sig2_lr, old_sig2(u))
    expect_identical(wlr_cumulative(u, "fh", return_cumulative = TRUE)$sig2.score, old_sig2(u, cumulative = TRUE))
    us <- ws_quiet(df_counting(d, "time", "status", "treat", strata.name = "s")); cs <- dcw(d, "one", strata.name = "s")
    expect_identical(us$z.score_stratified, cs$z.score_stratified); expect_identical(us$sig2_lr_stratified, cs$sig2_lr_stratified)
    # a bare table of counts taken from a weighted analysis: no subject-level data, earlier formula
    dc <- dcw(d)
    bare <- data.frame(at_points = dc$at_points, nbar0 = dc$nbar0, nbar1 = dc$nbar1, ybar0 = dc$ybar0, ybar1 = dc$ybar1, survP = dc$survP)
    expect_identical(wlr_cumulative(bare, "fh")$sig2.score, old_sig2(dc))
    expect_identical(wlr_cumulative(bare, "fh", return_cumulative = TRUE)$sig2.score, old_sig2(dc, cumulative = TRUE))
    expect_false(isTRUE(all.equal(wlr_cumulative(bare, "fh")$sig2.score, dc$sig2_lr)))
    expect_identical(wlr_cumulative(bare, "fh")$score, dc$lr)        # the score is the same either way
  }
})

test_that("L1: zero-weight rows contribute nothing to the weighted log-rank", {
  d <- ws_sim_ipw(200, 54); d0 <- d; d0$w[c(3, 50, 120)] <- 0
  full <- dcw(d0); sub <- dcw(d0[d0$w > 0, ])
  expect_equal(full$lr, sub$lr, tolerance = 1e-10)
  expect_equal(full$sig2_lr, sub$sig2_lr, tolerance = 1e-10)
})
