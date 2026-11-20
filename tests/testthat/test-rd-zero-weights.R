# Rd: the test of "weights not all equal" looks at the positive weights; zero-weight rows are inert.

rd_data <- function(ties) {
  old <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) get(".Random.seed", envir = globalenv()) else NULL
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else assign(".Random.seed", old, envir = globalenv()))
  set.seed(3); n <- 300; x <- rnorm(n); a <- rbinom(n, 1, plogis(-0.3 + x))
  tt <- rexp(n, 0.05 * exp(0.6 * x - 0.4 * a)); cc <- runif(n, 0, 60); keep <- runif(n) < 0.7
  d <- data.frame(time = if (ties) ceiling(pmin(tt, cc)) else pmin(tt, cc), status = as.numeric(tt <= cc), treat = a,
                  s = ifelse(x > 0, "hi", "lo"), stringsAsFactors = FALSE)
  d$w01 <- as.numeric(keep)
  d$wu <- ifelse(keep, ifelse(a == 1, 1 / plogis(-0.3 + x), 1 / (1 - plogis(-0.3 + x))), 0)
  list(d = d, sub = d[keep, ])
}
rel <- function(a, b) max(abs(a - b) / pmax(abs(b), 1e-300))

test_that("Rd: the helper decides on the positive weights", {
  expect_false(ws_weights_unequal(c(0, 1, 1, 0)))
  expect_false(ws_weights_unequal(c(2, 2, 2)))
  expect_false(ws_weights_unequal(numeric(0)))
  expect_true(ws_weights_unequal(c(0, 1, 2)))
  expect_true(ws_weights_unequal(c(1, 2)))
})

test_that("Rd: zero-weight rows are inert (0/1 weights and zeros among unequal weights; with and without ties)", {
  for (ties in c(FALSE, TRUE)) {
    D <- rd_data(ties); d <- D$d; sub <- D$sub
    for (wn in c("w01", "wu")) {
      wsub <- if (wn == "w01") NULL else wn
      f <- ws_quiet(df_counting(d, "time", "status", "treat", weight.name = wn, strata.name = "s"))
      p <- ws_quiet(df_counting(sub, "time", "status", "treat", weight.name = wsub, strata.name = "s"))
      j <- match(p$at_points, f$at_points)
      expect_false(anyNA(j))
      expect_lt(max(abs(f$sig2_surv1[j] - p$sig2_surv1)), 1e-12)
      expect_lt(max(abs(f$sig2_surv0[j] - p$sig2_surv0)), 1e-12)
      expect_lt(max(abs(f$sig2_survP[j] - p$sig2_survP)), 1e-12)
      expect_lt(abs(f$lr / p$lr - 1), 1e-10)
      expect_lt(abs(f$sig2_lr / p$sig2_lr - 1), 1e-10)
      expect_lt(abs(f$z.score / p$z.score - 1), 1e-10)
      expect_lt(abs(f$sig2_lr_stratified / p$sig2_lr_stratified - 1), 1e-10)
      expect_lt(abs(f$z.score_stratified / p$z.score_stratified - 1), 1e-10)
      for (sc in ws_schemes(sub)[c("fh00", "fh01")]) {
        cf <- ws_quiet(cox_rhogamma(f, scheme = sc[[1]], scheme_params = sc[[2]]))
        cp <- ws_quiet(cox_rhogamma(p, scheme = sc[[1]], scheme_params = sc[[2]]))
        expect_lt(abs(cf$fit$bhat - cp$fit$bhat), 1e-8)
        expect_lt(abs(cf$fit$sig_bhat_asy / cp$fit$sig_bhat_asy - 1), 1e-8)
        expect_lt(abs(cf$z.score / cp$z.score - 1), 1e-10)
      }
      kf <- ws_quiet(KM_diff(d, "time", "status", "treat", weight.name = wn, at_points = p$at_points))
      kp <- ws_quiet(KM_diff(sub, "time", "status", "treat", weight.name = wsub, at_points = p$at_points))
      expect_lt(max(abs(kf$sig2_dhat - kp$sig2_dhat)), 1e-12)
      expect_identical(kf$dhat, kp$dhat)
    }
  }
})

test_that("Rd: resampled quantities with zero weights estimate the same variance (not equal draw by draw)", {
  D <- rd_data(FALSE); d <- D$d; sub <- D$sub
  kf <- KM_diff(d, "time", "status", "treat", weight.name = "w01", at_points = sort(sub$time), draws.band = 2000)
  kp <- KM_diff(sub, "time", "status", "treat", at_points = sort(sub$time), draws.band = 2000)
  expect_identical(dim(kf$dhat_star), dim(kp$dhat_star))
  expect_lt(abs(kf$c_alpha_band / kp$c_alpha_band - 1), 0.1)
  f <- ws_quiet(df_counting(d, "time", "status", "treat", weight.name = "w01"))
  p <- ws_quiet(df_counting(sub, "time", "status", "treat"))
  cf <- ws_quiet(cox_rhogamma(f, scheme = "fh", scheme_params = list(rho = 0, gamma = 0), draws = 2000))
  cp <- ws_quiet(cox_rhogamma(p, scheme = "fh", scheme_params = list(rho = 0, gamma = 0), draws = 2000))
  expect_lt(abs(cf$fit$sig_bhat_star / cp$fit$sig_bhat_star - 1), 0.1)
})

test_that("Rd: without a zero weight nothing changes (unequal weights still take the robust formulas)", {
  d <- ws_sim_ipw(200, 5)
  f <- ws_quiet(df_counting(d, "time", "status", "treat", weight.name = "w"))
  cp <- survival::coxph(survival::Surv(time, status) ~ treat, data = d, weights = w, robust = TRUE, ties = "breslow")
  expect_lt(abs(f$z.score^2 / cp$rscore - 1), 1e-8)
})
