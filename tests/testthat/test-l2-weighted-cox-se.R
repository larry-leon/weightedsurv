# L2: cox_rhogamma() with subject weights -- SE and interval from the per-subject score contributions.

dcw <- function(d, w = "w", ...) ws_quiet(df_counting(d, "time", "status", "treat", weight.name = w, ...))
cxr <- function(dc, sc, draws = 0) ws_quiet(cox_rhogamma(dc, scheme = sc[[1]], scheme_params = sc[[2]], draws = draws))
rel <- function(a, b) abs(a / b - 1)

test_that("L2: FH(0,0) SE and z equal coxph(robust = TRUE, ties = 'breslow'), with and without ties", {
  for (ties in c(FALSE, TRUE)) {
    d <- ws_sim_ipw(300, 4, ties)
    cp <- survival::coxph(survival::Surv(time, status) ~ treat, data = d, weights = w, robust = TRUE, ties = "breslow")
    cx <- cxr(dcw(d), list("fh", list(rho = 0, gamma = 0)))
    expect_lt(rel(cx$fit$sig_bhat_asy, sqrt(cp$var[1, 1])), 1e-6)
    expect_lt(rel(cx$fit$z.score^2, cp$rscore), 1e-8)
    expect_lt(rel(cx$fit$bhat, unname(stats::coef(cp))), 1e-6)
    expect_equal(cx$hr_ci_asy$sig_bhat, cx$fit$sig_bhat_asy)
    expect_equal(log(cx$hr_ci_asy$upper), cx$fit$bhat + stats::qnorm(0.975) * cx$fit$sig_bhat_asy)
    expect_identical(cx$z.score, cx$fit$z.score)
  }
})

test_that("L2: every scheme's SE equals sqrt(sum psi^2) / info from the reference at bhat; z agrees with df_counting()", {
  for (ties in c(FALSE, TRUE)) {
    d <- ws_sim_ipw(300, 14, ties)
    for (sc in ws_schemes(d)) {
      dc <- dcw(d, scheme = sc[[1]], scheme_params = sc[[2]])
      cx <- cxr(dc, sc); f <- cx$fit
      ref <- ws_score_resid_ref(f$time, f$delta, f$z, f$w_hat, K = f$wt_rg, beta = f$bhat)
      ref0 <- ws_score_resid_ref(f$time, f$delta, f$z, f$w_hat, K = f$wt_rg, beta = 0)
      expect_lt(rel(f$sig_bhat_asy, sqrt(sum(ref$psi^2)) / ref$info), 1e-8)
      expect_lt(rel(f$sig2_score, sum(ref0$psi^2)), 1e-8)
      expect_lt(rel(f$u.zero, -ref0$U), 1e-8)
      expect_lt(rel(f$z.score, dc$z.score), 1e-8)
      # cox_score_rhogamma(score_only = FALSE) returns the same quantities
      cs <- cox_score_rhogamma(f$bhat, f$time, f$delta, f$z, w_hat = f$w_hat, wt_rg = f$wt_rg, score_only = FALSE)
      expect_lt(rel(cs$sig2_beta_asy, f$sig_bhat_asy^2), 1e-12)
      expect_lt(rel(cs$i_bhat, ref$info), 1e-8)
      expect_identical(cs$score, cox_score_rhogamma(f$bhat, f$time, f$delta, f$z, w_hat = f$w_hat, wt_rg = f$wt_rg))
    }
  }
})

test_that("L2: SE, z and the estimate do not depend on the scale of the weights", {
  for (ties in c(FALSE, TRUE)) {
    d <- ws_sim_ipw(300, 34, ties); d$w10 <- 10 * d$w; d$wn <- d$w / sum(d$w)
    for (sc in ws_schemes(d)) {
      a <- cxr(dcw(d, scheme = sc[[1]], scheme_params = sc[[2]]), sc)$fit
      for (wn in c("w10", "wn")) {
        b <- cxr(dcw(d, wn, scheme = sc[[1]], scheme_params = sc[[2]]), sc)$fit
        expect_true(all(is.finite(c(b$sig_bhat_asy, b$z.score))))
        expect_lt(rel(b$sig_bhat_asy, a$sig_bhat_asy), 1e-6)
        expect_lt(rel(b$z.score, a$z.score), 1e-9)
        expect_lt(abs(b$bhat - a$bhat), 1e-6)
      }
    }
  }
})

test_that("L2: no weights and equal weights keep the earlier variance", {
  for (ties in c(FALSE, TRUE)) {
    d <- ws_sim_ipw(200, 44, ties)
    for (sc in ws_schemes(d)) {
      u <- cxr(ws_quiet(df_counting(d, "time", "status", "treat", scheme = sc[[1]], scheme_params = sc[[2]])), sc)
      c1 <- cxr(dcw(d, "one", scheme = sc[[1]], scheme_params = sc[[2]]), sc)
      expect_identical(u, c1)
      f <- u$fit
      # the earlier formula: score_calculation() on the counts at bhat
      y0 <- colSums(outer(f$time[f$z == 0], f$time, ">=")); y1 <- colSums(outer(f$time[f$z == 1], f$time, ">=")) * exp(f$bhat)
      n0 <- diff(c(0, colSums(outer(f$time[f$z == 0 & f$delta == 1], f$time, "<="))))
      n1 <- diff(c(0, colSums(outer(f$time[f$z == 1 & f$delta == 1], f$time, "<="))))
      expect_equal(f$sig_bhat_asy^2, score_calculation(y1, y0, n1, n0, f$wt_rg)$sig2_beta_asy, tolerance = 1e-12)
    }
  }
})
