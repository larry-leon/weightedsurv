# W4: cox_rhogamma() on a df_counting() result made with weight_variance = "estimated".

w4_setup <- function(ties, type = "ipw", n = 400, seed = 11) {
  d <- ws_sim_ps(n, seed, ties = ties)
  pw <- ps_weights(treat ~ x1 + x2, data = d, type = type); d$w <- pw$weights
  list(d = d, pw = pw,
       de = ws_quiet(df_counting(d, "time", "status", "treat", weight.name = "w", weight_variance = "estimated", weight_model = pw)),
       dk = ws_quiet(df_counting(d, "time", "status", "treat", weight.name = "w")))
}
r1 <- function(a, b) abs(a / b - 1)
P00 <- list(rho = 0, gamma = 0)

test_that("W4: FH(0,0) SE and z equal the coxph dfbeta and score-residual routes (ties, both types)", {
  for (ties in c(FALSE, TRUE)) for (type in c("ipw", "stabilized")) {
    s <- w4_setup(ties, type); d <- s$d; pw <- s$pw
    adj <- function(r) { r <- as.numeric(r); as.numeric(r + pw$IF %*% colSums((r / d$w) * pw$dw)) }
    f0 <- survival::coxph(survival::Surv(time, status) ~ treat, data = d, weights = w, init = 0, iter.max = 0, robust = TRUE, ties = "breslow")
    r <- as.numeric(stats::residuals(f0, type = "score", weighted = TRUE)); chi <- sum(r)^2 / sum(adj(r)^2)
    fb <- survival::coxph(survival::Surv(time, status) ~ treat, data = d, weights = w, robust = TRUE, ties = "breslow")
    se_b <- sqrt(sum(adj(stats::residuals(fb, type = "dfbeta", weighted = TRUE))^2))
    cx <- ws_quiet(cox_rhogamma(s$de, scheme = "fh", scheme_params = P00)); ck <- ws_quiet(cox_rhogamma(s$dk, scheme = "fh", scheme_params = P00))
    expect_lt(r1(cx$fit$sig_bhat_asy, se_b), 1e-6)
    expect_lt(r1(cx$fit$z.score^2, chi), 1e-8)
    expect_identical(cx$z.score, cx$fit$z.score)
    expect_identical(cx$weight_variance, "estimated"); expect_identical(ck$weight_variance, "known")
    # unchanged: estimate, scores, information-based quantities, time weights
    for (f in c("bhat", "u.beta", "u.zero", "wt_rg", "time", "delta", "z", "w_hat")) expect_identical(cx$fit[[f]], ck$fit[[f]])
    expect_identical(cx$hr_ci_asy$hr, ck$hr_ci_asy$hr)
    expect_equal(cx$hr_ci_asy$lower, exp(cx$fit$bhat - stats::qnorm(0.975) * cx$fit$sig_bhat_asy))
    expect_identical(cx$zlogrank_text, paste0("logrank (1-sided) p = ", format_pval(1 - stats::pnorm(cx$z.score), eps = 0.001, digits = 4)))
    expect_identical(cx$cox_text_asy, paste0("HR = ", round(cx$hr_ci_asy$hr, 4), " (", round(cx$hr_ci_asy$lower, 4), ", ", round(cx$hr_ci_asy$upper, 4), ")"))
  }
})

test_that("W4: other schemes equal the reference; z agrees with df_counting()", {
  for (ties in c(FALSE, TRUE)) {
    d <- ws_sim_ps(400, 11, ties = ties); pw <- ps_weights(treat ~ x1 + x2, data = d); d$w <- pw$weights
    for (sc in ws_schemes(d)) {
      de <- ws_quiet(df_counting(d, "time", "status", "treat", weight.name = "w", scheme = sc[[1]], scheme_params = sc[[2]],
                                 weight_variance = "estimated", weight_model = pw))
      cx <- ws_quiet(cox_rhogamma(de, scheme = sc[[1]], scheme_params = sc[[2]]))
      ref <- ws_score_resid_ref(de$time, de$delta, de$z, de$w_hat, K = cx$fit$wt_rg, beta = cx$fit$bhat)
      se <- sqrt(sum(ws_score_ew_ref(ref$psi, de$w_hat, de$ew_parts$IF, de$ew_parts$dw)^2)) / ref$info
      expect_lt(r1(cx$fit$sig_bhat_asy, se), 1e-8)
      expect_lt(r1(cx$fit$z.score, de$z.score), 1e-8)
      r0 <- ws_score_resid_ref(de$time, de$delta, de$z, de$w_hat, K = cx$fit$wt_rg, beta = 0)
      expect_lt(r1(cx$fit$sig2_score, sum(ws_score_ew_ref(r0$psi, de$w_hat, de$ew_parts$IF, de$ew_parts$dw)^2)), 1e-8)
    }
  }
})

test_that("W4: the resampling standard error estimates the analytic one; resampled scores follow the contributions", {
  s <- w4_setup(FALSE)
  cs <- ws_quiet(cox_rhogamma(s$de, scheme = "fh", scheme_params = P00, draws = 4000))
  r <- cs$fit$sig_bhat_star / cs$fit$sig_bhat_asy
  expect_gt(r, 0.95); expect_lt(r, 1.05)
  ck <- ws_quiet(cox_rhogamma(s$dk, scheme = "fh", scheme_params = P00, draws = 0))
  expect_lt(cs$fit$sig_bhat_asy, ck$fit$sig_bhat_asy)
  rs <- stats::sd(cs$fit_resamples$score_star_null) / sqrt(cs$fit$sig2_score)
  expect_gt(rs, 0.95); expect_lt(rs, 1.05)
  de <- s$de; n <- length(de$time); B <- 50
  c50 <- ws_quiet(cox_rhogamma(de, scheme = "fh", scheme_params = P00, draws = B))
  Gd <- ws_with_seed(8316951, list(G0 = matrix(stats::rnorm(B * sum(de$z == 0)), ncol = B), G1 = matrix(stats::rnorm(B * sum(de$z == 1)), ncol = B)))
  G <- matrix(0, n, B); G[de$z == 0, ] <- Gd$G0; G[de$z == 1, ] <- Gd$G1
  psi <- ws_score_ew_ref(ws_score_resid_ref(de$time, de$delta, de$z, de$w_hat, beta = c50$fit$bhat)$psi, de$w_hat, de$ew_parts$IF, de$ew_parts$dw)
  expect_lt(max(abs(c50$fit_resamples$score_star + colSums(G * psi))), 1e-8 * max(abs(c50$fit_resamples$score_star)))
})

test_that("W4: with IF set to zero the results are the known-weight ones; a direct cox_score_rhogamma() call is unchanged", {
  for (ties in c(FALSE, TRUE)) {
    d <- ws_sim_ps(400, 11, ties = ties); pw <- ps_weights(treat ~ x1 + x2, data = d); d$w <- pw$weights
    pz <- pw; pz$IF <- pw$IF * 0
    d0 <- ws_quiet(df_counting(d, "time", "status", "treat", weight.name = "w", weight_variance = "estimated", weight_model = pz))
    dk <- ws_quiet(df_counting(d, "time", "status", "treat", weight.name = "w"))
    c0 <- ws_quiet(cox_rhogamma(d0, scheme = "fh", scheme_params = list(rho = 0, gamma = 1), draws = 100))
    ck <- ws_quiet(cox_rhogamma(dk, scheme = "fh", scheme_params = list(rho = 0, gamma = 1), draws = 100))
    for (f in c("sig_bhat_asy", "z.score", "sig2_score", "sig_bhat_star")) expect_lt(r1(c0$fit[[f]], ck$fit[[f]]), 1e-10)
    expect_lt(max(abs(c0$fit_resamples$score_star - ck$fit_resamples$score_star)), 1e-10 * max(abs(ck$fit_resamples$score_star)))
    de <- ws_quiet(df_counting(d, "time", "status", "treat", weight.name = "w", weight_variance = "estimated", weight_model = pw))
    a <- cox_score_rhogamma(0, de$time, de$delta, de$z, de$w_hat, score_only = FALSE)
    b <- cox_score_rhogamma(0, dk$time, dk$delta, dk$z, dk$w_hat, score_only = FALSE)
    expect_identical(a, b)
  }
})
