# W3: log-rank (and the survival-based outputs of df_counting()) with weight_variance = "estimated".

w3_setup <- function(ties, type = "ipw", n = 400, seed = 11) {
  d <- ws_sim_ps(n, seed, ties = ties)
  pw <- ps_weights(treat ~ x1 + x2, data = d, type = type); d$w <- pw$weights
  list(d = d, pw = pw)
}
w3_e <- function(d, pw, ...) ws_quiet(df_counting(d, "time", "status", "treat", weight.name = "w", weight_variance = "estimated", weight_model = pw, ...))
w3_k <- function(d, ...) ws_quiet(df_counting(d, "time", "status", "treat", weight.name = "w", ...))
r1 <- function(a, b) abs(a / b - 1)
P00 <- list(rho = 0, gamma = 0)

test_that("W3: FH(0,0) z^2 and logrank_results equal the coxph score-residual route (ties, strata, both types)", {
  for (ties in c(FALSE, TRUE)) for (type in c("ipw", "stabilized")) {
    s <- w3_setup(ties, type); d <- s$d; pw <- s$pw
    adj <- function(r) { r <- as.numeric(r); as.numeric(r + pw$IF %*% colSums((r / d$w) * pw$dw)) }
    sres <- function(f) { r <- as.numeric(stats::residuals(f, type = "score", weighted = TRUE)); sum(r)^2 / sum(adj(r)^2) }
    chi  <- sres(survival::coxph(survival::Surv(time, status) ~ treat, data = d, weights = w, init = 0, iter.max = 0, robust = TRUE, ties = "breslow"))
    chis <- sres(survival::coxph(survival::Surv(time, status) ~ treat + survival::strata(s), data = d, weights = w, init = 0, iter.max = 0, robust = TRUE, ties = "breslow"))
    chie <- sres(survival::coxph(survival::Surv(time, status) ~ treat, data = d, weights = w, init = 0, iter.max = 0, robust = TRUE))
    fit  <- survival::coxph(survival::Surv(time, status) ~ treat, data = d, weights = w, robust = TRUE)
    se_e <- sqrt(sum(adj(stats::residuals(fit, type = "dfbeta", weighted = TRUE))^2))
    de <- w3_e(d, pw); des <- w3_e(d, pw, strata.name = "s"); dk <- w3_k(d); dks <- w3_k(d, strata.name = "s")
    expect_lt(r1(de$z.score^2, chi), 1e-8); expect_lt(r1(des$z.score_stratified^2, chis), 1e-8)
    expect_lt(r1(de$logrank_results$chisq, chi), 1e-8); expect_lt(r1(des$logrank_results$chisq, chis), 1e-8)
    expect_identical(de$logrank_results$method, "coxph_robust_score")
    expect_identical(de$logrank_results$pval, 1 - stats::pchisq(de$logrank_results$chisq, 1))
    expect_lt(r1(de$cox_results$score, chie), 1e-8)
    expect_lt(r1(de$cox_results$se, se_e), 1e-8)
    expect_equal(as.numeric(de$cox_results$hr_ci), as.numeric(exp(stats::coef(fit) + c(-1, 1) * stats::qnorm(0.975) * se_e)), tolerance = 1e-10)
    expect_identical(dimnames(de$cox_results$hr_ci), dimnames(dk$cox_results$hr_ci))
    expect_equal(de$cox_results$pval, 2 * stats::pnorm(-abs(stats::coef(fit)[[1]] / se_e)), tolerance = 1e-10)
    expect_identical(names(de$cox_results), names(dk$cox_results))
    # consistency of the fields that follow
    expect_identical(de$sig2_lr, de$get_score$sig2.score); expect_identical(de$z.score, de$get_score$z.score)
    expect_equal(de$z.score, de$lr / sqrt(de$sig2_lr))
    expect_equal(des$z.score_stratified, des$lr_stratified / sqrt(des$sig2_lr_stratified))
    expect_identical(de$zlogrank_text, paste0("logrank (1-sided) p = ", format_pval(1 - stats::pnorm(de$z.score), eps = 0.001, digits = 2)))
    # estimates identical to the "known" call
    expect_identical(de$lr, dk$lr); expect_identical(de$get_score$score, dk$get_score$score)
    expect_identical(des$lr_stratified, dks$lr_stratified)
    expect_identical(de$cox_results$hr, dk$cox_results$hr); expect_identical(de$cox_results$cox_fit$coefficients, dk$cox_results$cox_fit$coefficients)
    expect_false(isTRUE(all.equal(de$sig2_lr, dk$sig2_lr)))
  }
})

test_that("W3: the stratified variance is one sum over all subjects, not a sum of per-stratum variances", {
  s <- w3_setup(FALSE); d <- s$d; pw <- s$pw
  des <- w3_e(d, pw, strata.name = "s")
  psi <- numeric(nrow(d)); o <- order(d$time)
  dd <- d[o, ]; IF <- pw$IF[o, ]; dw <- pw$dw[o, ]
  per <- 0
  for (st in c("hi", "lo")) {
    i <- which(dd$s == st)
    p_s <- ws_score_resid_ref(dd$time[i], dd$status[i], dd$treat[i], dd$w[i])$psi
    psi[i] <- p_s
    per <- per + sum(ws_score_ew_ref(p_s, dd$w[i], IF[i, ], dw[i, ])^2)
  }
  expect_lt(r1(des$sig2_lr_stratified, sum(ws_score_ew_ref(psi, dd$w, IF, dw)^2)), 1e-8)
  expect_gt(r1(des$sig2_lr_stratified, per), 1e-4)
})

test_that("W3: other schemes equal the reference with the scheme's time weights (with and without ties)", {
  for (ties in c(FALSE, TRUE)) {
    s <- w3_setup(ties); d <- s$d; pw <- s$pw
    for (sc in ws_schemes(d)[c("fh01", "MB", "schemper")]) {
      de <- w3_e(d, pw, scheme = sc[[1]], scheme_params = sc[[2]]); dk <- w3_k(d, scheme = sc[[1]], scheme_params = sc[[2]])
      K <- ws_quiet(cox_rhogamma(dk, scheme = sc[[1]], scheme_params = sc[[2]], draws = 0))$fit$wt_rg
      psi <- ws_score_resid_ref(de$time, de$delta, de$z, de$w_hat, K = K, beta = 0)$psi
      expect_lt(r1(de$sig2_lr, sum(ws_score_ew_ref(psi, de$w_hat, de$ew_parts$IF, de$ew_parts$dw)^2)), 1e-8)
      expect_identical(de$lr, dk$lr)
    }
  }
})

test_that("W3: wlr_cumulative() equals the direct subject-by-time sum; last value and censored-data values", {
  for (ties in c(FALSE, TRUE)) {
    s <- w3_setup(ties, n = 200, seed = 31); d <- s$d; pw <- s$pw
    for (P in list(P00, list(rho = 0, gamma = 1))) {
      de <- w3_e(d, pw, scheme = "fh", scheme_params = P)
      one <- wlr_cumulative(de, "fh", P); cum <- wlr_cumulative(de, "fh", P, return_cumulative = TRUE)
      expect_lt(r1(one$sig2.score, de$sig2_lr), 1e-10)
      expect_lt(r1(cum$sig2.score[length(cum$sig2.score)], de$sig2_lr), 1e-10)
      expect_identical(cum$score, wlr_cumulative(w3_k(d, scheme = "fh", scheme_params = P), "fh", P, return_cumulative = TRUE)$score)
      # direct: c_i(t) for every subject and grid time, then sum_i (c_i(t) + IF_i' B(t))^2
      U <- de$time; D <- de$delta; Z <- de$z; a <- de$w_hat; g <- de$at_points; IF <- de$ew_parts$IF; dw <- de$ew_parts$dw
      Kobs <- ws_quiet(cox_rhogamma(de, scheme = "fh", scheme_params = P, draws = 0))$fit$wt_rg
      sd <- sort(unique(U)); Ks <- Kobs[match(sd, U)]
      S0 <- vapply(sd, function(x) sum(a[U >= x]), 0); zb <- vapply(sd, function(x) sum(a[U >= x & Z == 1]), 0) / S0
      dL <- vapply(sd, function(x) sum(a[U == x & D == 1]), 0) / S0
      V <- vapply(g, function(t) {
        ci <- vapply(seq_along(U), function(i) {
          ev <- if (D[i] == 1 && U[i] <= t) Ks[match(U[i], sd)] * (Z[i] - zb[match(U[i], sd)]) else 0
          k <- sd <= min(t, U[i])
          a[i] * (ev - sum(Ks[k] * (Z[i] - zb[k]) * dL[k]))
        }, 0)
        sum((ci + as.numeric(IF %*% colSums((ci / a) * dw)))^2)
      }, 0)
      pos <- V > 1e-12
      expect_lt(max(abs(cum$sig2.score[pos] / V[pos] - 1)), 1e-10)
      expect_lt(max(abs(cum$sig2.score[!pos] - V[!pos])), 1e-12)
    }
    # FH(0,0): the value at t is sig2_lr of the data censored at t (same weights and model)
    de <- w3_e(d, pw); cum <- wlr_cumulative(de, "fh", P00, return_cumulative = TRUE)
    for (q in c(0.3, 0.5, 0.8)) {
      t0 <- unname(stats::quantile(d$time, q)); j0 <- max(which(de$at_points <= t0))
      dcen <- d; dcen$status <- d$status * (d$time <= t0); dcen$time <- pmin(d$time, t0)
      dn <- w3_e(dcen, pw)
      expect_lt(r1(cum$sig2.score[j0], dn$sig2_lr), 1e-8)
    }
  }
})

test_that("W3: wlr_dhat_estimates() agrees with KM_diff() and df_counting(); correlation within [-1, 1]", {
  for (ties in c(FALSE, TRUE)) {
    s <- w3_setup(ties); d <- s$d; pw <- s$pw
    de <- w3_e(d, pw); dk <- w3_k(d)
    km <- ws_quiet(KM_diff(d, "time", "status", "treat", weight.name = "w", at_points = de$at_points, weight_variance = "estimated", weight_model = pw))
    for (q in c(0.2, 0.5, 0.9)) {
      tz <- unname(stats::quantile(d$time, q)); j <- max(which(de$at_points <= tz))
      we <- wlr_dhat_estimates(de, tzero = tz); wk <- wlr_dhat_estimates(dk, tzero = tz)
      expect_lt(r1(we$sig2_dhat, km$sig2_dhat[j]), 1e-8)
      expect_lt(r1(we$sig2.score, de$sig2_lr), 1e-10); expect_equal(we$z.score, de$z.score, tolerance = 1e-10)
      expect_lte(abs(we$cor_wlr_dhat), 1)
      expect_equal(we$cor_wlr_dhat, we$cov_wlr_dhat / sqrt(we$sig2.score * we$sig2_dhat))
      expect_identical(we$score, wk$score); expect_identical(we$dhat, wk$dhat)
      # the covariance against the direct sum over subjects
      U <- de$time; a <- de$w_hat; IF <- de$ew_parts$IF; dw <- de$ew_parts$dw; g <- de$at_points
      psi <- ws_score_ew_ref(ws_score_resid_ref(U, de$delta, de$z, a)$psi, a, IF, dw)
      k1 <- ws_km_ew_ref(U, de$delta, a, IF, dw, which(de$z == 1), g); k0 <- ws_km_ew_ref(U, de$delta, a, IF, dw, which(de$z == 0), g)
      dd <- -(k1$S[j] * k1$Phi[, j] - k0$S[j] * k0$Phi[, j])
      expect_lt(r1(we$cov_wlr_dhat, sum(-psi * dd)), 1e-8)
    }
  }
})

test_that("W3: with IF set to zero the results are the known-weight ones; a bare table of counts is unchanged", {
  for (ties in c(FALSE, TRUE)) {
    s <- w3_setup(ties); d <- s$d; pw <- s$pw; pz <- pw; pz$IF <- pw$IF * 0
    d0 <- w3_e(d, pz, strata.name = "s"); dk <- w3_k(d, strata.name = "s")
    for (f in c("sig2_lr", "z.score", "sig2_lr_stratified", "z.score_stratified")) expect_lt(r1(d0[[f]], dk[[f]]), 1e-10)
    expect_lt(r1(d0$logrank_results$chisq, dk$logrank_results$chisq), 1e-8)
    expect_lt(r1(d0$cox_results$se, dk$cox_results$se), 1e-8)
    expect_lt(r1(d0$cox_results$score, dk$cox_results$score), 1e-8)
    expect_equal(as.numeric(d0$cox_results$hr_ci), as.numeric(dk$cox_results$hr_ci), tolerance = 1e-8)
    c0 <- wlr_cumulative(d0, "fh", P00, return_cumulative = TRUE); ck <- wlr_cumulative(dk, "fh", P00, return_cumulative = TRUE)
    expect_lt(max(abs(c0$sig2.score - ck$sig2.score)), 1e-10 * max(ck$sig2.score))
    tz <- stats::median(d$time)
    w0 <- wlr_dhat_estimates(d0, tzero = tz); wk <- wlr_dhat_estimates(dk, tzero = tz)
    for (f in c("sig2.score", "sig2_dhat", "cov_wlr_dhat", "cor_wlr_dhat")) expect_lt(r1(w0[[f]], wk[[f]]), 1e-10)
    de <- w3_e(d, pw)
    bare <- data.frame(at_points = de$at_points, nbar0 = de$nbar0, nbar1 = de$nbar1, ybar0 = de$ybar0, ybar1 = de$ybar1, survP = de$survP)
    bk <- data.frame(at_points = dk$at_points, nbar0 = dk$nbar0, nbar1 = dk$nbar1, ybar0 = dk$ybar0, ybar1 = dk$ybar1, survP = dk$survP)
    expect_identical(wlr_cumulative(bare, "fh", P00), wlr_cumulative(bk, "fh", P00))
  }
})
