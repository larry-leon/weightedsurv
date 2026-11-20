# W2: Kaplan-Meier variances, bands and RMST with weight_variance = "estimated".

w2_setup <- function(ties, type = "ipw", n = 400, seed = 11) {
  d <- ws_sim_ps(n, seed, ties = ties)
  pw <- ps_weights(treat ~ x1 + x2, data = d, type = type)
  d$w <- pw$weights
  list(d = d, pw = pw, g = sort(unique(d$time)))
}
kmd_e <- function(d, pw, ...) ws_quiet(KM_diff(d, "time", "status", "treat", weight.name = "w", weight_variance = "estimated", weight_model = pw, ...))
kmd_k <- function(d, ...) ws_quiet(KM_diff(d, "time", "status", "treat", weight.name = "w", ...))
dfc_e <- function(d, pw, ...) ws_quiet(df_counting(d, "time", "status", "treat", weight.name = "w", weight_variance = "estimated", weight_model = pw, ...))

test_that("W2: KM_diff() and df_counting() variances equal the n x m reference (with and without ties)", {
  for (ties in c(FALSE, TRUE)) {
    s <- w2_setup(ties); d <- s$d; pw <- s$pw; g <- s$g
    tf <- survival::aeqSurv(survival::Surv(d$time, d$status)); tm <- tf[, "time"]; dl <- tf[, "status"]
    k1 <- ws_km_ew_ref(tm, dl, d$w, pw$IF, pw$dw, which(d$treat == 1), g)
    k0 <- ws_km_ew_ref(tm, dl, d$w, pw$IF, pw$dw, which(d$treat == 0), g)
    kp <- ws_km_ew_ref(tm, dl, d$w, pw$IF, pw$dw, seq_len(nrow(d)), g)
    vd <- ws_km_ew_diff_ref(k1, k0)
    fe <- kmd_e(d, pw, at_points = g); fk <- kmd_k(d, at_points = g)
    expect_lt(ws_rel(fe$sig2_surv1, k1$sig2), 1e-8); expect_lt(ws_rel(fe$sig2_surv0, k0$sig2), 1e-8)
    expect_lt(ws_rel(fe$sig2_dhat, vd), 1e-8)
    expect_false(isTRUE(all.equal(fe$sig2_dhat, fe$sig2_surv0 + fe$sig2_surv1, tolerance = 1e-6)))   # not the sum of the arm variances
    expect_identical(fe$lower, fe$dhat - stats::qnorm(0.975) * sqrt(fe$sig2_dhat))
    # estimates identical to the "known" call
    for (f in c("surv0", "surv1", "dhat", "at_points")) expect_identical(fe[[f]], fk[[f]])
    expect_identical(fe$weight_variance, "estimated"); expect_identical(fk$weight_variance, "known")
    # df_counting(): arms and pooled at common time points
    de <- dfc_e(d, pw); dk <- ws_quiet(df_counting(d, "time", "status", "treat", weight.name = "w"))
    jd <- match(g, de$at_points)
    expect_lt(ws_rel(de$sig2_surv1[jd], k1$sig2), 1e-8); expect_lt(ws_rel(de$sig2_surv0[jd], k0$sig2), 1e-8)
    expect_lt(ws_rel(de$sig2_survP[jd], kp$sig2), 1e-8)
    for (f in c("surv0", "surv1", "survP", "survG", "ybar0", "ybar1", "nbar0", "nbar1", "riskpoints0", "riskpoints1", "at_points", "time", "w_hat"))
      expect_identical(de[[f]], dk[[f]])
    expect_identical(de$quantile_results$quantile, dk$quantile_results$quantile)
    expect_identical(de$weight_variance, "estimated"); expect_identical(dk$weight_variance, "known")
    expect_null(dk$ew_parts)
    o <- order(tm)
    expect_equal(de$ew_parts$IF, pw$IF[o, , drop = FALSE]); expect_equal(de$ew_parts$dw, pw$dw[o, , drop = FALSE])
  }
})

test_that("W2: df_counting() stratum variances equal the reference for each stratum group", {
  for (ties in c(FALSE, TRUE)) {
    s <- w2_setup(ties); d <- s$d; pw <- s$pw
    tf <- survival::aeqSurv(survival::Surv(d$time, d$status)); tm <- tf[, "time"]; dl <- tf[, "status"]
    de <- dfc_e(d, pw, strata.name = "s")
    o <- order(tm); ts <- tm[o]; ss <- d$s[o]
    for (st in c("hi", "lo")) {
      for (grp in c("0", "1", "P")) {
        rows <- which(d$s == st & (if (grp == "P") TRUE else d$treat == as.numeric(grp)))
        ref <- ws_km_ew_ref(tm, dl, d$w, pw$IF, pw$dw, rows, de$at_points)
        got <- de[[paste0("sig2_surv", grp, "_stratum")]]
        sel <- which(ss == st)
        expect_lt(ws_rel(got[sel], ref$sig2[match(ts[sel], de$at_points)]), 1e-8)
      }
    }
  }
})

test_that("W2: resampled curves follow the reference and their variance matches the analytic one", {
  s <- w2_setup(FALSE); d <- s$d; pw <- s$pw; g <- s$g; n <- nrow(d)
  B <- 4000
  f4 <- kmd_e(d, pw, at_points = g, draws.band = B)
  pos <- f4$sig2_dhat > 0
  r <- stats::median(apply(f4$dhat_star[pos, ], 1, stats::var) / f4$sig2_dhat[pos])
  expect_gt(r, 0.95); expect_lt(r, 1.05)
  # dhat_star = - t(S1 Phi1 - S0 Phi0) %*% G, G holding each subject's row of the multipliers generated today
  tf <- survival::aeqSurv(survival::Surv(d$time, d$status)); tm <- tf[, "time"]; dl <- tf[, "status"]
  k1 <- ws_km_ew_ref(tm, dl, d$w, pw$IF, pw$dw, which(d$treat == 1), g)
  k0 <- ws_km_ew_ref(tm, dl, d$w, pw$IF, pw$dw, which(d$treat == 0), g)
  o <- order(tm); zs <- d$treat[o]
  G <- matrix(0, n, B)
  Gd <- ws_with_seed(8316951, list(G0 = matrix(stats::rnorm(B * sum(zs == 0)), ncol = B), G1 = matrix(stats::rnorm(B * sum(zs == 1)), ncol = B)))
  G[o[zs == 0], ] <- Gd$G0; G[o[zs == 1], ] <- Gd$G1
  ref <- -t(k1$Phi * rep(k1$S, each = n) - k0$Phi * rep(k0$S, each = n)) %*% G
  expect_lt(ws_rel(f4$dhat_star, ref), 1e-8)
  expect_identical(f4$Zdhat_star, f4$dhat_star / sqrt(f4$sig2_dhat))
  expect_identical(f4$sb_upper, f4$dhat + f4$c_alpha_band * sqrt(f4$sig2_dhat))
})

test_that("W2: ipw and stabilised objects from the same model give the same results", {
  for (ties in c(FALSE, TRUE)) {
    s <- w2_setup(ties); d <- s$d; pw <- s$pw; g <- s$g
    pst <- ps_weights(treat ~ x1 + x2, data = d, type = "stabilized"); ds <- d; ds$w <- pst$weights
    fe <- kmd_e(d, pw, at_points = g, draws.band = 300); fs <- kmd_e(ds, pst, at_points = g, draws.band = 300)
    for (f in c("sig2_surv0", "sig2_surv1", "sig2_dhat", "dhat_star")) expect_lt(ws_rel(fs[[f]], fe[[f]]), 1e-8)
    expect_lt(abs(fs$c_alpha_band / fe$c_alpha_band - 1), 1e-8)
  }
})

test_that("W2: with IF set to zero the results are the known-weight ones", {
  for (ties in c(FALSE, TRUE)) {
    s <- w2_setup(ties); d <- s$d; pw <- s$pw; g <- s$g
    pz <- pw; pz$IF <- pw$IF * 0
    f0 <- kmd_e(d, pz, at_points = g, draws.band = 300); fk <- kmd_k(d, at_points = g, draws.band = 300)
    for (f in c("sig2_surv0", "sig2_surv1", "sig2_dhat", "dhat_star", "surv0_star", "surv1_star")) expect_lt(ws_rel(f0[[f]], fk[[f]]), 1e-10)
    d0 <- dfc_e(d, pz, strata.name = "s"); dk <- ws_quiet(df_counting(d, "time", "status", "treat", weight.name = "w", strata.name = "s"))
    for (f in c("sig2_surv0", "sig2_surv1", "sig2_survP", "sig2_surv0_stratum", "sig2_surv1_stratum", "sig2_survP_stratum"))
      expect_lt(ws_rel(d0[[f]], dk[[f]]), 1e-10)
    expect_equal(d0$quantile_results, dk$quantile_results, tolerance = 1e-10)
  }
})

test_that("W2: the influence-based formulas run whatever the weights (no all-equal test)", {
  d <- ws_sim_ps(300, 14)
  pw <- ps_weights(treat ~ 1, data = d); d$w <- pw$weights      # intercept only: the weights are constant within arm
  g <- sort(unique(d$time))
  tf <- survival::aeqSurv(survival::Surv(d$time, d$status)); tm <- tf[, "time"]; dl <- tf[, "status"]
  k1 <- ws_km_ew_ref(tm, dl, d$w, pw$IF, pw$dw, which(d$treat == 1), g)
  k0 <- ws_km_ew_ref(tm, dl, d$w, pw$IF, pw$dw, which(d$treat == 0), g)
  fe <- kmd_e(d, pw, at_points = g)
  expect_lt(ws_rel(fe$sig2_dhat, ws_km_ew_diff_ref(k1, k0)), 1e-8)
})

test_that("W2: cumulative_rmst_bands() and plotKM.band_subgroups() pass the option on; estimates unchanged", {
  s <- w2_setup(FALSE); d <- s$d; pw <- s$pw
  fe <- kmd_e(d, pw, draws.band = 200); fk <- kmd_k(d, draws.band = 200)
  re <- ws_quiet(cumulative_rmst_bands(d, fe, "time", "status", "treat", weight.name = "w", draws_sb = 200, plot = FALSE,
                                       weight_variance = "estimated", weight_model = pw))
  rk <- ws_quiet(cumulative_rmst_bands(d, fk, "time", "status", "treat", weight.name = "w", draws_sb = 200, plot = FALSE))
  expect_identical(re$rmst_time, rk$rmst_time)
  expect_identical(re$weight_variance, "estimated"); expect_identical(rk$weight_variance, "known")
  expect_false(isTRUE(all.equal(re$sig2_rmst_time, rk$sig2_rmst_time)))
  # the second draw set is the estimated one too
  f2 <- kmd_e(d, pw, at_points = fe$at_points, draws.band = 200, seedstart = 99999)
  cs <- apply(f2$dhat_star, 2, ws_cumulative_mean, dt = diff(fe$at_points), rule = "trapezoid") / sqrt(re$sig2_rmst_time)
  expect_equal(unname(re$c_alpha_band), unname(stats::quantile(apply(abs(cs), 2, max, na.rm = TRUE), 0.95)), tolerance = 1e-12)
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  pe <- plotKM.band_subgroups(d, "time", "status", "treat", weight.name = "w", draws.band = 100, sg_labels = "x2 == 1", sg_colors = "blue",
                              weight_variance = "estimated", weight_model = pw)
  pk <- plotKM.band_subgroups(d, "time", "status", "treat", weight.name = "w", draws.band = 100, sg_labels = "x2 == 1", sg_colors = "blue")
  ke <- kmd_e(d, pw, at_points = pe$fit_itt$at_points, draws.band = 100)
  expect_identical(pe$fit_itt$sig2_dhat, ke$sig2_dhat); expect_identical(pe$fit_itt$sb_lower, ke$sb_lower)
  expect_identical(pe$fit_itt$weight_variance, "estimated")
  expect_identical(pe$mean, pk$mean); expect_identical(pe$Dhat_subgroups, pk$Dhat_subgroups)
})

test_that("W2: each misuse is an error that names the argument and the problem", {
  s <- w2_setup(FALSE); d <- s$d; pw <- s$pw
  pst <- ps_weights(treat ~ x1 + x2, data = d, type = "stabilized"); d$wst <- pst$weights
  est <- function(...) KM_diff(d, "time", "status", "treat", weight_variance = "estimated", ...)
  expect_error(est(weight_model = pw), "weight.name must name")
  expect_error(est(weight.name = "w"), "weight_model is missing")
  expect_error(est(weight.name = "w", weight_model = pw$model), "weight_model must be the object returned by ps_weights")
  expect_error(est(weight.name = "wst", weight_model = pw), "weight column 'wst' differs")
  expect_error(KM_diff(d[-1, ], "time", "status", "treat", weight.name = "w", weight_variance = "estimated", weight_model = pw), "399 rows.*fitted on 400")
  d2 <- d[c(2, 1, 3:nrow(d)), ]; d2$treat[1:2] <- c(0, 1) * (d$treat[1] == 1) + c(1, 0) * (d$treat[1] == 0)
  expect_error(KM_diff(d2, "time", "status", "treat", weight.name = "w", weight_variance = "estimated", weight_model = pw), "treatment column 'treat' differs|weight column")
  expect_error(est(weight.name = "w", weight_model = pw, draws = 50), "draws must be 0")
  expect_error(KM_diff(d, "time", "status", "treat", weight.name = "w", weight_model = pw, weight_variance = "bootstrap"), "should be one of")
  expect_error(df_counting(d, "time", "status", "treat", weight.name = "wst", weight_variance = "estimated", weight_model = pw), "df_counting\\(\\).*weight column 'wst' differs")
  expect_error(df_counting(d, "time", "status", "treat", weight.name = "w", draws = 20, weight_variance = "estimated", weight_model = pw), "draws must be 0")
  fk <- kmd_k(d, draws.band = 50); fe <- kmd_e(d, pw, draws.band = 50)
  expect_error(cumulative_rmst_bands(d, fk, "time", "status", "treat", weight.name = "w", draws_sb = 50, plot = FALSE, weight_variance = "estimated", weight_model = pw),
               "fit was made with weight_variance = \"known\"")
  expect_error(cumulative_rmst_bands(d, fe, "time", "status", "treat", weight.name = "w", draws_sb = 50, plot = FALSE), "fit was made with weight_variance = \"estimated\"")
  fold <- fk; fold$weight_variance <- NULL                          # a fit without the element counts as "known"
  expect_silent(ws_quiet(cumulative_rmst_bands(d, fold, "time", "status", "treat", weight.name = "w", draws_sb = 50, plot = FALSE)))
  # with "known", weight_model is ignored
  expect_identical(KM_diff(d, "time", "status", "treat", weight.name = "w", weight_model = "anything")$sig2_dhat, kmd_k(d)$sig2_dhat)
})

test_that("W2: check.seKM draws the comparison plot and raises no standard-error discrepancy warning", {
  s <- w2_setup(FALSE, n = 200); d <- s$d; pw <- s$pw
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  ws <- character(0)
  r <- withCallingHandlers(df_counting(d, "time", "status", "treat", weight.name = "w", check.seKM = TRUE, verbose = TRUE,
                                       weight_variance = "estimated", weight_model = pw),
                           warning = function(w) { ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning") })
  expect_false(any(grepl("se\\(KM\\)", ws)))
})
