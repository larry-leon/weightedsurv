#!/usr/bin/env Rscript
# oc_sim.R -- operating characteristics of weighted inference (TASK_km_oc_evaluation.md), Parts A-D
# and the section 9 diagnostics. One R process per code version.
#
#   Rscript dev/audit/oc/oc_sim.R --version=cran|proposed --pkg=<path> [--cells=1:19] [--reps=2000]
#                                 [--cores=12] [--out=raw] [--check] [--gates]
#
#   --version=cran      the package at --pkg is weightedsurv 0.1.0 (worktree of 166919e)
#   --version=proposed  the package at --pkg is the branch at HEAD; also runs every script-side estimator
#   --check             re-runs replicates 1-10 of cells 1 and 5 and compares with the stored per-replicate
#                       results (within 1e-10); exit status 0 only if they agree
#   --gates             (proposed only) the single-dataset checks of Gate G1: items 1, 2, 5a (first part), 5b, 5c, 5d
#
# Per-replicate results: <out>/<version>_cellNN.rds, a list with `res` (replicates x fields numeric matrix),
# `chk` (the per-replicate data checksum) and timing. Field names are <estimator>.<quantity>.<field>[.<month>].

suppressPackageStartupMessages(library(survival))
.args <- commandArgs(FALSE)
.here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", .args, value = TRUE)[1])))
source(file.path(.here, "oc_common.R"))
opt <- parse_args(list(version = NA, pkg = NA, cells = "1:19", reps = "2000", cores = "12", out = "raw",
                       check = FALSE, gates = FALSE))
stopifnot(opt$version %in% c("cran", "proposed"), !is.na(opt$pkg))
VERSION <- opt$version
suppressMessages(pkgload::load_all(opt$pkg, quiet = TRUE, export_all = TRUE, helpers = FALSE, attach_testthat = FALSE))
PKG_VERSION <- as.character(utils::packageVersion("weightedsurv"))
stopifnot((VERSION == "cran") == (PKG_VERSION == "0.1.0"))
grDevices::pdf(NULL)                                   # no plot file, whatever a call draws
Q <- qnorm(0.975)
NT <- length(tt)

named <- function(prefix, x) stats::setNames(as.numeric(x), paste0(prefix, ".", if (length(x) == NT) tt else seq_along(x)))
na_t <- function(prefix) named(prefix, rep(NA_real_, NT))

# Band summary over the grid: a grid point counts only where the band has positive width (both limits
# finite, upper above lower). cov_in: truth inside at every such point in [3, 42]; cov_any: at every such
# point; width: mean width over such points in [3, 42]. nonfin: grid points with a non-finite limit.
band_stats <- function(prefix, grid, lo, hi, truth, calpha) {
  pos <- is.finite(lo) & is.finite(hi) & hi > lo
  inr <- grid >= band_range[1] & grid <= band_range[2]
  cov <- truth >= lo & truth <= hi
  v <- c(cov_in = as.numeric(all(cov[pos & inr])), cov_any = as.numeric(all(cov[pos])),
         width = if (any(pos & inr)) mean((hi - lo)[pos & inr]) else NA_real_,
         npos_in = sum(pos & inr), npos = sum(pos), nonfin = sum(!is.finite(lo) | !is.finite(hi)),
         calpha = as.numeric(calpha), fail = 0)
  stats::setNames(v, paste0(prefix, ".", names(v)))
}
band_fail <- function(prefix) {
  v <- c(cov_in = 0, cov_any = 0, width = NA, npos_in = NA, npos = NA, nonfin = NA, calpha = NA, fail = 1)
  stats::setNames(v, paste0(prefix, ".", names(v)))
}

# ---- exact calls (same arguments under both versions, except arguments 0.1.0 does not have) ----
call_km <- function(df, wn, r, draws = 0, draws.band = 500)
  KM_diff(df, tte.name = "time", event.name = "status", treat.name = "treat", weight.name = wn, alpha = 0.05,
          draws = draws, draws.band = draws.band, seedstart = 8316951 + r, show_resamples = FALSE)
call_rmst <- function(df, fit, wn, ...)
  cumulative_rmst_bands(df, fit, tte.name = "time", event.name = "status", treat.name = "treat", weight.name = wn,
                        draws_sb = 500, plot = FALSE, ...)
call_dc <- function(df, wn, g)
  df_counting(df, tte.name = "time", event.name = "status", treat.name = "treat", weight.name = wn,
              scheme = "fh", scheme_params = list(rho = 0, gamma = g))
call_cox <- function(dc, g, draws = 0)
  cox_rhogamma(dc, scheme = "fh", scheme_params = list(rho = 0, gamma = g), draws = draws)

# ---- Part A / B blocks ----
km_fields <- function(est, fit, jj, k, band = TRUE) {
  if (inherits(fit, "error")) {
    return(c(stats::setNames(1, paste0(est, ".km.fail")),
             unlist(lapply(c("S0.est", "S0.se", "S1.est", "S1.se", "d.est", "d.se", "d.lo", "d.hi"),
                           function(f) na_t(paste0(est, ".", f)))),
             if (band) band_fail(paste0(est, ".dband"))))
  }
  c(stats::setNames(0, paste0(est, ".km.fail")),
    named(paste0(est, ".S0.est"), fit$surv0[jj]), named(paste0(est, ".S0.se"), sqrt(fit$sig2_surv0[jj])),
    named(paste0(est, ".S1.est"), fit$surv1[jj]), named(paste0(est, ".S1.se"), sqrt(fit$sig2_surv1[jj])),
    named(paste0(est, ".d.est"), fit$dhat[jj]), named(paste0(est, ".d.se"), sqrt(fit$sig2_dhat[jj])),
    named(paste0(est, ".d.lo"), fit$lower[jj]), named(paste0(est, ".d.hi"), fit$upper[jj]),
    if (band) band_stats(paste0(est, ".dband"), fit$at_points, fit$sb_lower, fit$sb_upper,
                         delta_true(fit$at_points, k), fit$c_alpha_band))
}
rmst_fields <- function(est, rm, jj, Rgrid) {
  if (inherits(rm, "error")) {
    return(c(stats::setNames(1, paste0(est, ".rmst.fail")),
             unlist(lapply(c("R.est", "R.se", "R.lo", "R.hi"), function(f) na_t(paste0(est, ".", f)))),
             band_fail(paste0(est, ".rband"))))
  }
  c(stats::setNames(0, paste0(est, ".rmst.fail")),
    named(paste0(est, ".R.est"), rm$rmst_time[jj]), named(paste0(est, ".R.se"), sqrt(rm$sig2_rmst_time[jj])),
    named(paste0(est, ".R.lo"), rm$rmst_time_lower[jj]), named(paste0(est, ".R.hi"), rm$rmst_time_upper[jj]),
    band_stats(paste0(est, ".rband"), rm$at_points, rm$rmst_time_sb_lower, rm$rmst_time_sb_upper, Rgrid, rm$c_alpha_band))
}
cox_fields <- function(est, g, b = NA, se = NA, lo = NA, hi = NA, z = NA, fail = 0) {
  v <- c(b = b, se = se, lo = lo, hi = hi, z = z, fail = fail)
  stats::setNames(as.numeric(v), paste0(est, ".cox", g, ".", names(v)))
}

# ---- script-side pieces for Part D and the diagnostics (HEAD only: use ws_km_psi_parts()) ----
# n_a x m matrix of psi_i(t_j) for one arm (rows in the order of U) at at_points, computed on the arm's
# grid G = sort(unique(c(U, at_points))) exactly as KM_diff() does (km_arm_on_grid()).
psi_mat <- function(U, D, W, at_points) {
  G <- sort(unique(c(U, at_points)))
  ybar <- colSums(outer(U, G, FUN = ">=") * W)
  nbar <- colSums(outer(U[D == 1], G, FUN = "<=") * W[D == 1])
  parts <- ws_km_psi_parts(U, D, G, ybar, nbar)
  m <- length(G); n <- length(U)
  P <- matrix(rep(-parts$H, each = n), n, m)           # -H[j] while j < k_i
  ex <- col(P) >= parts$k                              # from j = k_i on: the exit value
  P[ex] <- parts$exit[row(P)[ex]]
  P[parts$k < 1, ] <- 0
  P[, match(at_points, G), drop = FALSE]
}
# the multipliers KM_diff() uses (seed, draws), as an n x B matrix with rows in data order
km_multipliers <- function(time_fix, treat, seed, B) {
  ord <- if (is.unsorted(time_fix)) order(time_fix) else seq_along(time_fix)
  zs <- treat[ord]
  n0 <- sum(zs == 0); n1 <- sum(zs == 1)
  set.seed(seed)
  G0 <- matrix(rnorm(B * n0), ncol = B)
  G1 <- matrix(rnorm(B * n1), ncol = B)
  G <- matrix(0, length(treat), B)
  G[ord[zs == 0], ] <- G0
  G[ord[zs == 1], ] <- G1
  G
}
# RMST pointwise limits and band from the estimate and two sets of resampled difference curves: the
# steps of cumulative_rmst_bands() at HEAD (copied).
rmst_from_stars <- function(at_points, dhat, star1, star2, mean_rule = "trapezoid") {
  dt <- diff(at_points)
  rmst_time <- ws_cumulative_mean(dhat, dt, mean_rule)
  cum_rmst_stars <- apply(star1, 2, ws_cumulative_mean, dt = dt, rule = mean_rule)
  sig2_rmst_time <- apply(cum_rmst_stars, 1, var, na.rm = TRUE)
  cum_rmst_stars2 <- apply(star2, 2, ws_cumulative_mean, dt = dt, rule = mean_rule)
  rmst_time_star <- cum_rmst_stars2 / sqrt(sig2_rmst_time)
  sups <- apply(abs(rmst_time_star), 2, max, na.rm = TRUE)
  c_alpha_band <- quantile(sups, c(0.95))
  list(at_points = at_points, rmst_time = rmst_time, sig2_rmst_time = sig2_rmst_time,
       rmst_time_lower = rmst_time - 1.96 * sqrt(sig2_rmst_time), rmst_time_upper = rmst_time + 1.96 * sqrt(sig2_rmst_time),
       c_alpha_band = c_alpha_band,
       rmst_time_sb_lower = rmst_time - c_alpha_band * sqrt(sig2_rmst_time),
       rmst_time_sb_upper = rmst_time + c_alpha_band * sqrt(sig2_rmst_time))
}
# band for delta from resampled difference curves and a variance: the steps of KM_diff() at HEAD (copied)
dband_from_stars <- function(dhat, dhat_star, sig2_dhat, alpha = 0.05) {
  Zdhat_star <- dhat_star / sqrt(sig2_dhat)
  sups <- apply(abs(Zdhat_star), 2, max, na.rm = TRUE)
  c_alpha_band <- quantile(sups, 1 - alpha)
  list(c_alpha_band = c_alpha_band, sb_lower = dhat - c_alpha_band * sqrt(sig2_dhat),
       sb_upper = dhat + c_alpha_band * sqrt(sig2_dhat))
}
# Everything of section 7 for the Kaplan-Meier quantities, for a given IF (n x p). `fit` and `fit2` are the
# KM_diff() fits with seedstart = 8316951 + r and 99999 (the two draw sets of cumulative_rmst_bands()).
ew_km_all <- function(df, fit, fit2, r, IF, dw, B = 500) {
  tf <- aeqSurv(Surv(df$time, df$status))
  time_fix <- tf[, "time"]; delta_fix <- tf[, "status"]
  idx0 <- which(df$treat == 0); idx1 <- which(df$treat == 1)
  P0 <- psi_mat(time_fix[idx0], delta_fix[idx0], df$w[idx0], fit$at_points)
  P1 <- psi_mat(time_fix[idx1], delta_fix[idx1], df$w[idx1], fit$at_points)
  ek <- ew_km(P0, P1, fit$surv0, fit$surv1, idx0, idx1, df$w, IF, dw)
  G <- km_multipliers(time_fix, df$treat, 8316951 + r, B)
  G2 <- km_multipliers(time_fix, df$treat, 99999, B)
  # the script's own resampled curves (Gate 5b), and the adjusted ones
  own0 <- -fit$surv0 * crossprod(P0, df$w[idx0] * G[idx0, , drop = FALSE])
  own1 <- -fit$surv1 * crossprod(P1, df$w[idx1] * G[idx1, , drop = FALSE])
  adj <- function(star, S, Dg, Gm) star - (S * Dg) %*% (t(IF) %*% Gm)
  s0 <- adj(fit$surv0_star, fit$surv0, ek$Dg0, G);   s1 <- adj(fit$surv1_star, fit$surv1, ek$Dg1, G)
  t0 <- adj(fit2$surv0_star, fit$surv0, ek$Dg0, G2); t1 <- adj(fit2$surv1_star, fit$surv1, ek$Dg1, G2)
  dstar <- s1 - s0; dstar2 <- t1 - t0
  list(ek = ek, P0 = P0, P1 = P1, idx0 = idx0, idx1 = idx1, own0 = own0, own1 = own1,
       band = dband_from_stars(fit$dhat, dstar, ek$sig2_dhat),
       rmst = rmst_from_stars(fit$at_points, fit$dhat, dstar, dstar2, "trapezoid"))
}
# section 9: largest single share of the influence sum at month t, one arm
infl_share <- function(U, D, W, t) {
  G <- sort(unique(U))
  ybar <- colSums(outer(U, G, FUN = ">=") * W)
  nbar <- colSums(outer(U[D == 1], G, FUN = "<=") * W[D == 1])
  parts <- ws_km_psi_parts(U, D, G, ybar, nbar)
  j <- findInterval(t, G)
  if (j < 1) return(NA_real_)
  psi <- ifelse(parts$k < 1, 0, ifelse(parts$k <= j, parts$exit, -parts$H[j]))
  v <- (W * psi)^2
  if (sum(v) > 0) max(v) / sum(v) else NA_real_
}

# ---- one replicate ----
one_rep <- function(k, r) {
  cl <- cells[k, ]
  d <- make_data(k, r); df <- d$df; wn <- d$wname
  main <- if (VERSION == "cran") "cran" else "proposed"
  tmax <- c(max(df$time[df$treat == 0], -Inf), max(df$time[df$treat == 1], -Inf))
  out <- named("evaluable", as.numeric(min(tmax) >= tt))

  # Part A and the band for delta
  fit <- tryCatch(call_km(df, wn, r), error = function(e) e)
  jj <- if (inherits(fit, "error")) rep(NA_integer_, NT) else last_le(fit$at_points, tt)
  out <- c(out, km_fields(main, fit, jj, k))
  if (VERSION == "proposed" && cl$d20) {
    fit20 <- tryCatch(call_km(df, wn, r, draws = 20), error = function(e) e)
    out <- c(out, km_fields("proposed_d20", fit20, jj, k))
  }

  # RMST (not cell 12)
  Rgrid <- NULL
  if (cl$do_rmst) {
    if (inherits(fit, "error")) {
      out <- c(out, rmst_fields(main, fit, jj, NULL), na_t("R.truth"))
      if (VERSION == "proposed") out <- c(out, rmst_fields("proposed_rect", fit, jj, NULL))
    } else {
      Rgrid <- R_true(fit$at_points, k)
      out <- c(out, named("R.truth", Rgrid[jj]))
      rm <- tryCatch(call_rmst(df, fit, wn), error = function(e) e)
      out <- c(out, rmst_fields(main, rm, jj, Rgrid))
      if (VERSION == "proposed") {
        rmr <- tryCatch(call_rmst(df, fit, wn, mean_rule = "rectangle"), error = function(e) e)
        out <- c(out, rmst_fields("proposed_rect", rmr, jj, Rgrid))
      }
    }
  }

  # Part D, Kaplan-Meier quantities (proposed_ew)
  ep <- NULL
  if (VERSION == "proposed" && cl$ew) {
    ep <- ew_parts(d$g, df$treat)
    ew <- tryCatch({
      fit2 <- KM_diff(df, tte.name = "time", event.name = "status", treat.name = "treat", weight.name = wn,
                      at_points = fit$at_points, alpha = 0.05, risk.points = NULL, modify_tau = FALSE,
                      draws.band = 500, seedstart = 99999, show_resamples = FALSE)
      ew_km_all(df, fit, fit2, r, ep$IF, ep$dw)
    }, error = function(e) e)
    if (inherits(ew, "error")) {
      out <- c(out, km_fields("proposed_ew", ew, jj, k), rmst_fields("proposed_ew", ew, jj, NULL))
    } else {
      se0 <- sqrt(ew$ek$sig2_surv0); se1 <- sqrt(ew$ek$sig2_surv1); sed <- sqrt(ew$ek$sig2_dhat)
      fe <- list(at_points = fit$at_points, surv0 = fit$surv0, surv1 = fit$surv1, dhat = fit$dhat,
                 sig2_surv0 = ew$ek$sig2_surv0, sig2_surv1 = ew$ek$sig2_surv1, sig2_dhat = ew$ek$sig2_dhat,
                 lower = fit$dhat - Q * sed, upper = fit$dhat + Q * sed,
                 sb_lower = ew$band$sb_lower, sb_upper = ew$band$sb_upper, c_alpha_band = ew$band$c_alpha_band)
      out <- c(out, km_fields("proposed_ew", fe, jj, k), rmst_fields("proposed_ew", ew$rmst, jj, Rgrid))
    }
  }

  # Part C (not cell 12)
  if (cl$do_cox) {
    for (g in c(0, 1)) {
      res <- tryCatch({
        dc <- suppressWarnings(call_dc(df, wn, g))
        cx <- suppressWarnings(call_cox(dc, g))
        if (is.null(cx$fit) || !is.finite(cx$fit$bhat)) stop("root finding failed")
        list(dc = dc, cx = cx)
      }, error = function(e) e)
      pk <- if (VERSION == "cran") "cran" else "pkg"   # the package's own Cox output
      if (inherits(res, "error")) {
        out <- c(out, cox_fields(pk, g, fail = 1))
        if (VERSION == "proposed") {
          out <- c(out, cox_fields("proposed", g, fail = 1))
          if (cl$ew) out <- c(out, cox_fields("proposed_ew", g, fail = 1))
        }
      } else {
        cx <- res$cx; dc <- res$dc
        out <- c(out, cox_fields(pk, g, cx$fit$bhat, cx$fit$sig_bhat_asy, log(cx$hr_ci_asy$lower), log(cx$hr_ci_asy$upper), dc$z.score))
        if (VERSION == "proposed") {
          cr <- cox_robust(cx$fit); cr0 <- cox_robust(cx$fit, 0)
          out <- c(out, cox_fields("proposed", g, cx$fit$bhat, cr$se, cx$fit$bhat - Q * cr$se, cx$fit$bhat + Q * cr$se, cr0$z),
                   stats::setNames(abs(cr$U) / sqrt(sum(cr$psi^2)), paste0("proposed.cox", g, ".rootcheck")))
          if (cl$ew) {
            ord <- order(df$time)
            aligned <- isTRUE(all.equal(as.numeric(cx$fit$time), df$time[ord])) && all(cx$fit$z == df$treat[ord])
            if (aligned) {
              ewpsi <- function(c0) c0$psi + as.numeric(ep$IF[ord, , drop = FALSE] %*% colSums((c0$psi / cx$fit$w_hat) * ep$dw[ord, , drop = FALSE]))
              se_ew <- sqrt(sum(ewpsi(cr)^2)) / cr$info
              z_ew <- cr0$U / sqrt(sum(ewpsi(cr0)^2))
              out <- c(out, cox_fields("proposed_ew", g, cx$fit$bhat, se_ew, cx$fit$bhat - Q * se_ew, cx$fit$bhat + Q * se_ew, z_ew))
            } else out <- c(out, cox_fields("proposed_ew", g, fail = 1))
          }
        }
      }
      if (VERSION == "cran" && g == 0) {
        # cran_mean1: weights rescaled to mean 1 (U cells: the weights are already 1, same calls)
        df1 <- df; df1$w1 <- df$w * nrow(df) / sum(df$w)
        r1 <- tryCatch({
          dc1 <- suppressWarnings(call_dc(df1, if (is.null(wn)) NULL else "w1", g))
          cx1 <- suppressWarnings(call_cox(dc1, g))
          if (is.null(cx1$fit) || !is.finite(cx1$fit$bhat)) stop("root finding failed")
          cox_fields("cran_mean1", g, cx1$fit$bhat, cx1$fit$sig_bhat_asy, log(cx1$hr_ci_asy$lower), log(cx1$hr_ci_asy$upper), dc1$z.score)
        }, error = function(e) cox_fields("cran_mean1", g, fail = 1))
        out <- c(out, r1)
        # cran_star: the package's resampling SE, first 500 replicates
        if (r <= 500) {
          rs <- if (inherits(res, "error")) cox_fields("cran_star", g, fail = 1) else tryCatch({
            set.seed(3e7 + 1e6 * k + r)
            cxs <- suppressWarnings(call_cox(res$dc, g, draws = 200))
            cox_fields("cran_star", g, log(cxs$hr_ci_star$hr), cxs$fit$sig_bhat_star, log(cxs$hr_ci_star$lower), log(cxs$hr_ci_star$upper), NA)
          }, error = function(e) cox_fields("cran_star", g, fail = 1))
          out <- c(out, rs)
        } else out <- c(out, cox_fields("cran_star", g, fail = NA))
      }
    }
  }

  # section 9 diagnostics (HEAD process)
  if (VERSION == "proposed") {
    tf <- aeqSurv(Surv(df$time, df$status))
    for (a in 0:1) {
      i <- df$treat == a; w <- df$w[i]
      out <- c(out, stats::setNames(c(sum(w)^2 / (sum(i) * sum(w^2)), infl_share(tf[i, "time"], tf[i, "status"], w, 18)),
                                    paste0("diag.", c("ess", "share18"), a)))
    }
  }
  attr(out, "chk") <- d$chk
  out
}

run_reps <- function(k, reps, cores) {
  f <- function(r) tryCatch(one_rep(k, r), error = function(e) structure(list(msg = conditionMessage(e), r = r), class = "rep_error"))
  L <- if (cores > 1) parallel::mclapply(reps, f, mc.cores = cores, mc.preschedule = TRUE) else lapply(reps, f)
  bad <- vapply(L, function(x) !is.numeric(x), logical(1))
  if (any(bad)) {                                      # a worker died or a replicate failed outside the caught calls
    for (i in which(bad)) L[[i]] <- f(reps[i])
    bad <- vapply(L, function(x) !is.numeric(x), logical(1))
    if (any(bad)) stop(sprintf("cell %d: replicate %d failed outside the caught calls: %s", k, reps[which(bad)[1]], L[[which(bad)[1]]]$msg))
  }
  nm <- names(L[[1]])
  stopifnot(all(vapply(L, function(x) identical(names(x), nm), logical(1))) || k > 0)
  allnm <- unique(unlist(lapply(L, names)))            # cran_star fields exist for every replicate; names are fixed per cell
  res <- t(vapply(L, function(x) as.numeric(x[allnm]), numeric(length(allnm))))
  colnames(res) <- allnm; rownames(res) <- reps
  list(res = res, chk = vapply(L, function(x) attr(x, "chk"), character(1)))
}

# ---- Gate G1 single-dataset checks ----
if (isTRUE(opt$gates)) {
  stopifnot(VERSION == "proposed")
  fails <- character(0)
  say <- function(...) cat(sprintf(...), "\n", sep = "")
  # 1. Truth, on one sample of 10^6 per scenario with both potential outcomes
  for (sc in c("ph", "null", "delay")) {
    set.seed(match(sc, c("ph", "null", "delay")) + 20261119)
    N <- 1e6; X <- rnorm(N); H <- -log(pnorm(0.6 * X + 0.8 * rnorm(N)))
    T0 <- LamInv(H, 0, sc); T1 <- LamInv(H, 1, sc)
    kk <- match(sc, c("ph", "null", "delay")); kcell <- c(1, 13, 17)[kk]
    d1 <- max(sapply(tt, function(t) max(abs(mean(T0 > t) - S_true(t, 0, kcell)), abs(mean(T1 > t) - S_true(t, 1, kcell)))))
    say("G1.1 truth, %s: max |mean(T^a > t) - S_a(t)| = %.5f (limit 0.002)", sc, d1)
    if (d1 > 0.002) fails <- c(fails, "G1.1 marginal")
    if (sc == "ph") {
      ws <- exp(0.5 * X - 0.125)
      d3 <- max(sapply(tt, function(t) max(abs(weighted.mean(T0 > t, ws) - S_true(t, 0, 3)), abs(weighted.mean(T1 > t, ws) - S_true(t, 1, 3)))))
      say("G1.1 truth, cell 3: max |weighted mean - formula| = %.5f (limit 0.003)", d3)
      if (d3 > 0.003) fails <- c(fails, "G1.1 cell 3")
    }
    e <- plogis(-0.5 + X); A <- rbinom(N, 1, e); w <- ifelse(A == 1, 1 / e, 1 / (1 - e))
    dc <- max(sapply(tt, function(t) max(abs(weighted.mean((T0 > t)[A == 0], w[A == 0]) - S_true(t, 0, kcell)),
                                         abs(weighted.mean((T1 > t)[A == 1], w[A == 1]) - S_true(t, 1, kcell)))))
    say("G1.1 truth, %s, conf with true weights: max |weighted mean among treat == a - S_a(t)| = %.5f (limit 0.004)", sc, dc)
    if (dc > 0.004) fails <- c(fails, "G1.1 conf")
  }
  # 2. cox_robust(), one cell-5 dataset
  d <- make_data(5, 1); df <- d$df
  for (g in c(0, 1)) {
    cx <- call_cox(call_dc(df, "w", g), g); cr <- cox_robust(cx$fit)
    rc <- abs(cr$U) / sqrt(sum(cr$psi^2))
    say("G1.2 cox_robust, FH(0,%d): |U| / sqrt(sum(psi^2)) at bhat = %.3g (limit 1e-6)", g, rc)
    if (rc > 1e-6) fails <- c(fails, "G1.2 root")
    if (g == 0) {
      se_ref <- sqrt(coxph(Surv(time, status) ~ treat, data = df, weights = w, robust = TRUE, ties = "breslow")$var)
      rel <- abs(cr$se - se_ref) / se_ref
      say("G1.2 cox_robust, FH(0,0): SE %.10f, coxph(robust = TRUE) %.10f, relative difference %.3g (limit 1e-6)", cr$se, se_ref, rel)
      if (rel > 1e-6) fails <- c(fails, "G1.2 se")
    }
  }
  # 5a (first part), 5b, 5c, 5d: one cell-8 dataset
  r <- 1; d <- make_data(8, r); df <- d$df
  fit <- call_km(df, "w", r)
  fit2 <- KM_diff(df, tte.name = "time", event.name = "status", treat.name = "treat", weight.name = "w",
                  at_points = fit$at_points, alpha = 0.05, risk.points = NULL, modify_tau = FALSE,
                  draws.band = 500, seedstart = 99999, show_resamples = FALSE)
  ep <- ew_parts(d$g, df$treat)
  z0 <- ew_km_all(df, fit, fit2, r, ep$IF * 0, ep$dw)
  relmax <- function(a, b) { s <- b > 0; max(abs(a[s] - b[s]) / b[s]) }
  v0 <- fit$surv0^2 * colSums((df$w[z0$idx0] * z0$P0)^2); v1 <- fit$surv1^2 * colSums((df$w[z0$idx1] * z0$P1)^2)
  a5 <- max(relmax(v0, fit$sig2_surv0), relmax(v1, fit$sig2_surv1))
  say("G1.5a (dataset): S^2 sum w^2 psi^2 from P_a against sig2_surv0 / sig2_surv1: max relative difference %.3g (limit 1e-8)", a5)
  g5a <- a5 <= 1e-8
  b5 <- max(abs(z0$own0 - fit$surv0_star), abs(z0$own1 - fit$surv1_star))
  say("G1.5b: script's -S sum G w psi against surv0_star / surv1_star: max absolute difference %.3g (limit 1e-8)", b5)
  fin <- function(a, b) { s <- is.finite(a) | is.finite(b); if (!any(s)) 0 else max(abs(a[s] - b[s])) }
  c5 <- max(abs(z0$band$c_alpha_band - fit$c_alpha_band), fin(z0$band$sb_lower, fit$sb_lower), fin(z0$band$sb_upper, fit$sb_upper))
  say("G1.5c: IF = 0, script band against c_alpha_band, sb_lower, sb_upper: max absolute difference %.3g (limit 1e-8)", c5)
  rm <- call_rmst(df, fit, "w")
  d5 <- max(fin(z0$rmst$rmst_time_lower, rm$rmst_time_lower), fin(z0$rmst$rmst_time_upper, rm$rmst_time_upper),
            fin(z0$rmst$rmst_time_sb_lower, rm$rmst_time_sb_lower), fin(z0$rmst$rmst_time_sb_upper, rm$rmst_time_sb_upper))
  say("G1.5d: IF = 0, script RMST pointwise and band limits against cumulative_rmst_bands(): max absolute difference %.3g (limit 1e-8)", d5)
  say("Part D gates: 5a(dataset) %s, 5b %s, 5c %s, 5d %s", if (g5a) "PASS" else "FAIL", if (b5 <= 1e-8) "PASS" else "FAIL",
      if (c5 <= 1e-8) "PASS" else "FAIL", if (d5 <= 1e-8) "PASS" else "FAIL")
  if (length(fails)) { say("G1 FAIL: %s", paste(fails, collapse = ", ")); quit(status = 1) }
  say("G1 items 1 and 2: PASS")
  quit(status = 0)
}

# ---- --check: reproducibility of replicates 1-10 of cells 1 and 5 ----
if (isTRUE(opt$check)) {
  okall <- TRUE
  for (k in c(1, 5)) {
    f <- file.path(.here, opt$out, sprintf("%s_cell%02d.rds", VERSION, k))
    if (!file.exists(f)) { cat("missing", f, "\n"); quit(status = 1) }
    st <- readRDS(f)
    now <- run_reps(k, 1:10, 1)
    a <- st$res[1:10, , drop = FALSE]; b <- now$res
    same <- identical(colnames(a), colnames(b)) && identical(is.na(a), is.na(b)) &&
      all(abs(a - b)[!is.na(a)] <= 1e-10) && identical(unname(st$chk[1:10]), unname(now$chk))
    cat(sprintf("check cell %d: %s (max abs difference %.3g)\n", k, if (same) "OK" else "DIFFERENT",
                if (identical(dim(a), dim(b))) max(c(0, abs(a - b)), na.rm = TRUE) else NA))
    okall <- okall && same
  }
  quit(status = if (okall) 0 else 1)
}

# ---- run ----
cell_ids <- eval(parse(text = paste0("c(", opt$cells, ")")))
nrep <- as.integer(opt$reps); cores <- as.integer(opt$cores)
outdir <- file.path(.here, opt$out); dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("oc_sim.R version=%s (weightedsurv %s, %s) cells=%s reps=%d cores=%d out=%s  started %s\n", VERSION, PKG_VERSION,
            normalizePath(opt$pkg), opt$cells, nrep, cores, opt$out, format(Sys.time())))
t_all <- proc.time()[["elapsed"]]
for (k in cell_ids) {
  t0 <- proc.time()[["elapsed"]]
  rr <- run_reps(k, seq_len(nrep), cores)
  secs <- proc.time()[["elapsed"]] - t0
  saveRDS(c(rr, list(cell_id = k, version = VERSION, pkg_version = PKG_VERSION, secs = secs, cores = cores)),
          file.path(outdir, sprintf("%s_cell%02d.rds", VERSION, k)))
  failcols <- grep("\\.fail$", colnames(rr$res), value = TRUE)
  nf <- colSums(rr$res[, failcols, drop = FALSE] == 1, na.rm = TRUE)
  cat(sprintf("  cell %2d %-12s %d replicates %7.1f s  failures: %s\n", k, cells$cell[k], nrep, secs,
              if (any(nf > 0)) paste(sprintf("%s=%d", names(nf)[nf > 0], nf[nf > 0]), collapse = " ") else "none"))
}
cat(sprintf("oc_sim.R version=%s done: wall-clock %.1f s (%.1f min)\n", VERSION, proc.time()[["elapsed"]] - t_all,
            (proc.time()[["elapsed"]] - t_all) / 60))
