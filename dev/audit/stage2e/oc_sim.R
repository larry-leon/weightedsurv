#!/usr/bin/env Rscript
# oc_sim.R -- Stage 2e, M7 (TASK_km_fixes_stage2e.md, section 6): Monte Carlo with the package code.
# Adapted copy of dev/audit/oc/oc_sim.R. Kept from it, unchanged: the time grid and summaries (named, band_stats,
# band_fail, km_fields, rmst_fields, cox_fields), the exact calls (call_km, call_rmst, call_dc, call_cox; the two
# new arguments are passed on), and the script-side estimated-weight code of the evaluation (psi_mat,
# km_multipliers, rmst_from_stars, dband_from_stars, ew_km_all), which requirement 1 compares the package with.
# Changed: one process, the package at HEAD only; the cells are those of section 6 (oc_common.R, table m7);
# every analysis is run with weight_variance = "known" and "estimated" on the same data and weight column.
#
#   Rscript dev/audit/stage2e/oc_sim.R --pkg=<path> [--cells=7,8,...] [--reps=<n, default the cell's>] [--cores=12]
#                                      [--out=raw] [--pilot]
#
#   --pilot   times replicates on one core (10 of cells 31, 36, 38, 40; 3 of every other cell) and prints the
#             projected wall-clock of the full run on --cores cores; nothing is stored
#
# Per-replicate results: <out>/m7_cellNN.rds, a list with `res` (replicates x fields numeric matrix), `chk`
# (the per-replicate data checksum) and timing. Field names are <estimator>.<quantity>.<field>[.<month>], with
# estimator = known / estimated (A, B, C in cell 40), plus req1.*, req2.*, req3.*, req4.* (the largest
# difference of each requirement in that replicate), ess.ratio0 / ess.ratio1, wd.* and ps.fail.

suppressPackageStartupMessages(library(survival))
.args <- commandArgs(FALSE)
.here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", .args, value = TRUE)[1])))
source(file.path(.here, "oc_common.R"))
opt <- parse_args(list(pkg = NA, cells = paste(m7$k, collapse = ","), reps = NA, cores = "12", out = "raw", pilot = FALSE))
stopifnot(!is.na(opt$pkg))
suppressMessages(pkgload::load_all(opt$pkg, quiet = TRUE, export_all = TRUE, helpers = FALSE, attach_testthat = FALSE))
PKG_VERSION <- as.character(utils::packageVersion("weightedsurv"))
PKG_SHA <- tryCatch(system(sprintf("git -C %s rev-parse HEAD", shQuote(opt$pkg)), intern = TRUE), error = function(e) NA_character_)
grDevices::pdf(NULL)                                   # no plot file, whatever a call draws
Q <- qnorm(0.975)
NT <- length(tt)

named <- function(prefix, x) stats::setNames(as.numeric(x), paste0(prefix, ".", if (length(x) == NT) tt else seq_along(x)))
na_t <- function(prefix) named(prefix, rep(NA_real_, NT))

# Band summary over the grid (as in the evaluation): a grid point counts only where the band has positive
# width. cov_in: truth inside at every such point in [3, 42]; cov_any: at every such point.
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

# ---- exact calls (the evaluation's, with the two new arguments) ----
call_km <- function(df, wn, r, draws = 0, draws.band = 500, wv = "known", pw = NULL)
  KM_diff(df, tte.name = "time", event.name = "status", treat.name = "treat", weight.name = wn, alpha = 0.05,
          draws = draws, draws.band = draws.band, seedstart = 8316951 + r, show_resamples = FALSE,
          weight_variance = wv, weight_model = pw)
call_rmst <- function(df, fit, wn, wv = "known", pw = NULL, ...)
  cumulative_rmst_bands(df, fit, tte.name = "time", event.name = "status", treat.name = "treat", weight.name = wn,
                        draws_sb = 500, plot = FALSE, weight_variance = wv, weight_model = pw, ...)
call_dc <- function(df, wn, g, wv = "known", pw = NULL, strata = NULL)
  df_counting(df, tte.name = "time", event.name = "status", treat.name = "treat", weight.name = wn, strata.name = strata,
              scheme = "fh", scheme_params = list(rho = 0, gamma = g), weight_variance = wv, weight_model = pw)
call_cox <- function(dc, g, draws = 0)
  cox_rhogamma(dc, scheme = "fh", scheme_params = list(rho = 0, gamma = g), draws = draws)

# ---- blocks (as in the evaluation; the truth is taken by scenario) ----
km_fields <- function(est, fit, jj, sc, band = TRUE) {
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
                         delta_sc(fit$at_points, sc), fit$c_alpha_band))
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

# ---- the evaluation's script-side estimated-weight code (dev/audit/oc/oc_sim.R, unchanged) ----
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
dband_from_stars <- function(dhat, dhat_star, sig2_dhat, alpha = 0.05) {
  Zdhat_star <- dhat_star / sqrt(sig2_dhat)
  sups <- apply(abs(Zdhat_star), 2, max, na.rm = TRUE)
  c_alpha_band <- quantile(sups, 1 - alpha)
  list(c_alpha_band = c_alpha_band, sb_lower = dhat - c_alpha_band * sqrt(sig2_dhat),
       sb_upper = dhat + c_alpha_band * sqrt(sig2_dhat))
}
ew_km_all <- function(df, fit, fit2, r, IF, dw, B = 500) {
  tf <- aeqSurv(Surv(df$time, df$status))
  time_fix <- tf[, "time"]; delta_fix <- tf[, "status"]
  idx0 <- which(df$treat == 0); idx1 <- which(df$treat == 1)
  P0 <- psi_mat(time_fix[idx0], delta_fix[idx0], df$w[idx0], fit$at_points)
  P1 <- psi_mat(time_fix[idx1], delta_fix[idx1], df$w[idx1], fit$at_points)
  ek <- ew_km(P0, P1, fit$surv0, fit$surv1, idx0, idx1, df$w, IF, dw)
  G <- km_multipliers(time_fix, df$treat, 8316951 + r, B)
  G2 <- km_multipliers(time_fix, df$treat, 99999, B)
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

# ---- M7 helpers ----
# largest difference relative to the largest reference value (the definition of Gate G2, check 6); finite entries
relv <- function(a, b) {
  s <- is.finite(a) | is.finite(b)
  if (!any(s)) return(0)
  if (!identical(is.finite(a), is.finite(b))) return(Inf)
  max(abs(a[s] - b[s])) / max(abs(b[s]))
}
same <- function(a, b) as.numeric(identical(a, b))

# One option ("known" or "estimated") on one dataset and weight column: the evaluation's analyses.
# Returns the fields (prefix `est`) and the objects that the requirement checks need.
run_option <- function(est, df, wn, wv, pw, r, cl, jj0 = NULL) {
  pwm <- if (wv == "estimated") pw else NULL
  out <- numeric(0); keep <- list()
  fit <- tryCatch(call_km(df, wn, r, draws.band = if (cl$band) 500 else 0, wv = wv, pw = pwm), error = function(e) e)
  jj <- if (inherits(fit, "error")) rep(NA_integer_, NT) else last_le(fit$at_points, tt)
  out <- c(out, km_fields(est, fit, jj, cl$sc, band = cl$band))
  keep$fit <- fit; keep$jj <- jj
  if (cl$rmst) {
    if (inherits(fit, "error")) out <- c(out, rmst_fields(est, fit, jj, NULL))
    else {
      rm <- tryCatch(call_rmst(df, fit, wn, wv = wv, pw = pwm), error = function(e) e)
      out <- c(out, rmst_fields(est, rm, jj, R_sc(fit$at_points, cl$sc)))
      keep$rm <- rm
    }
  }
  for (g in if (cl$fh01) c(0, 1) else 0) {
    res <- tryCatch({
      dc <- suppressWarnings(call_dc(df, wn, g, wv = wv, pw = pwm))
      cx <- suppressWarnings(call_cox(dc, g))
      if (is.null(cx$fit) || !is.finite(cx$fit$bhat)) stop("root finding failed")
      list(dc = dc, cx = cx)
    }, error = function(e) e)
    if (inherits(res, "error")) out <- c(out, cox_fields(est, g, fail = 1))
    else out <- c(out, cox_fields(est, g, res$cx$fit$bhat, res$cx$fit$sig_bhat_asy, log(res$cx$hr_ci_asy$lower),
                                  log(res$cx$hr_ci_asy$upper), res$dc$z.score))
    keep[[paste0("cox", g)]] <- res
  }
  list(out = out, keep = keep)
}
# identical estimates under the two options (requirement 3): 1 = identical (or not evaluated in the cell)
req3_fields <- function(a, b, cl) {
  ok <- function(x, y, f) if (inherits(x, "error") || inherits(y, "error") || is.null(x) || is.null(y)) NA_real_ else f(x, y)
  v <- c(dhat = ok(a$fit, b$fit, function(x, y) same(x$dhat, y$dhat)),
         rmst = if (cl$rmst) ok(a$rm, b$rm, function(x, y) same(x$rmst_time, y$rmst_time)) else NA_real_,
         bhat0 = ok(a$cox0, b$cox0, function(x, y) same(x$cx$fit$bhat, y$cx$fit$bhat)),
         lr0 = ok(a$cox0, b$cox0, function(x, y) same(x$dc$lr, y$dc$lr)),
         bhat1 = if (cl$fh01) ok(a$cox1, b$cox1, function(x, y) same(x$cx$fit$bhat, y$cx$fit$bhat)) else NA_real_,
         lr1 = if (cl$fh01) ok(a$cox1, b$cox1, function(x, y) same(x$dc$lr, y$dc$lr)) else NA_real_)
  stats::setNames(v, paste0("req3.", names(v)))
}

# ---- one replicate ----
one_rep <- function(k, r) {
  cl <- m7[m7$k == k, ]
  d <- m7_data(k, r); df <- d$df
  tmax <- c(max(df$time[df$treat == 0], -Inf), max(df$time[df$treat == 1], -Inf))
  out <- named("evaluable", as.numeric(min(tmax) >= tt))
  fin <- function(v) { attr(v, "chk") <- d$chk; v }

  # ---- cell 40: analyses of the subgroup X > 0 ----
  if (cl$sub) {
    sg <- df$X > 0; dsub <- df[sg, ]
    tmax <- c(max(dsub$time[dsub$treat == 0], -Inf), max(dsub$time[dsub$treat == 1], -Inf))
    out <- named("evaluable", as.numeric(min(tmax) >= tt))
    pf <- tryCatch(ps_weights(treat ~ X, data = df, type = cl$type), error = function(e) e)       # model on all 800 rows
    ps <- tryCatch(ps_weights(treat ~ X, data = dsub, type = cl$type), error = function(e) e)     # model refitted on the subgroup
    out <- c(out, ps.fail = as.numeric(inherits(pf, "error") || inherits(ps, "error")))
    if (out[["ps.fail"]] == 1) return(fin(out))
    dA <- dsub; dA$w <- pf$weights[sg]
    dB <- dsub; dB$w <- ps$weights
    A <- run_option("A", dA, "w", "known", NULL, r, cl)        # full-sample weights, "known"
    B <- run_option("B", dB, "w", "estimated", ps, r, cl)      # refitted weights, "estimated"
    C <- run_option("C", dB, "w", "known", NULL, r, cl)        # refitted weights, "known"
    ess <- function(dd) { e <- suppressWarnings(call_dc(dd, "w", 0))$ess; c(e$ratio[2], e$ratio[1]) }
    eA <- tryCatch(ess(dA), error = function(e) c(NA, NA)); eB <- tryCatch(ess(dB), error = function(e) c(NA, NA))
    out <- c(out, A$out, B$out, C$out, req3_fields(C$keep, B$keep, cl),
             stats::setNames(c(eA, eB), c("A.ess.ratio0", "A.ess.ratio1", "ess.ratio0", "ess.ratio1")))
    return(fin(out))
  }

  pw <- tryCatch(if (!is.null(d$g)) ps_weights(model = d$g, type = cl$type) else ps_weights(treat ~ X, data = df, type = cl$type),
                 error = function(e) e)
  out <- c(out, ps.fail = as.numeric(inherits(pw, "error")))
  if (inherits(pw, "error")) return(fin(out))
  if (!cl$repro) df$w <- pw$weights            # cells 7, 8, 9, 15, 19 keep the evaluation's own weight column

  # ---- cell 33: stratified log-rank only ----
  if (cl$strat) {
    for (g in c(0, 1)) {
      zk <- tryCatch(suppressWarnings(call_dc(df, "w", g, strata = "s")), error = function(e) e)
      ze <- tryCatch(suppressWarnings(call_dc(df, "w", g, wv = "estimated", pw = pw, strata = "s")), error = function(e) e)
      zz <- function(x) if (inherits(x, "error")) NA_real_ else x$z.score_stratified
      out <- c(out, stats::setNames(c(zz(zk), zz(ze), as.numeric(inherits(zk, "error")), as.numeric(inherits(ze, "error"))),
                                    paste0(c("known", "estimated", "known", "estimated"), ".slr", g, c(".z", ".z", ".fail", ".fail"))),
               stats::setNames(if (inherits(zk, "error") || inherits(ze, "error")) NA_real_ else same(zk$lr_stratified, ze$lr_stratified),
                               paste0("req3.lrs", g)))
      if (g == 0 && !inherits(zk, "error")) out <- c(out, ess.ratio0 = zk$ess$ratio[2], ess.ratio1 = zk$ess$ratio[1])
    }
    return(fin(out))
  }

  K <- run_option("known", df, "w", "known", NULL, r, cl)
  E <- run_option("estimated", df, "w", "estimated", pw, r, cl)
  out <- c(out, K$out, E$out, req3_fields(K$keep, E$keep, cl))
  dc0 <- K$keep$cox0
  out <- c(out, if (inherits(dc0, "error")) c(ess.ratio0 = NA_real_, ess.ratio1 = NA_real_)
                else c(ess.ratio0 = dc0$dc$ess$ratio[2], ess.ratio1 = dc0$dc$ess$ratio[1]))

  # ---- requirement 1 (cells 7, 8, 9, 15, 19): the package under "estimated" against the evaluation's script ----
  if (cl$repro) {
    v <- c(km = NA_real_, band = NA_real_, rmst = NA_real_, cox0.se = NA_real_, cox0.z = NA_real_, cox1.se = NA_real_, cox1.z = NA_real_,
           aligned0 = NA_real_, aligned1 = NA_real_)
    fit <- K$keep$fit; fe <- E$keep$fit; jj <- K$keep$jj
    ep <- ew_parts(d$g, df$treat)
    sw <- tryCatch({
      fit2 <- KM_diff(df, tte.name = "time", event.name = "status", treat.name = "treat", weight.name = "w",
                      at_points = fit$at_points, alpha = 0.05, risk.points = NULL, modify_tau = FALSE,
                      draws.band = 500, seedstart = 99999, show_resamples = FALSE)
      ew_km_all(df, fit, fit2, r, ep$IF, ep$dw)
    }, error = function(e) e)
    if (!inherits(sw, "error") && !inherits(fe, "error")) {
      ok <- !is.na(jj)
      v["km"] <- max(relv(fe$sig2_surv0[jj[ok]], sw$ek$sig2_surv0[jj[ok]]), relv(fe$sig2_surv1[jj[ok]], sw$ek$sig2_surv1[jj[ok]]),
                     relv(fe$sig2_dhat[jj[ok]], sw$ek$sig2_dhat[jj[ok]]))
      v["band"] <- max(abs(fe$c_alpha_band / sw$band$c_alpha_band - 1), relv(fe$sb_lower, sw$band$sb_lower), relv(fe$sb_upper, sw$band$sb_upper))
      rme <- E$keep$rm
      if (!is.null(rme) && !inherits(rme, "error")) {
        v["rmst"] <- max(relv(rme$sig2_rmst_time, sw$rmst$sig2_rmst_time), abs(rme$c_alpha_band / sw$rmst$c_alpha_band - 1),
                         relv(rme$rmst_time_sb_lower, sw$rmst$rmst_time_sb_lower), relv(rme$rmst_time_sb_upper, sw$rmst$rmst_time_sb_upper))
      }
    }
    for (g in c(0, 1)) {
      rk <- K$keep[[paste0("cox", g)]]; re <- E$keep[[paste0("cox", g)]]
      if (inherits(rk, "error") || inherits(re, "error")) next
      cx <- rk$cx
      cr <- cox_robust(cx$fit); cr0 <- cox_robust(cx$fit, 0)
      ord <- order(df$time)
      aligned <- isTRUE(all.equal(as.numeric(cx$fit$time), df$time[ord])) && all(cx$fit$z == df$treat[ord])
      v[paste0("aligned", g)] <- as.numeric(aligned)
      if (aligned) {
        ewpsi <- function(c0) c0$psi + as.numeric(ep$IF[ord, , drop = FALSE] %*% colSums((c0$psi / cx$fit$w_hat) * ep$dw[ord, , drop = FALSE]))
        se_ew <- sqrt(sum(ewpsi(cr)^2)) / cr$info
        z_ew <- cr0$U / sqrt(sum(ewpsi(cr0)^2))
        v[paste0("cox", g, ".se")] <- abs(re$cx$fit$sig_bhat_asy / se_ew - 1)
        v[paste0("cox", g, ".z")] <- max(abs(abs(re$cx$z.score) / abs(z_ew) - 1), abs(abs(re$dc$z.score) / abs(z_ew) - 1))
      }
    }
    out <- c(out, stats::setNames(v, paste0("req1.", names(v))))
  }

  # ---- requirement 2 and wlr_dhat_estimates() (cells 31, 32): "stabilized" against "ipw" from the same model ----
  if (cl$both) {
    pwi <- ps_weights(model = d$g, type = "ipw"); di <- df; di$w <- pwi$weights
    for (wv in c("known", "estimated")) {
      st <- if (wv == "known") K$keep else E$keep
      v <- c(sig2d = NA_real_, calpha = NA_real_, rmstvar = NA_real_)
      fi <- tryCatch(call_km(di, "w", r, wv = wv, pw = if (wv == "estimated") pwi else NULL), error = function(e) e)
      if (!inherits(fi, "error") && !inherits(st$fit, "error")) {
        v["sig2d"] <- relv(st$fit$sig2_dhat, fi$sig2_dhat)
        v["calpha"] <- abs(st$fit$c_alpha_band / fi$c_alpha_band - 1)
        ri <- tryCatch(call_rmst(di, fi, "w", wv = wv, pw = if (wv == "estimated") pwi else NULL), error = function(e) e)
        if (!inherits(ri, "error") && !inherits(st$rm, "error")) v["rmstvar"] <- relv(st$rm$sig2_rmst_time, ri$sig2_rmst_time)
      }
      out <- c(out, stats::setNames(v, paste0("req2.", wv, ".", names(v))))
    }
    e0 <- E$keep$cox0
    wd <- if (inherits(e0, "error")) NULL else tryCatch(wlr_dhat_estimates(e0$dc, tzero = 18), error = function(e) NULL)
    out <- c(out, if (is.null(wd)) c(wd.score = NA_real_, wd.dhat = NA_real_, wd.cor = NA_real_, wd.cov = NA_real_, wd.sig2s = NA_real_, wd.sig2d = NA_real_)
                  else c(wd.score = wd$score, wd.dhat = wd$dhat, wd.cor = wd$cor_wlr_dhat, wd.cov = wd$cov_wlr_dhat, wd.sig2s = wd$sig2.score, wd.sig2d = wd$sig2_dhat))
  }

  # ---- requirement 4 (cell 34, tied times): FH(0,0) under "estimated" against the coxph residual routes (Breslow) ----
  if (cl$coxph) {
    v <- c(z = NA_real_, se = NA_real_)
    e0 <- E$keep$cox0
    if (!inherits(e0, "error")) {
      adj <- function(rr) { rr <- as.numeric(rr); as.numeric(rr + pw$IF %*% colSums((rr / df$w) * pw$dw)) }
      f0 <- coxph(Surv(time, status) ~ treat, data = df, weights = w, init = 0, iter.max = 0, robust = TRUE, ties = "breslow")
      rs <- as.numeric(residuals(f0, type = "score", weighted = TRUE))
      fb <- coxph(Surv(time, status) ~ treat, data = df, weights = w, robust = TRUE, ties = "breslow")
      se_b <- sqrt(sum(adj(residuals(fb, type = "dfbeta", weighted = TRUE))^2))
      v["z"] <- abs(e0$dc$z.score^2 / (sum(rs)^2 / sum(adj(rs)^2)) - 1)
      v["se"] <- abs(e0$cx$fit$sig_bhat_asy / se_b - 1)
    }
    out <- c(out, stats::setNames(v, paste0("req4.", names(v))))
  }
  fin(out)
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
  allnm <- unique(unlist(lapply(L, names)))            # a replicate whose ps_weights() failed has fewer fields
  res <- t(vapply(L, function(x) as.numeric(x[allnm]), numeric(length(allnm))))
  colnames(res) <- allnm; rownames(res) <- reps
  list(res = res, chk = vapply(L, function(x) attr(x, "chk"), character(1)))
}

cell_ids <- eval(parse(text = paste0("c(", opt$cells, ")")))
stopifnot(all(cell_ids %in% m7$k))
cores <- as.integer(opt$cores)

# ---- pilot ----
if (isTRUE(opt$pilot)) {
  cat(sprintf("M7 pilot, one core (weightedsurv %s at %s)  %s\n", PKG_VERSION, substr(PKG_SHA, 1, 7), format(Sys.time())))
  per <- numeric(0)
  for (k in m7$k) {
    nr <- if (k %in% c(31, 36, 38, 40)) 10 else 3
    t0 <- proc.time()[["elapsed"]]
    invisible(run_reps(k, seq_len(nr), 1))
    per[as.character(k)] <- (proc.time()[["elapsed"]] - t0) / nr
    cat(sprintf("  cell %2d %-20s %2d replicates: %6.2f s per replicate\n", k, m7$cell[m7$k == k], nr, per[as.character(k)]))
  }
  tot <- sum(per * m7$reps)
  cat(sprintf("simulation: %.0f core-seconds; projected wall-clock on %d cores: %.1f min\n", tot, cores, tot / cores / 60))
  saveRDS(list(per = per, total_core_secs = tot, cores = cores), file.path(.here, "m7_pilot_sim.rds"))
  quit(status = 0)
}

# ---- run ----
outdir <- file.path(.here, opt$out); dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("oc_sim.R (M7) weightedsurv %s at %s, cells=%s cores=%d out=%s  started %s\n", PKG_VERSION, substr(PKG_SHA, 1, 7),
            opt$cells, cores, opt$out, format(Sys.time())))
t_all <- proc.time()[["elapsed"]]
for (k in cell_ids) {
  nrep <- if (is.na(opt$reps)) m7$reps[m7$k == k] else as.integer(opt$reps)
  t0 <- proc.time()[["elapsed"]]
  rr <- run_reps(k, seq_len(nrep), cores)
  secs <- proc.time()[["elapsed"]] - t0
  saveRDS(c(rr, list(cell_id = k, pkg_version = PKG_VERSION, pkg_sha = PKG_SHA, secs = secs, cores = cores)),
          file.path(outdir, sprintf("m7_cell%02d.rds", k)))
  failcols <- grep("\\.fail$", colnames(rr$res), value = TRUE)
  nf <- colSums(rr$res[, failcols, drop = FALSE] == 1, na.rm = TRUE)
  cat(sprintf("  cell %2d %-20s %d replicates %7.1f s  failures: %s\n", k, m7$cell[m7$k == k], nrep, secs,
              if (any(nf > 0)) paste(sprintf("%s=%d", names(nf)[nf > 0], nf[nf > 0]), collapse = " ") else "none"))
}
cat(sprintf("oc_sim.R (M7) done: wall-clock %.1f s (%.1f min)\n", proc.time()[["elapsed"]] - t_all,
            (proc.time()[["elapsed"]] - t_all) / 60))
