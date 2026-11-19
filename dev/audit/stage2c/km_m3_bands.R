#!/usr/bin/env Rscript
# =============================================================================
# km_m3_bands.R -- M3: simultaneous-band coverage (TASK_km_fixes_stage2c.md, section 6)
#
#   Rscript dev/audit/stage2c/km_m3_bands.R run <version> <pkgdir> [reps]
#       version: "before" (the START_SHA worktree) or "after" (HEAD); pkgdir: package source to load
#       -> dev/audit/stage2c/km_m3_<version>.rds (per-replicate records and a digest of every output)
#   Rscript dev/audit/stage2c/km_m3_bands.R summary
#       -> km_m3_results.csv and the table on stdout; checks that cell U is identical in every replicate
#
# Design: two arms, N = 800, C ~ U(0, 40), X ~ N(0, 1), T ~ Exp(0.08 exp(0.7 X - 0.4 A)).
#   U: A ~ Bernoulli(0.5), no weights;  D: A ~ Bernoulli(expit(-0.5 + X)), true propensity weights;
#   E: as D, weights from glm(A ~ X) in each replicate;  F: A ~ Bernoulli(expit(-0.5 + 2 X)), true weights.
# Replicate r of cell k (k = 1..4 for U, D, E, F) uses set.seed(200000 * k + r).
# Per replicate: fit <- KM_diff(df, ..., weight.name, draws.band = 500, show_resamples = FALSE) on the
# default grid; rm <- cumulative_rmst_bands(df, fit, ..., weight.name, draws_sb = 500, plot = FALSE).
# Truth: delta(t) = S_1(t) - S_0(t), S_a(t) = E_X[exp(-0.08 exp(0.7 X - 0.4 a) t)], on a grid of step
# 0.01 by numerical integration over X in [-12, 12]; R(t) = integral_0^t delta by the trapezoid rule on that grid.
# Read-off: the estimate at t is the value at the last grid point <= t (right-continuous step), compared
# with delta(t); the cumulative RMST at grid point t_j is compared with R(t_j).
# Simultaneous coverage is recorded over all grid points ("all") and over the grid points where the
# standard error is positive ("se>0"): where it is 0 (before the first event, or after S reaches 0)
# the band has width 0 and cannot cover a true value that is not exactly the estimate.
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[1] else "summary"
if (!file.exists("DESCRIPTION")) stop("Run from the repo root.")
OUT <- file.path("dev", "audit", "stage2c")
CELLS <- c("U", "D", "E", "F")
TT <- c(5, 10, 20)
expit <- function(x) 1 / (1 + exp(-x))

gen <- function(k, r) {
  set.seed(200000 * k + r)
  cell <- CELLS[k]; N <- 800
  X <- stats::rnorm(N)
  e <- switch(cell, U = rep(0.5, N), D = , E = expit(-0.5 + X), F = expit(-0.5 + 2 * X))
  A <- stats::rbinom(N, 1, e)
  Tt <- stats::rexp(N, 0.08 * exp(0.7 * X - 0.4 * A)); C <- stats::runif(N, 0, 40)
  if (cell == "E") e <- stats::fitted(stats::glm(A ~ X, family = stats::binomial))
  data.frame(time = pmin(Tt, C), status = as.numeric(Tt <= C), treat = A, w = ifelse(A == 1, 1 / e, 1 / (1 - e)))
}

truth <- local({
  g <- seq(0, 40, by = 0.01)
  Sa <- function(a) vapply(g, function(t) stats::integrate(function(x) exp(-0.08 * exp(0.7 * x - 0.4 * a) * t) * stats::dnorm(x),
                                                          -12, 12, rel.tol = 1e-10)$value, 0)
  d <- Sa(1) - Sa(0)
  R <- c(0, cumsum((d[-1] + d[-length(d)]) / 2 * diff(g)))
  list(g = g, delta = function(t) stats::approx(g, d, t, rule = 2)$y, R = function(t) stats::approx(g, R, t, rule = 2)$y)
})

one_rep <- function(k, r) {
  df <- gen(k, r)
  wn <- if (CELLS[k] == "U") NULL else "w"
  fit <- KM_diff(df, "time", "status", "treat", weight.name = wn, draws.band = 500, show_resamples = FALSE)
  rm <- cumulative_rmst_bands(df, fit, "time", "status", "treat", weight.name = wn, draws_sb = 500, plot = FALSE)
  at <- fit$at_points
  se <- sqrt(fit$sig2_dhat)
  dtrue <- truth$delta(at)
  i_t <- findInterval(TT, at)
  win <- at >= 2 & at <= 25
  cov_sb <- fit$sb_lower <= dtrue & dtrue <= fit$sb_upper
  pos <- se > 0
  rse <- sqrt(rm$sig2_rmst_time)
  Rtrue <- truth$R(at)
  cov_rsb <- rm$rmst_time_sb_lower <= Rtrue & Rtrue <= rm$rmst_time_sb_upper
  rpos <- rse > 0
  j20 <- findInterval(20, at)
  c(cell_k = k, r = r,
    setNames(as.numeric(fit$lower[i_t] <= truth$delta(TT) & truth$delta(TT) <= fit$upper[i_t]), paste0("pw_d_", TT)),
    sb_all = all(cov_sb), sb_all_pos = all(cov_sb[pos]), sb_2_25 = all(cov_sb[win]), sb_2_25_pos = all(cov_sb[win & pos]),
    rmst_pw_20 = as.numeric(rm$rmst_time_lower[j20] <= Rtrue[j20] & Rtrue[j20] <= rm$rmst_time_upper[j20]),
    rmst_sb_all = all(cov_rsb), rmst_sb_all_pos = all(cov_rsb[rpos]), rmst_sb_2_25 = all(cov_rsb[win]), rmst_sb_2_25_pos = all(cov_rsb[win & rpos]),
    c_alpha = unname(fit$c_alpha_band), rmst_c_alpha = unname(rm$c_alpha_band),
    hw_sb_10 = unname(fit$c_alpha_band) * se[i_t[2]], se_dhat_10 = se[i_t[2]], dhat_10 = fit$dhat[i_t[2]],
    n_se0 = sum(!pos), n_rse0 = sum(!rpos)) -> rec
  list(rec = rec, digest = digest::digest(list(fit, rm[setdiff(names(rm), "mean_rule")])))
}

if (mode == "run") {
  version <- args[2]; pkgdir <- args[3]; R_REPS <- if (length(args) >= 4) as.integer(args[4]) else 1000L
  stopifnot(version %in% c("before", "after"))
  suppressMessages(pkgload::load_all(pkgdir, quiet = TRUE, export_all = FALSE, helpers = FALSE))
  grDevices::pdf(NULL)
  CORES <- max(1L, parallel::detectCores() - 2L)
  sha <- system(sprintf("git -C '%s' rev-parse HEAD", pkgdir), intern = TRUE)
  cat(sprintf("M3 %s: package from %s (%s); %d cores; %d replicates per cell; started %s\n", version, pkgdir, sha, CORES, R_REPS, format(Sys.time())))
  t0 <- proc.time()[["elapsed"]]; out <- list()
  for (k in seq_along(CELLS)) {
    tk <- proc.time()[["elapsed"]]
    res <- parallel::mclapply(seq_len(R_REPS), function(r) one_rep(k, r), mc.cores = CORES, mc.preschedule = TRUE)
    bad <- vapply(res, inherits, logical(1), "try-error")
    if (any(bad)) stop("cell ", CELLS[k], ": ", as.character(res[[which(bad)[1]]]))
    out[[k]] <- list(rec = do.call(rbind, lapply(res, `[[`, "rec")), digest = vapply(res, `[[`, "", "digest"))
    el <- proc.time()[["elapsed"]] - t0
    cat(sprintf("cell %s: %d replicates in %.0f s (total %.1f min)\n", CELLS[k], R_REPS, proc.time()[["elapsed"]] - tk, el / 60))
    if (el > 60 * 60 && k < length(CELLS)) { cat("STOP: wall-clock passed 60 min\n"); break }
  }
  saveRDS(list(version = version, sha = sha, cells = CELLS[seq_along(out)], out = out,
               minutes = (proc.time()[["elapsed"]] - t0) / 60), file.path(OUT, sprintf("km_m3_%s.rds", version)))
  cat(sprintf("wall-clock: %.1f min\n", (proc.time()[["elapsed"]] - t0) / 60))
  quit(status = 0)
}

if (mode == "summary") {
  B <- readRDS(file.path(OUT, "km_m3_before.rds")); A <- readRDS(file.path(OUT, "km_m3_after.rds"))
  cat(sprintf("before: %s (%.1f min); after: %s (%.1f min)\n", B$sha, B$minutes, A$sha, A$minutes))
  rows <- list()
  for (v in list(B, A)) for (k in seq_along(v$out)) {
    m <- as.data.frame(v$out[[k]]$rec)
    rows[[length(rows) + 1]] <- data.frame(cell = CELLS[k], version = v$version, n_reps = nrow(m),
      pw_d_5 = mean(m$pw_d_5), pw_d_10 = mean(m$pw_d_10), pw_d_20 = mean(m$pw_d_20),
      sb_all = mean(m$sb_all), sb_all_se_pos = mean(m$sb_all_pos), sb_2_25 = mean(m$sb_2_25), sb_2_25_se_pos = mean(m$sb_2_25_pos),
      rmst_pw_20 = mean(m$rmst_pw_20), rmst_sb_all = mean(m$rmst_sb_all), rmst_sb_all_se_pos = mean(m$rmst_sb_all_pos),
      rmst_sb_2_25 = mean(m$rmst_sb_2_25), rmst_sb_2_25_se_pos = mean(m$rmst_sb_2_25_pos),
      mean_c_alpha = mean(m$c_alpha), mean_rmst_c_alpha = mean(m$rmst_c_alpha), mean_hw_sb_10 = mean(m$hw_sb_10),
      se_over_sd_dhat_10 = mean(m$se_dhat_10) / stats::sd(m$dhat_10),
      mean_n_se0 = mean(m$n_se0), mean_n_rmst_se0 = mean(m$n_rse0), stringsAsFactors = FALSE)
  }
  res <- do.call(rbind, rows)
  res <- res[order(match(res$cell, CELLS), match(res$version, c("before", "after"))), ]
  utils::write.csv(res, file.path(OUT, "km_m3_results.csv"), row.names = FALSE)
  options(width = 250); print(res, row.names = FALSE, digits = 4)
  kU <- which(CELLS == "U")
  same_dig <- identical(B$out[[kU]]$digest, A$out[[kU]]$digest)
  same_rec <- identical(B$out[[kU]]$rec, A$out[[kU]]$rec)
  cat(sprintf("\ncell U: before and after identical in every replicate: digests of all KM_diff and cumulative_rmst_bands output %s (%d replicates); records %s\n",
              if (same_dig) "IDENTICAL" else "DIFFER", length(A$out[[kU]]$digest), if (same_rec) "IDENTICAL" else "DIFFER"))
  for (k in setdiff(seq_along(CELLS), kU)) cat(sprintf("cell %s: replicates whose output differs before/after: %d of %d\n", CELLS[k],
                                                      sum(B$out[[k]]$digest != A$out[[k]]$digest), length(A$out[[k]]$digest)))
  quit(status = if (same_dig && same_rec) 0 else 1)
}
stop("usage: km_m3_bands.R run <before|after> <pkgdir> [reps] | summary")
