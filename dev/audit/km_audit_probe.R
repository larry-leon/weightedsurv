#!/usr/bin/env Rscript
# =============================================================================
# km_audit_probe.R  --  seed probes for TASK_weightedsurv_km_audit.md
#
# Purpose : re-run, against the LOCAL working tree, every finding that was
#           confirmed in chat on GitHub HEAD edd61a2 (2026-02-12).
# Run     : from the repo root
#               Rscript dev/audit/km_audit_probe.R
#           --installed  test an installed copy instead of the tree
#           --no-write   print the table only; leave probe_results.csv alone
# Output  : a table on stdout and dev/audit/probe_results.csv
#           (the P5 row holds a timing, so that one row differs run to run)
# Status  : REPRODUCED      defect is present in this tree
#           NOT_REPRODUCED  defect is absent in this tree
#           PASS / FAIL     control checks (things that should agree with survival)
#           INFO            measurement only
#           PROBE_ERROR     the probe itself failed -- read the 'observed' column
#
# survival is deliberately NOT attached (probe R6 depends on that).  Every
# survival call below is namespaced.  This script changes nothing in R/.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if ("--installed" %in% args) {
  library(weightedsurv)
  src <- paste0("INSTALLED weightedsurv ", utils::packageVersion("weightedsurv"))
} else {
  if (!file.exists("DESCRIPTION")) stop("Run from the repo root (DESCRIPTION not found).")
  if (!requireNamespace("pkgload", quietly = TRUE)) stop("pkgload is required to load the local tree (or pass --installed).")
  pkgload::load_all(".", quiet = TRUE, export_all = FALSE, helpers = FALSE)
  src <- "LOCAL TREE via pkgload::load_all"
}
out_dir <- if (dir.exists("dev")) file.path("dev", "audit") else "."
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
grDevices::pdf(NULL)   # swallow any plots

# ---- helpers ---------------------------------------------------------------
rows <- list()
add <- function(id, area, probe, observed, status) {
  rows[[length(rows) + 1]] <<- data.frame(id = id, area = area, probe = probe,
                                          observed = observed, status = status, stringsAsFactors = FALSE)
}
# run expr; return list(value, error, warnings)
run <- function(expr) {
  ws <- character(0)
  val <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning")
  }), error = function(e) structure(conditionMessage(e), class = "probe_error"))
  list(value = val, error = if (inherits(val, "probe_error")) unclass(val) else NULL, warnings = ws)
}
probe <- function(id, area, title, code) {
  r <- tryCatch(code(), error = function(e) list(observed = paste("probe failed:", conditionMessage(e)), status = "PROBE_ERROR"))
  add(id, area, title, r$observed, r$status)
}
rep_if <- function(cond) if (isTRUE(cond)) "REPRODUCED" else "NOT_REPRODUCED"
f <- function(x, d = 4) formatC(x, digits = d, format = "f")

sim <- function(n, seed, ties = FALSE, cens = 30) {
  set.seed(seed)
  d <- data.frame(treat = rbinom(n, 1, 0.5))
  d$time <- rexp(n, 0.08 * exp(-0.3 * d$treat)); cc <- runif(n, 0, cens)
  d$status <- as.numeric(d$time <= cc); d$time <- pmin(d$time, cc)
  if (ties) d$time <- ceiling(d$time * 2) / 2
  d$w <- runif(n, 0.3, 3)
  d
}
# survfit for one arm, evaluated at 'times'
sf_at <- function(d, g, times, w = NULL, robust = NULL) {
  dd <- d[d$treat == g, ]
  y <- survival::Surv(dd$time, dd$status)
  fit <- if (is.null(w)) survival::survfit(y ~ 1)
         else if (is.null(robust)) survival::survfit(y ~ 1, weights = dd[[w]])
         else survival::survfit(y ~ 1, weights = dd[[w]], robust = robust)
  summary(fit, times = times, extend = TRUE)
}
km1 <- function(U, D, W, at) {          # package single-arm KM + SE through the exported helpers
  re <- calculate_risk_event_counts(U, D, W, at)
  km <- KM_estimates(re$ybar, re$nbar, re$sig2w_multiplier)
  list(S = km$S_KM, se = sqrt(km$sig2_KM))
}

# ---- controls: these should PASS ---------------------------------------------
probe("OK1", "control", "Unweighted S and Greenwood SE equal survfit (veteran; n=300 heavy ties)", function() {
  v <- survival::veteran; v$treat <- as.numeric(v$trt) - 1
  dd <- sim(300, 1, ties = TRUE)
  mx <- 0
  for (d in list(v, dd)) {
    kd <- KM_diff(d, "time", "status", "treat")
    for (g in 0:1) {
      ref <- sf_at(d, g, kd$at_points)
      S <- if (g == 0) kd$surv0 else kd$surv1; v2 <- if (g == 0) kd$sig2_surv0 else kd$sig2_surv1
      mx <- max(mx, abs(S - ref$surv), abs(sqrt(v2) - ref$std.err), na.rm = TRUE)
    }
  }
  list(observed = sprintf("max abs diff (S and SE) = %.1e", mx), status = if (mx < 1e-10) "PASS" else "FAIL")
})
probe("OK2", "control", "Weighted S equals survfit(weights=) (no ties and heavy ties)", function() {
  mx <- 0
  for (tie in c(FALSE, TRUE)) {
    d <- sim(300, 2, ties = tie); kd <- KM_diff(d, "time", "status", "treat", weight.name = "w")
    for (g in 0:1) mx <- max(mx, abs((if (g == 0) kd$surv0 else kd$surv1) - sf_at(d, g, kd$at_points, "w")$surv))
  }
  list(observed = sprintf("max abs diff S = %.1e", mx), status = if (mx < 1e-10) "PASS" else "FAIL")
})
probe("OK3", "control", "Weighted SE, natural-scale weights, no ties: within 10% of survfit robust where S>0.1", function() {
  d <- sim(300, 3); g <- d[d$treat == 1, ]; at <- sort(unique(g$time))
  k <- km1(g$time, g$status, g$w, at); ref <- sf_at(d, 1, at, "w", robust = TRUE)
  keep <- k$S > 0.1 & ref$std.err > 0
  rr <- range(k$se[keep] / ref$std.err[keep])
  list(observed = sprintf("pkg SE / survfit robust SE in [%.3f, %.3f]", rr[1], rr[2]),
       status = if (rr[1] > 0.9 && rr[2] < 1.1) "PASS" else "FAIL")
})

# ---- K: Kaplan-Meier estimation --------------------------------------------
probe("K1", "KM estimate", "KM_diff(modify_tau=TRUE): KM on a grid that omits event times (continuous times, n=60)", function() {
  d <- sim(60, 11)
  kd <- KM_diff(d, "time", "status", "treat", draws.band = 50, modify_tau = TRUE, show_resamples = FALSE)
  e0 <- max(abs(kd$surv0 - sf_at(d, 0, kd$at_points)$surv)); e1 <- max(abs(kd$surv1 - sf_at(d, 1, kd$at_points)$surv))
  list(observed = sprintf("max|S-survfit| arm0 %s, arm1 %s; min package S = %s", f(e0), f(e1), f(min(kd$surv0, kd$surv1))),
       status = rep_if(max(e0, e1) > 1e-6))
})
probe("K1b", "KM estimate", "KM_diff(at_points = seq(0,30,5)): user grid without event times", function() {
  d <- sim(120, 5); ap <- seq(0, 30, 5)
  kd <- KM_diff(d, "time", "status", "treat", at_points = ap)
  e <- max(abs(kd$surv0 - sf_at(d, 0, ap)$surv))
  list(observed = sprintf("max|S0-survfit| = %s", f(e)), status = rep_if(e > 1e-6))
})
probe("K2", "KM quantile", "Median when S(t) equals 0.5 on an interval (no censoring, n even 2..60)", function() {
  res <- t(sapply(seq(2, 60, by = 2), function(n) {
    km <- KM_estimates(n:1, 1:n)
    c(pkg = km_quantile(1:n, km$S_KM, sqrt(km$sig2_KM))$qhat,
      ref = unname(summary(survival::survfit(survival::Surv(1:n, rep(1, n)) ~ 1))$table["median"]))
  }))
  nbad <- sum(abs(res[, "pkg"] - res[, "ref"]) > 1e-8)
  list(observed = sprintf("%d of %d medians differ from survfit (package returns the left end, survfit the midpoint)", nbad, nrow(res)),
       status = rep_if(nbad > 0))
})
probe("K2b", "KM quantile", "Same case through df_counting() with default arguments: is anything signalled?", function() {
  dd <- data.frame(time = c(1:10, 1:10 + 0.5), status = 1, treat = rep(0:1, each = 10))
  r <- run(df_counting(dd, "time", "status", "treat", by.risk = 2))
  if (!is.null(r$error)) return(list(observed = paste("error:", r$error), status = "PROBE_ERROR"))
  q <- r$value$quantile_results; chk <- r$value$quantile_check
  differs <- abs(q$quantile[q$group == "control"] - chk$time[1]) > 1e-8
  qw <- grep("quantile", r$warnings, value = TRUE)
  list(observed = sprintf("control median: package %s, survfit %s; quantile warnings raised: %d",
                          f(q$quantile[q$group == "control"], 2), f(chk$time[1], 2), length(qw)),
       status = rep_if(differs && length(qw) == 0))
})
probe("K3", "KM variance", "Weighted SE is not invariant to rescaling the weights (w vs w/sum(w))", function() {
  d <- sim(300, 3); g <- d[d$treat == 1, ]; at <- sort(unique(g$time)); q <- round(quantile(seq_along(at), c(.25, .5, .75, .95)))
  a <- km1(g$time, g$status, g$w, at)$se[q]; b <- km1(g$time, g$status, g$w / sum(g$w), at)$se[q]
  list(observed = sprintf("SE(w/sum w)/SE(w) at 25/50/75/95%% of grid: %s", paste(f(b / a, 2), collapse = ", ")),
       status = rep_if(max(abs(b / a - 1)) > 0.01))
})
probe("K4", "KM variance", "Weighted SE with tied event times vs survfit robust", function() {
  d <- sim(300, 3, ties = TRUE); g <- d[d$treat == 1, ]; at <- sort(unique(g$time)); q <- round(quantile(seq_along(at), c(.2, .4, .6, .8)))
  a <- km1(g$time, g$status, g$w, at)$se[q]; ref <- sf_at(d, 1, at, "w", robust = TRUE)$std.err[q]
  list(observed = sprintf("pkg SE / survfit robust SE at 4 times: %s", paste(f(a / ref, 2), collapse = ", ")),
       status = rep_if(max(a / ref) > 1.1))
})
probe("K5", "KM variance", "All weights equal to 2 vs weights 2 with one perturbed by 1e-9 (branch discontinuity)", function() {
  d <- sim(300, 3); g <- d[d$treat == 1, ]; at <- sort(unique(g$time)); i <- round(length(at) / 2)
  one <- km1(g$time, g$status, rep(1, nrow(g)), at)$se[i]; two <- km1(g$time, g$status, rep(2, nrow(g)), at)$se[i]
  near <- km1(g$time, g$status, 2 + c(1e-9, rep(0, nrow(g) - 1)), at)$se[i]
  list(observed = sprintf("SE: w=1 %s | w=2 %s | w~2 %s", f(one), f(two), f(near)), status = rep_if(abs(two / near - 1) > 0.05))
})
probe("K6", "KM estimate", "Weighted KM whose last observation is an event: floating-point residue, negative S", function() {
  neg <- 0; R <- 400
  for (s in seq_len(R)) {
    set.seed(s); n <- 40; tt <- sort(rexp(n)); st <- rbinom(n, 1, .7); st[n] <- 1; w <- runif(n, .2, 4)
    if (min(km1(tt, st, w, tt)$S) < 0) neg <- neg + 1
  }
  # end to end: df_counting() then the plotting function with its defaults
  fail <- 0; tried <- 0
  for (s in 1:60) {
    d <- sim(80, 100 + s)
    for (g in 0:1) { i <- which(d$treat == g); d$status[i[which.max(d$time[i])]] <- 1 }
    dc <- run(df_counting(d, "time", "status", "treat", weight.name = "w"))
    if (!is.null(dc$error)) next
    tried <- tried + 1
    p <- run(KM_plot_2sample_weighted_counting(dc$value))
    if (!is.null(p$error) && grepl("outside", p$error)) fail <- fail + 1
  }
  list(observed = sprintf("min(S) < 0 in %d of %d single-arm sims; KM_plot_2sample_weighted_counting() stops in %d of %d datasets", neg, R, fail, tried),
       status = rep_if(neg > 0))
})
probe("K7", "KM estimate", "Unsorted at_points are accepted and give values outside [0,1]", function() {
  d <- sim(120, 5)
  r <- run(KM_diff(d, "time", "status", "treat", at_points = rev(sort(unique(d$time)))))
  if (!is.null(r$error)) return(list(observed = paste("error:", r$error), status = "NOT_REPRODUCED"))
  rg <- range(r$value$surv0, na.rm = TRUE)
  list(observed = sprintf("no error; S0 ranges from %s to %s", f(rg[1], 2), f(rg[2], 2)), status = rep_if(rg[1] < 0 || rg[2] > 1))
})

# ---- R: robustness and input handling ---------------------------------------
err_probe <- function(expr_fun, pattern) function() {
  r <- run(expr_fun())
  if (is.null(r$error)) return(list(observed = sprintf("no error (%d warnings)", length(r$warnings)), status = "NOT_REPRODUCED"))
  list(observed = paste("error:", r$error), status = rep_if(grepl(pattern, r$error)))
}
d5 <- sim(120, 5)
probe("R1", "robustness", "KM_diff() when one arm is empty", err_probe(function() KM_diff(transform(d5, treat = 1), "time", "status", "treat"), "sig2w_multiplier"))
probe("R2", "robustness", "df_counting() when one arm has no events", err_probe(function() {
  d0 <- d5; d0$status[d0$treat == 1] <- 0; df_counting(d0, "time", "status", "treat", get.cox = FALSE) }, "subscript out of bounds"))
probe("R3a", "robustness", "KM_diff() with an NA time", err_probe(function() { dn <- d5; dn$time[3] <- NA; KM_diff(dn, "time", "status", "treat") }, "TRUE/FALSE"))
probe("R3b", "robustness", "df_counting() with an NA time", err_probe(function() { dn <- d5; dn$time[3] <- NA; df_counting(dn, "time", "status", "treat") }, "TRUE/FALSE"))
probe("R3c", "robustness", "KM_diff() with an NA weight: silent NA survival", function() {
  dw <- d5; dw$w[5] <- NA
  r <- run(KM_diff(dw, "time", "status", "treat", weight.name = "w"))
  if (!is.null(r$error)) return(list(observed = paste("error:", r$error), status = "NOT_REPRODUCED"))
  na <- anyNA(c(r$value$surv0, r$value$surv1))
  list(observed = sprintf("no error, %d warnings, NA in survival: %s", length(r$warnings), na), status = rep_if(na && length(r$warnings) == 0))
})
probe("R4a", "robustness", "KM_diff() accepts negative weights", function() {
  dn <- d5; dn$w[1:5] <- -1
  r <- run(KM_diff(dn, "time", "status", "treat", weight.name = "w"))
  list(observed = if (is.null(r$error)) sprintf("no error, %d warnings", length(r$warnings)) else paste("error:", r$error),
       status = rep_if(is.null(r$error) && length(r$warnings) == 0))
})
probe("R4b", "robustness", "df_counting() with zero weights: false 'Discrepancy in KM curve fit'", function() {
  dz <- d5; dz$w[dz$treat == 0][1:10] <- 0
  r <- run(df_counting(dz, "time", "status", "treat", weight.name = "w", get.cox = FALSE))
  kd <- KM_diff(dz, "time", "status", "treat", weight.name = "w")
  g <- dz[dz$treat == 0 & dz$w > 0, ]
  ref <- summary(survival::survfit(survival::Surv(g$time, g$status) ~ 1, weights = g$w), times = kd$at_points, extend = TRUE)$surv
  km_ok <- max(abs(kd$surv0 - ref)) < 1e-10
  w <- grep("Discrepancy in KM", r$warnings, value = TRUE)
  list(observed = sprintf("KM itself correct: %s; discrepancy warnings: %d%s", km_ok, length(w), if (!is.null(r$error)) paste("; error:", r$error) else ""),
       status = rep_if(km_ok && length(w) > 0))
})
probe("R5", "robustness", "KM_diff(draws.band = 10) with default show_resamples", err_probe(function() KM_diff(d5, "time", "status", "treat", draws.band = 10), "subscript out of bounds"))
probe("R6", "robustness", "Stratified df_counting() when survival is not attached", function() {
  if ("package:survival" %in% search()) return(list(observed = "survival is attached in this session; run in a clean session", status = "PROBE_ERROR"))
  v <- survival::veteran; v$treat <- as.numeric(v$trt) - 1
  r <- run(df_counting(v, "time", "status", "treat", strata.name = "celltype"))
  list(observed = if (is.null(r$error)) "no error" else paste("error:", r$error), status = rep_if(!is.null(r$error) && grepl("strata", r$error)))
})
probe("R7", "robustness", "Event coded 1/2: df_counting() accepts, KM_diff() rejects", function() {
  d12 <- d5; d12$status <- d12$status + 1
  a <- run(df_counting(d12, "time", "status", "treat")); b <- run(KM_diff(d12, "time", "status", "treat"))
  list(observed = sprintf("df_counting: %s; KM_diff: %s", if (is.null(a$error)) "ok" else "error", if (is.null(b$error)) "ok" else "error"),
       status = rep_if(is.null(a$error) != is.null(b$error)))
})

# ---- I: inference built on the KM curves ------------------------------------
probe("I1", "inference", "Simultaneous-band critical value ignores alpha", function() {
  a <- KM_diff(d5, "time", "status", "treat", draws.band = 500, alpha = 0.05, show_resamples = FALSE)$c_alpha_band
  b <- KM_diff(d5, "time", "status", "treat", draws.band = 500, alpha = 0.20, show_resamples = FALSE)$c_alpha_band
  list(observed = sprintf("c_alpha: alpha=0.05 -> %s, alpha=0.20 -> %s", f(a), f(b)), status = rep_if(isTRUE(all.equal(unname(a), unname(b)))))
})
probe("I2a", "inference", "wlr_dhat_estimates(): pooled S taken one grid point after tzero", function() {
  dc <- df_counting(d5, "time", "status", "treat"); t0 <- 12
  got <- wlr_dhat_estimates(dc, tzero = t0)$sig2_dhat
  dN <- diff(c(0, dc$nbar0 + dc$nbar1)); yb <- dc$ybar0 + dc$ybar1
  dJ <- ifelse(yb == 1, 0, (dN - 1) / (yb - 1)); dL <- ifelse(yb == 0, 0, dN / yb)
  h2 <- ifelse(dc$ybar0 * dc$ybar1 > 0, yb / (dc$ybar0 * dc$ybar1), 0) * (dc$at_points <= t0)
  S0 <- dc$survP[max(which(dc$at_points <= t0))]
  want <- S0^2 * sum(h2 * (1 - dJ) * dL)
  list(observed = sprintf("sig2_dhat returned %.6f; with S_pooled(tzero) it is %.6f (ratio %s)", got, want, f(got / want)),
       status = rep_if(abs(got / want - 1) > 1e-8))
})
probe("I2b", "inference", "wlr_dhat_estimates() with tzero beyond the last time", function() {
  dc <- df_counting(d5, "time", "status", "treat")
  r <- run(wlr_dhat_estimates(dc, tzero = 1e6))
  if (!is.null(r$error)) return(list(observed = paste("error:", r$error), status = "NOT_REPRODUCED"))
  list(observed = sprintf("no error; length(dhat) = %d", length(r$value$dhat)), status = rep_if(length(r$value$dhat) == 0))
})
probe("I3", "inference", "cumulative_rmst_bands(): trapezoid rule on the KM step function vs exact area", function() {
  d <- sim(200, 21)
  kd <- KM_diff(d, "time", "status", "treat", draws.band = 200, show_resamples = FALSE)
  rm <- cumulative_rmst_bands(d, kd, "time", "status", "treat", draws_sb = 200, plot = FALSE)
  tau <- max(kd$at_points)
  ex <- sapply(0:1, function(g) { dd <- d[d$treat == g, ]
    unname(summary(survival::survfit(survival::Surv(dd$time, dd$status) ~ 1), rmean = tau)$table["rmean"]) })
  got <- utils::tail(rm$rmst_time, 1); want <- ex[2] - ex[1]
  list(observed = sprintf("tau=%s: package %s, exact step-function area %s", f(tau, 2), f(got), f(want)), status = rep_if(abs(got - want) > 1e-6))
})
probe("I4", "inference", "Quantile self-check uses survfit's 0.95 limits whatever conf_level is", function() {
  d <- sim(200, 9)
  a <- run(df_counting(d, "time", "status", "treat", conf_level = 0.95, verbose = TRUE))
  b <- run(df_counting(d, "time", "status", "treat", conf_level = 0.90, verbose = TRUE))
  na <- length(grep("quantile", a$warnings)); nb <- length(grep("quantile", b$warnings))
  list(observed = sprintf("quantile-discrepancy warnings: conf_level 0.95 -> %d, 0.90 -> %d", na, nb), status = rep_if(na == 0 && nb > 0))
})

probe("I5", "inference", "qprob != 0.5: quantile_results and the survfit cross-check use opposite conventions", function() {
  d <- sim(200, 9)
  r <- run(df_counting(d, "time", "status", "treat", qprob = 0.25, verbose = TRUE))
  if (!is.null(r$error)) return(list(observed = paste("error:", r$error), status = "PROBE_ERROR"))
  pk <- r$value$quantile_results$quantile[r$value$quantile_results$group == "control"]
  ck <- r$value$quantile_check$time[1]
  g <- d[d$treat == 0, ]; fit <- survival::survfit(survival::Surv(g$time, g$status) ~ 1)
  q25 <- unname(quantile(fit, 0.25)$quantile); q75 <- unname(quantile(fit, 0.75)$quantile)
  list(observed = sprintf("control, qprob=0.25: quantile_results %s (survfit probs=0.75 is %s); quantile_check %s (survfit probs=0.25 is %s)",
                          f(pk, 2), f(q75, 2), f(ck, 2), f(q25, 2)),
       status = rep_if(abs(pk - q75) < 1e-8 && abs(ck - q25) < 1e-8 && abs(pk - ck) > 1e-8))
})

# ---- P: package behaviour ----------------------------------------------------
probe("P2", "package", "Functions overwrite the caller's RNG state (set.seed inside)", function() {
  set.seed(123); invisible(KM_diff(d5, "time", "status", "treat", draws.band = 50, show_resamples = FALSE)); x1 <- runif(1)
  set.seed(999); invisible(KM_diff(d5, "time", "status", "treat", draws.band = 50, show_resamples = FALSE)); x2 <- runif(1)
  list(observed = sprintf("next runif after user seeds 123 / 999: %.6f / %.6f", x1, x2), status = rep_if(x1 == x2))
})
probe("P5", "package", "df_counting() cost at n = 2000 (O(n^2) matrices)", function() {
  d <- sim(2000, 7); invisible(gc(reset = TRUE))
  tm <- system.time(df_counting(d, "time", "status", "treat"))[["elapsed"]]
  list(observed = sprintf("elapsed %.1fs, peak R memory %.0f MB", tm, sum(gc()[, 6])), status = "INFO")
})

# ---- report ------------------------------------------------------------------
res <- do.call(rbind, rows)
no_write <- "--no-write" %in% args
if (!no_write) utils::write.csv(res, file.path(out_dir, "probe_results.csv"), row.names = FALSE)
cat("\nweightedsurv KM audit probe\n")
cat("source   :", src, "\n")
cat("R        :", R.version.string, "| survival", as.character(utils::packageVersion("survival")), "\n")
cat("survival attached:", "package:survival" %in% search(), "\n\n")
cat(sprintf("| %-4s | %-14s | %-11s | %s\n", "id", "status", "area", "probe -> observed"))
for (i in seq_len(nrow(res))) with(res[i, ], cat(sprintf("| %-4s | %-14s | %-11s | %s -> %s\n", id, status, area, probe, observed)))
cat("\n", paste(names(table(res$status)), table(res$status), sep = ": ", collapse = " | "), "\n", sep = "")
if (!no_write) cat("written:", file.path(out_dir, "probe_results.csv"), "\n")
if (any(res$status == "PROBE_ERROR")) quit(status = 1)
