#!/usr/bin/env Rscript
# =============================================================================
# km_m1_montecarlo.R -- M1: calibration of the weighted KM standard error
# (TASK_km_audit.md section 3, family M1; run only because the kickoff said "M1: GO")
#
# Budget: 6 cells x 2,000 replicates, n = 200, <= 5 min (stop and record if > 10 min).
# Data   : one arm, T ~ Exp(0.08), C ~ U(0, 30); truth S(t) = exp(-0.08 t).
# Cells  : ties {continuous; times rounded UP to 0.5} x weights {none; w ~ U(0.3,3)
#          independent of T; the same w normalised to sum 1}.
# Times  : t = 5, 10, 15.
# SEs    : package (calculate_risk_event_counts + KM_estimates);
#          survfit(robust = TRUE) (with weights when weighted);
#          documented formula S(t) * sqrt(sum_{s<=t} sum_i w_i^2 dN_i(s) / Y_w(s)^2).
# Output : dev/audit/km_m1_results.csv (unless --no-write). Seeded and deterministic.
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
NO_WRITE <- "--no-write" %in% args
if (!file.exists("DESCRIPTION")) stop("Run from the repo root.")
suppressMessages(pkgload::load_all(".", quiet = TRUE, export_all = FALSE, helpers = FALSE))
t0 <- proc.time()[["elapsed"]]
R <- 2000; n <- 200; TT <- c(5, 10, 15); truth <- exp(-0.08 * TT)

doc <- function(time, status, w, at) {
  ut <- sort(unique(time[status == 1]))
  Yw <- colSums(outer(time, ut, ">=") * w)
  E <- outer(time, ut, "==") * (status == 1)
  dNw <- colSums(E * w); dN2 <- colSums(E * w^2)
  S <- cumprod(1 - dNw / Yw); v <- cumsum(dN2 / Yw^2)
  idx <- findInterval(at, ut)
  Sa <- c(1, S)[idx + 1]; list(S = Sa, se = Sa * sqrt(c(0, v)[idx + 1]))
}
out <- list()
for (ties in c("continuous", "ceiling 0.5")) for (wt in c("none", "U(0.3,3)", "U(0.3,3)/sum")) {
  est <- array(NA_real_, c(R, 3, 3), dimnames = list(NULL, c("package", "survfit robust", "documented"), TT))
  se <- est
  for (r in seq_len(R)) {
    set.seed(31000 + r)                       # same T, C and w across cells within replicate r
    Tt <- stats::rexp(n, 0.08); C <- stats::runif(n, 0, 30); u <- stats::runif(n, 0.3, 3)
    time <- pmin(Tt, C); st <- as.numeric(Tt <= C)
    if (ties == "ceiling 0.5") time <- ceiling(time * 2) / 2
    w <- switch(wt, "none" = rep(1, n), "U(0.3,3)" = u, "U(0.3,3)/sum" = u / sum(u))
    at <- sort(unique(c(time, TT))); m <- match(TT, at)
    re <- calculate_risk_event_counts(time, st, w, at)
    km <- KM_estimates(re$ybar, re$nbar, re$sig2w_multiplier)
    est[r, "package", ] <- km$S_KM[m]; se[r, "package", ] <- sqrt(km$sig2_KM[m])
    y <- survival::Surv(time, st)
    fit <- if (wt == "none") survival::survfit(y ~ 1, robust = TRUE) else survival::survfit(y ~ 1, weights = w, robust = TRUE)
    sm <- summary(fit, times = TT, extend = TRUE)
    est[r, "survfit robust", ] <- sm$surv; se[r, "survfit robust", ] <- sm$std.err
    dd <- doc(time, st, w, TT)
    est[r, "documented", ] <- dd$S; se[r, "documented", ] <- dd$se
  }
  for (e in dimnames(est)[[2]]) for (k in seq_along(TT)) {
    S <- est[, e, k]; s <- se[, e, k]
    out[[length(out) + 1]] <- data.frame(ties = ties, weights = wt, estimator = e, t = TT[k], truth = round(truth[k], 6),
      mean_S = mean(S), emp_sd = stats::sd(S), mean_se = mean(s), se_over_sd = mean(s) / stats::sd(S),
      coverage = mean(abs(S - truth[k]) <= 1.96 * s), stringsAsFactors = FALSE)
  }
  el <- proc.time()[["elapsed"]] - t0
  message(sprintf("cell ties=%s weights=%s done, %.0fs elapsed", ties, wt, el))
  if (el > 600) { message("M1 exceeded 10 min; stopping"); break }
}
res <- do.call(rbind, out)
num <- c("mean_S", "emp_sd", "mean_se", "se_over_sd", "coverage"); res[num] <- lapply(res[num], round, 4)
if (!NO_WRITE) utils::write.csv(res, file.path("dev", "audit", "km_m1_results.csv"), row.names = FALSE)
options(width = 200); print(res, row.names = FALSE)
cat(sprintf("\nM1: %d replicates per cell, n = %d; elapsed %.0fs\n", R, n, proc.time()[["elapsed"]] - t0))
