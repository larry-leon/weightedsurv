#!/usr/bin/env Rscript
# m7_req1_diagnosis.R -- M7 requirement 1, the replicates in which the package under "estimated" and the
# evaluation's script differ by more than 1e-8 in the Cox / log-rank part, or in which the script's alignment
# confirmation does not pass (found in raw/m7_cellNN.rds), and one ordinary replicate for comparison.
# For each: whether the time fix (survival::aeqSurv, as df_counting() applies it) changed a time; the package,
# the script (cox_robust() + the section-7 adjustment, no handling of tied times), the tie-aware reference of the
# tests (ws_score_resid_ref() + ws_score_ew_ref()) and the coxph() residual routes (Breslow) on the time-fixed data.
# Also counts, over all 2,500 reproduction replicates, those in which the time fix changes a time.
#   Rscript dev/audit/stage2e/m7_req1_diagnosis.R   -> dev/audit/stage2e/m7_req1_diagnosis.txt
suppressPackageStartupMessages(library(survival))
source("dev/audit/stage2e/oc_common.R")
suppressMessages(pkgload::load_all(".", quiet = TRUE, export_all = TRUE, helpers = FALSE))
source("tests/testthat/helper-score-resid.R"); source("tests/testthat/helper-ew-ref.R")
q <- function(expr) { r <- NULL; invisible(utils::capture.output(r <- suppressWarnings(expr))); r }
sink("dev/audit/stage2e/m7_req1_diagnosis.txt", split = TRUE)
KR <- list(); tab <- NULL
for (k in c(7, 8, 9, 15, 19)) {
  x <- readRDS(sprintf("dev/audit/stage2e/raw/m7_cell%02d.rds", k))$res
  cx <- pmax(x[, "req1.cox0.se"], x[, "req1.cox0.z"], x[, "req1.cox1.se"], x[, "req1.cox1.z"])
  al <- x[, "req1.aligned0"] == 1 & x[, "req1.aligned1"] == 1
  nfix <- vapply(seq_len(nrow(x)), function(r) { df <- make_data(k, r)$df; sum(aeqSurv(Surv(df$time, df$status))[, "time"] != df$time) }, numeric(1))
  km <- pmax(x[, "req1.km"], x[, "req1.band"], x[, "req1.rmst"])
  tab <- rbind(tab, data.frame(cell = k, replicates = nrow(x), time_fix_changes_a_time = sum(nfix > 0), alignment_not_passed = sum(!al, na.rm = TRUE),
                               cox_differs_alignment_passed = sum(al & cx > 1e-8, na.rm = TRUE), cox_differs_and_no_time_changed = sum(al & cx > 1e-8 & nfix == 0, na.rm = TRUE),
                               largest_cox_difference_no_time_changed = max(cx[al & nfix == 0], na.rm = TRUE), largest_km_band_rmst_difference = max(km, na.rm = TRUE)))
  for (r in which((al & cx > 1e-8) | !al | is.na(al))) KR[[length(KR) + 1]] <- c(k, r)
}
cat("Requirement 1 by cell (replicates 1-500):\n"); print(tab, row.names = FALSE, digits = 3)
KR[[length(KR) + 1]] <- c(8, 1)
for (kr in KR) {
  k <- kr[1]; r <- kr[2]
  d <- make_data(k, r); df <- d$df; pw <- ps_weights(model = d$g)
  tf <- aeqSurv(Surv(df$time, df$status)); tfix <- tf[, "time"]
  o <- order(df$time); st <- df$time[o]
  gaps <- diff(st); nclose <- sum(gaps / st[-1] < 1e-7)
  cat(sprintf("\ncell %d replicate %d: n = %d; raw times distinct: %s; time-fixed times with a tie: %d; raw pairs closer than 1e-7 relative: %d; rows whose time aeqSurv changed: %d\n",
              k, r, nrow(df), !anyDuplicated(df$time), sum(duplicated(tfix)), nclose, sum(tfix != df$time)))
  if (nclose) { i <- which(gaps / st[-1] < 1e-7)[1]; cat(sprintf("   closest pair: %.12f and %.12f (difference %.3g), status %d and %d, arm %d and %d\n", st[i], st[i + 1], gaps[i], df$status[o][i], df$status[o][i + 1], df$treat[o][i], df$treat[o][i + 1])) }
  for (g in c(0, 1)) {
    sp <- list(rho = 0, gamma = g)
    dk <- q(df_counting(df, "time", "status", "treat", weight.name = "w", scheme = "fh", scheme_params = sp))
    de <- q(df_counting(df, "time", "status", "treat", weight.name = "w", scheme = "fh", scheme_params = sp, weight_variance = "estimated", weight_model = pw))
    ck <- q(cox_rhogamma(dk, scheme = "fh", scheme_params = sp)); ce <- q(cox_rhogamma(de, scheme = "fh", scheme_params = sp))
    ep <- ew_parts(d$g, df$treat); ord <- order(df$time)
    cr <- cox_robust(ck$fit); cr0 <- cox_robust(ck$fit, 0)
    ewpsi <- function(c0) c0$psi + as.numeric(ep$IF[ord, , drop = FALSE] %*% colSums((c0$psi / ck$fit$w_hat) * ep$dw[ord, , drop = FALSE]))
    se_script <- sqrt(sum(ewpsi(cr)^2)) / cr$info; z_script <- cr0$U / sqrt(sum(ewpsi(cr0)^2))
    # the Stage 2d / 2e test reference (handles tied times), on the package's time-fixed data
    rb <- ws_score_resid_ref(de$time, de$delta, de$z, de$w_hat, K = ce$fit$wt_rg, beta = ce$fit$bhat); r0 <- ws_score_resid_ref(de$time, de$delta, de$z, de$w_hat, K = ce$fit$wt_rg, beta = 0)
    se_ref <- sqrt(sum(ws_score_ew_ref(rb$psi, de$w_hat, de$ew_parts$IF, de$ew_parts$dw)^2)) / rb$info
    z_ref <- -r0$U / sqrt(sum(ws_score_ew_ref(r0$psi, de$w_hat, de$ew_parts$IF, de$ew_parts$dw)^2))
    cat(sprintf("   FH(0,%d) estimated: SE package %.10f, script %.10f (rel %.2e), tie-aware reference %.10f (rel %.2e); |z| package %.8f, script %.8f (rel %.2e), tie-aware reference %.8f (rel %.2e)\n",
                g, ce$fit$sig_bhat_asy, se_script, abs(ce$fit$sig_bhat_asy / se_script - 1), se_ref, abs(ce$fit$sig_bhat_asy / se_ref - 1),
                abs(ce$z.score), abs(z_script), abs(abs(ce$z.score) / abs(z_script) - 1), abs(z_ref), abs(abs(ce$z.score) / abs(z_ref) - 1)))
    cat(sprintf("            known:     SE package %.10f, script cox_robust %.10f (rel %.2e); |z| package %.8f, script %.8f (rel %.2e)\n",
                ck$fit$sig_bhat_asy, cr$se, abs(ck$fit$sig_bhat_asy / cr$se - 1), abs(ck$z.score), abs(cr0$z), abs(abs(ck$z.score) / abs(cr0$z) - 1)))
    if (g == 0) {
      df$tfix <- tfix; df$sfix <- tf[, "status"]
      adj <- function(rr) { rr <- as.numeric(rr); as.numeric(rr + pw$IF %*% colSums((rr / df$w) * pw$dw)) }
      f0 <- coxph(Surv(tfix, sfix) ~ treat, data = df, weights = w, init = 0, iter.max = 0, robust = TRUE, ties = "breslow", timefix = FALSE)
      rs <- as.numeric(residuals(f0, type = "score", weighted = TRUE))
      fb <- coxph(Surv(tfix, sfix) ~ treat, data = df, weights = w, robust = TRUE, ties = "breslow", timefix = FALSE)
      se_b <- sqrt(sum(adj(residuals(fb, type = "dfbeta", weighted = TRUE))^2))
      cat(sprintf("            coxph residual routes (Breslow): z^2 package / coxph - 1 = %.2e; SE package / coxph - 1 = %.2e\n", ce$z.score^2 / (sum(rs)^2 / sum(adj(rs)^2)) - 1, ce$fit$sig_bhat_asy / se_b - 1))
    }
  }
}
sink()
