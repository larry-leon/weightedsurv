#!/usr/bin/env Rscript
# m6_summarise.R -- tables and required checks of M6 from m6_raw/cellNN.rds (TASK_km_fixes_stage2d.md, section 6).
#   Rscript dev/audit/stage2d/m6_summarise.R  -> m6_results.csv, m6_tables.md, m6_checks.txt (exit 1 if a required check fails)
if (!file.exists("DESCRIPTION")) stop("Run from the repo root.")
S <- file.path("dev", "audit", "stage2d")
source(file.path(S, "m6_common.R"))
oc <- utils::read.csv(file.path("dev", "audit", "oc", "oc_summary.csv"), stringsAsFactors = FALSE)
ocv <- function(k, q, f) { r <- oc[oc$cell_id == k & oc$quantity == q & oc$estimator == "proposed", ]; if (nrow(r) == 1) r[[f]] else NA_real_ }
rows <- list(); l4 <- list(); chk <- character(0); fails <- character(0)
say <- function(...) chk <<- c(chk, sprintf(...))
for (i in seq_len(nrow(m6_cells))) {
  cl <- m6_cells[i, ]; k <- cl$k
  x <- readRDS(file.path(S, "m6_raw", sprintf("cell%02d.rds", k)))
  R <- x$res; n <- nrow(R)
  for (g in c("fh00", "fh01")) {
    f <- function(nm) R[, paste0(g, ".", nm)]
    ok <- f("fail") == 0
    b <- f("bhat"); se <- f("se"); z <- f("z")
    cov <- f("lo") <= cl$truth & cl$truth <= f("hi")
    s5 <- seq_len(min(500, n)); st <- !is.na(f("star_se")[s5])
    scov <- (f("star_lo") <= cl$truth & cl$truth <= f("star_hi"))[s5]
    row <- data.frame(k = k, cell = cl$cell, scheme = if (g == "fh00") "FH(0,0)" else "FH(0,1)", n_rep = n, n_fail = sum(!ok),
      n_nonfinite_z = sum(ok & !is.finite(z)), truth = cl$truth,
      bias = if (cl$strata) NA else mean(b, na.rm = TRUE) - cl$truth, bias_mcse = if (cl$strata) NA else stats::sd(b, na.rm = TRUE) / sqrt(sum(!is.na(b))),
      emp_sd = if (cl$strata) NA else stats::sd(b, na.rm = TRUE), mean_se = if (cl$strata) NA else mean(se, na.rm = TRUE),
      se_sd = if (cl$strata) NA else mean(se, na.rm = TRUE) / stats::sd(b, na.rm = TRUE),
      coverage = if (cl$strata) NA else mean(cov, na.rm = TRUE), mean_length = if (cl$strata) NA else mean(f("hi") - f("lo"), na.rm = TRUE),
      star_n = if (g == "fh00" && !cl$strata) sum(st) else NA,
      star_se_sd = if (g == "fh00" && !cl$strata) mean(f("star_se")[s5], na.rm = TRUE) / stats::sd(b[s5], na.rm = TRUE) else NA,
      star_coverage = if (g == "fh00" && !cl$strata) mean(scov, na.rm = TRUE) else NA,
      star_mean_length = if (g == "fh00" && !cl$strata) mean((f("star_hi") - f("star_lo"))[s5], na.rm = TRUE) else NA,
      z_sd = stats::sd(z, na.rm = TRUE), reject = if (cl$null) mean(abs(z) > stats::qnorm(0.975), na.rm = TRUE) else NA,
      oc_bias = ocv(k, paste0("loghr_", g), "bias"), oc_se_sd = ocv(k, paste0("loghr_", g), "se_sd"),
      oc_coverage = ocv(k, paste0("loghr_", g), "coverage"), oc_mean_length = ocv(k, paste0("loghr_", g), "mean_length"),
      oc_reject = if (cl$null) ocv(k, paste0("logrank_", g), "reject") else NA, stringsAsFactors = FALSE)
    rows[[length(rows) + 1]] <- row
    # required: cox_robust agreement where there is no tied time (cells 2, 5, 8, 10, 14, 16)
    if (k %in% c(2, 5, 8, 10, 14, 16)) {
      nt <- ok & f("ties") == 0
      m <- c(se = max(f("chk_se")[nt]), z = max(f("chk_z")[nt]), zcox = max(f("chk_zcox")[nt]))
      pass <- all(is.finite(m)) && all(m < 1e-8) && sum(nt) > 0
      say("REQUIRED cell %2d %-12s %s: %d replicates without tied times (%d with); max relative difference from cox_robust(): sig_bhat_asy %.3g, |z| df_counting %.3g, |z| cox_rhogamma %.3g (limit 1e-8): %s",
          k, cl$cell, row$scheme, sum(nt), sum(ok & f("ties") == 1), m["se"], m["z"], m["zcox"], if (pass) "PASS" else "FAIL")
      if (!pass) fails <- c(fails, sprintf("cox_robust cell %d %s", k, g))
    }
    if (k %in% c(12, 21) && g == "fh00") {
      m <- max(f("cph_se")[ok]); pass <- is.finite(m) && m < 1e-6 && all(ok)
      say("REQUIRED cell %2d %-12s FH(0,0): %d replicates (%d with tied times); max relative difference of sig_bhat_asy from coxph(robust = TRUE, ties = 'breslow') %.3g (limit 1e-6): %s",
          k, cl$cell, sum(ok), sum(f("ties") == 1), m, if (pass) "PASS" else "FAIL")
      if (!pass) fails <- c(fails, sprintf("coxph cell %d", k))
    }
    if (g == "fh00" && k != 1) {
      say("evidence cell %2d %-18s FH(0,0): max relative difference from coxph robust: z^2 against rscore %.3g%s", k, cl$cell, max(f("cph_z")[ok]),
          if (cl$strata) " (stratified)" else sprintf(", sig_bhat_asy %.3g", max(f("cph_se")[ok])))
    }
    if (k %in% c(5, 14)) {
      sc <- f("l4_score"); dh <- f("l4_dhat")
      l4[[length(l4) + 1]] <- data.frame(k = k, cell = cl$cell, scheme = row$scheme,
        mean_cor = mean(f("l4_cor")), emp_cor = stats::cor(sc, dh),
        mean_sig2_score = mean(f("l4_sig2")), emp_var_score = stats::var(sc),
        mean_sig2_dhat = mean(f("l4_sig2_dhat")), emp_var_dhat = stats::var(dh),
        mean_cov = mean(f("l4_cov")), emp_cov = stats::cov(sc, dh),
        max_rel_sig2_dhat_vs_KM_diff = max(abs(f("l4_sig2_dhat") / f("km_sig2_dhat") - 1)), max_abs_cor = max(abs(f("l4_cor"))))
    }
  }
}
res <- do.call(rbind, rows); L4 <- do.call(rbind, l4)
utils::write.csv(res, file.path(S, "m6_results.csv"), row.names = FALSE)
utils::write.csv(L4, file.path(S, "m6_l4_results.csv"), row.names = FALSE)
f3 <- function(x, d = 3) ifelse(is.na(x), "–", formatC(x, digits = d, format = "f"))
md <- c("| k | cell | scheme | bias (MCSE) | SE/SD | coverage | mean length | resampling, first 500: SE/SD, coverage, length | reject (null cells) | evaluation `proposed`: bias, SE/SD, coverage, length, reject |",
        "|---|---|---|---|---|---|---|---|---|---|")
for (i in seq_len(nrow(res))) { r <- res[i, ]
  md <- c(md, sprintf("| %d | %s | %s | %s | %s | %s | %s | %s | %s | %s |", r$k, r$cell, r$scheme,
    if (is.na(r$bias)) "–" else sprintf("%s (%s)", f3(r$bias, 4), f3(r$bias_mcse, 4)), f3(r$se_sd), f3(r$coverage), f3(r$mean_length),
    if (is.na(r$star_se_sd)) "–" else sprintf("%s, %s, %s", f3(r$star_se_sd), f3(r$star_coverage), f3(r$star_mean_length)),
    f3(r$reject), if (is.na(r$oc_se_sd) && is.na(r$oc_reject)) "not run" else
      sprintf("%s, %s, %s, %s, %s", f3(r$oc_bias, 4), f3(r$oc_se_sd), f3(r$oc_coverage), f3(r$oc_mean_length), f3(r$oc_reject)))) }
md <- c(md, "", "| k | cell | scheme | correlation: mean estimated / empirical | var(score): mean estimated / empirical | var(dhat at 18): mean estimated / empirical | covariance: mean estimated / empirical | max rel. diff. of sig2_dhat from KM_diff() | max abs. cor |",
        "|---|---|---|---|---|---|---|---|---|")
for (i in seq_len(nrow(L4))) { r <- L4[i, ]
  md <- c(md, sprintf("| %d | %s | %s | %.3f / %.3f | %.1f / %.1f | %.5f / %.5f | %.4f / %.4f | %.2g | %.3f |", r$k, r$cell, r$scheme, r$mean_cor, r$emp_cor,
                      r$mean_sig2_score, r$emp_var_score, r$mean_sig2_dhat, r$emp_var_dhat, r$mean_cov, r$emp_cov, r$max_rel_sig2_dhat_vs_KM_diff, r$max_abs_cor)) }
writeLines(md, file.path(S, "m6_tables.md"))
chk <- c(chk, if (length(fails)) paste("M6 required checks: FAIL:", paste(fails, collapse = "; ")) else "M6 required checks: PASS")
writeLines(chk, file.path(S, "m6_checks.txt"))
cat(md, sep = "\n"); cat("\n"); cat(chk, sep = "\n")
quit(status = if (length(fails)) 1 else 0)
