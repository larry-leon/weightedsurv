#!/usr/bin/env Rscript
# =============================================================================
# vignette_option.R -- Step V, parts 2-4 of TASK_km_fixes_stage2e.md. No vignette is edited.
# Uses the data frames of the weighted examples as the vignettes build them (saved by
# vignette_impact.R knit, under HEAD: <knit_after>/data_<file>.rds) and the package at HEAD.
#   part 2: does ps_weights() with the vignette's formula and type = "stabilized" reproduce the weight column?
#   part 3: effective sample size by arm (df_counting()$ess)
#   part 4: "known" against "estimated": log-rank z and p, hazard ratio and interval (cox_rhogamma, FH(0,0)),
#           SE of the survival difference at the median follow-up time, band critical value, RMST difference with
#           its interval at the last time; and the weighted plotKM.band_subgroups() and cumulative_rmst_bands()
#           figures under the two options, side by side (vignette_known_vs_estimated.pdf).
#
#   Rscript dev/audit/stage2e/vignette_option.R <knit_after>
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
if (!file.exists("DESCRIPTION")) stop("Run from the repo root.")
suppressMessages(pkgload::load_all(".", quiet = TRUE, export_all = FALSE, helpers = FALSE))
S2 <- file.path("dev", "audit", "stage2e"); KA <- args[1]
q <- function(expr) { r <- NULL; invisible(utils::capture.output(r <- suppressWarnings(expr))); r }
F8 <- treat ~ age + meno + size + grade3 + nodes + pgr + chemo + er
F7 <- treat ~ age + meno + size + grade3 + nodes + pgr + er
EX <- list(
  list(id = "examples / Rotterdam", file = "vignettes/weightedsurv_examples.Rmd", line = 299, data = c("examples", "df_rotterdam"), formula = F8, band = TRUE),
  list(id = "examples / cross-trial (CTC)", file = "vignettes/weightedsurv_examples.Rmd", line = 444, data = c("examples", "df_CTC"), formula = F7, band = TRUE),
  list(id = "martingale / Rotterdam", file = "vignettes/articles/martingale_resampling_applications.Rmd", line = 174, data = c("martingale", "df_rotterdam"), formula = F8, band = TRUE),
  list(id = "coxmeth / Rotterdam", file = "vignettes/articles/weightedcox_methodology.Rmd", line = 329, data = c("coxmeth", "df_rotterdam"), formula = F8, band = FALSE))
sink(file.path(S2, "vignette_option.txt"), split = TRUE)
cat("Step V, parts 2-4. Package at HEAD (", system("git rev-parse --short HEAD", intern = TRUE), "). Data frames as the vignettes build them.\n", sep = "")
grDevices::pdf(file.path(S2, "vignette_known_vs_estimated.pdf"), width = 13, height = 10)
for (ex in EX) {
  d <- readRDS(file.path(KA, paste0("data_", ex$data[1], ".rds")))[[ex$data[2]]]
  cat(sprintf("\n==== %s (%s:%d), n = %d, %d treated ====\n", ex$id, ex$file, ex$line, nrow(d), sum(d$treat == 1)))
  pw <- ps_weights(ex$formula, data = d, type = "stabilized")
  cat(sprintf("part 2: ps_weights(%s, type = \"stabilized\") against the column sw.weights: largest absolute difference %.3g (weights range %.4g to %.4g)\n",
              paste(deparse(ex$formula), collapse = ""), max(abs(pw$weights - d$sw.weights)), min(d$sw.weights), max(d$sw.weights)))
  cat(sprintf("        rows of the model = rows of the data frame: %s; missing values in the model's variables: %d\n",
              length(pw$weights) == nrow(d), sum(!stats::complete.cases(d[, all.vars(ex$formula)]))))
  dfc <- function(wv, sp) q(df_counting(d, "tte", "event", "treat", arms = c("treat", "control"), by.risk = 24, weight.name = "sw.weights",
                                         scheme = "fh", scheme_params = sp, weight_variance = wv, weight_model = if (wv == "estimated") pw else NULL))
  P00 <- list(rho = 0, gamma = 0)
  dk <- dfc("known", P00); de <- dfc("estimated", P00)
  cat("part 3: effective sample size by arm (df_counting()$ess)\n"); print(dk$ess, row.names = FALSE, digits = 4)
  ck <- q(cox_rhogamma(dk, scheme = "fh", scheme_params = P00, draws = 0)); ce <- q(cox_rhogamma(de, scheme = "fh", scheme_params = P00, draws = 0))
  tm <- stats::median(d$tte)
  row <- function(lab, k, e) cat(sprintf("  %-58s %-28s %-28s\n", lab, k, e))
  cat("part 4: known against estimated\n"); row("", "known", "estimated")
  f4 <- function(x) formatC(x, digits = 4, format = "g")
  row("log-rank z (df_counting, FH(0,0))", f4(dk$z.score), f4(de$z.score))
  row("  one-sided p = 1 - pnorm(z)", f4(1 - stats::pnorm(dk$z.score)), f4(1 - stats::pnorm(de$z.score)))
  row("  logrank_results: chisq, two-sided p", sprintf("%s, %s", f4(dk$logrank_results$chisq), f4(dk$logrank_results$pval)), sprintf("%s, %s", f4(de$logrank_results$chisq), f4(de$logrank_results$pval)))
  ci <- function(cx) sprintf("%.4f (%.4f, %.4f)", cx$hr_ci_asy$hr, cx$hr_ci_asy$lower, cx$hr_ci_asy$upper)
  row("hazard ratio and 95% interval (cox_rhogamma, FH(0,0))", ci(ck), ci(ce))
  row("  SE of the log hazard ratio", f4(ck$fit$sig_bhat_asy), f4(ce$fit$sig_bhat_asy))
  row("  coxph: cox_results$cox_text; $se", sprintf("%s; %s", dk$cox_results$cox_text, f4(dk$cox_results$se)), sprintf("%s; %s", de$cox_results$cox_text, f4(de$cox_results$se)))
  stopifnot(identical(dk$lr, de$lr), identical(ck$fit$bhat, ce$fit$bhat), identical(dk$surv1, de$surv1))
  if (ex$band) {
    pk_args <- list(df = d, xlabel = "Months", ylabel = "Survival difference", yseq_length = 5, cex_Yaxis = 0.7, risk_cex = 0.7, tau_add = 42, by.risk = 12,
                    risk_delta = 0.075, risk.pad = 0.03, ymax.pad = 0.125, tte.name = "tte", treat.name = "treat", event.name = "event",
                    weight.name = "sw.weights", draws.band = 1000, qtau = 0.025, show_resamples = FALSE)
    grDevices::pdf(NULL)
    tk <- do.call(plotKM.band_subgroups, pk_args)
    te <- do.call(plotKM.band_subgroups, c(pk_args, list(weight_variance = "estimated", weight_model = pw)))
    rk <- cumulative_rmst_bands(df = d, fit = tk$fit_itt, tte.name = "tte", event.name = "event", treat.name = "treat", weight.name = "sw.weights", draws_sb = 1000, plot = FALSE)
    re <- cumulative_rmst_bands(df = d, fit = te$fit_itt, tte.name = "tte", event.name = "event", treat.name = "treat", weight.name = "sw.weights", draws_sb = 1000, plot = FALSE,
                                weight_variance = "estimated", weight_model = pw)
    invisible(grDevices::dev.off())
    stopifnot(identical(tk$fit_itt$dhat, te$fit_itt$dhat), identical(rk$rmst_time, re$rmst_time))
    j <- findInterval(tm, tk$fit_itt$at_points)
    row(sprintf("SE of the survival difference at the median time (%.1f)", tm), f4(sqrt(tk$fit_itt$sig2_dhat[j])), f4(sqrt(te$fit_itt$sig2_dhat[j])))
    row("  difference at that time", f4(tk$fit_itt$dhat[j]), f4(te$fit_itt$dhat[j]))
    row("band critical value (draws.band = 1000)", f4(tk$fit_itt$c_alpha_band), f4(te$fit_itt$c_alpha_band))
    row("  mean width of the band", f4(mean(tk$fit_itt$sb_upper - tk$fit_itt$sb_lower)), f4(mean(te$fit_itt$sb_upper - te$fit_itt$sb_lower)))
    rr <- function(r) sprintf("%.3f (%.3f, %.3f)", r$rmst_maxtau_ci[1], r$rmst_maxtau_ci[2], r$rmst_maxtau_ci[3])
    row(sprintf("RMST difference and 95%% interval at the last time (%.1f)", max(rk$at_points)), rr(rk), rr(re))
    row("  RMST band critical value", f4(rk$c_alpha_band), f4(re$c_alpha_band))
    # figures: the vignette's weighted band plot and RMST plot under the two options, side by side (common y range for the band plot)
    yl <- range(c(tk$fit_itt$sb_lower, tk$fit_itt$sb_upper, te$fit_itt$sb_lower, te$fit_itt$sb_upper))
    graphics::par(mfrow = c(2, 2), oma = c(0, 0, 2.5, 0))
    for (wv in c("known", "estimated")) {
      a <- c(pk_args, list(ymin = round(yl[1], 2), ymax = round(yl[2], 2) + 0.02), if (wv == "estimated") list(weight_variance = "estimated", weight_model = pw))
      do.call(plotKM.band_subgroups, a)
      graphics::legend("topleft", c("Difference estimate", "95% pointwise CIs"), lwd = c(2, 1), col = "black", lty = c(1, 2), bty = "n", cex = 0.7)
      graphics::title(main = sprintf("plotKM.band_subgroups(), weight_variance = \"%s\"\nband critical value %.3f", wv, (if (wv == "known") tk else te)$fit_itt$c_alpha_band), cex.main = 0.95)
    }
    for (wv in c("known", "estimated")) {
      cumulative_rmst_bands(df = d, fit = (if (wv == "known") tk else te)$fit_itt, tte.name = "tte", event.name = "event", treat.name = "treat", weight.name = "sw.weights",
                            draws_sb = 1000, xlab = "months", rmst_max_cex = 0.7, weight_variance = wv, weight_model = if (wv == "estimated") pw else NULL)
      graphics::title(main = sprintf("cumulative_rmst_bands(), weight_variance = \"%s\"", wv), cex.main = 0.95)
    }
    graphics::mtext(sprintf("%s (%s:%d): weights treated as known (left) and as estimated (right)", ex$id, ex$file, ex$line), outer = TRUE, cex = 1.05, font = 2)
  } else {
    P01 <- list(rho = 0, gamma = 1); dk1 <- dfc("known", P01); de1 <- dfc("estimated", P01)
    c1k <- q(cox_rhogamma(dk1, scheme = "fh", scheme_params = P01, draws = 0)); c1e <- q(cox_rhogamma(de1, scheme = "fh", scheme_params = P01, draws = 0))
    row("FH(0,1): z; hazard ratio and 95% interval", sprintf("%s; %s", f4(c1k$z.score), ci(c1k)), sprintf("%s; %s", f4(c1e$z.score), ci(c1e)))
    cat("  (this example has no band or RMST call)\n")
  }
}
invisible(grDevices::dev.off())
sink()
