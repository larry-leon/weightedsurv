#!/usr/bin/env Rscript
# =============================================================================
# report_numbers.R -- table B of REPORT_km_fixes_stage2c.md (R1 and R2 in numbers)
#   Per golden dataset (weights wv): weighted SE of KM_diff(draws = 0) at the 25th / 50th / 75th percentile
#   of the time grid, Stage 2b (golden_s2b.rds) -> Stage 2c (current code), with the ratio to
#   survfit(weights = wv, robust = TRUE) on the time-fixed data; band critical value of
#   KM_diff(draws.band = 200) before -> after; RMST pointwise limits and band limits at the last grid time
#   (cumulative_rmst_bands on that fit) before -> after. Same datasets and calls as golden.R.
# Writes dev/audit/stage2c/report_numbers.txt.
# =============================================================================
if (!file.exists("DESCRIPTION")) stop("Run from the repo root.")
suppressMessages(pkgload::load_all(".", quiet = TRUE, export_all = FALSE, helpers = FALSE))
S2 <- file.path("dev", "audit", "stage2c")
ex <- parse(file.path(S2, "golden.R"))
for (e in ex) if (is.call(e) && identical(e[[1]], as.name("<-")) && as.character(e[[2]]) %in% c("sim", "add_w", "add_wv", "make_data")) eval(e)
D <- make_data()
gold <- readRDS(file.path(S2, "golden_s2b.rds"))$values
out <- character(0); say <- function(...) out <<- c(out, sprintf(...))
say("B1. Weighted SE, KM_diff(weight.name = \"wv\", draws = 0): Stage 2b -> Stage 2c; ratio to survfit(weights, robust = TRUE)")
say("| dataset | arm | t (grid pct) | SE 2b | SE 2c | ratio 2b | ratio 2c |")
say("|---|---|---|---|---|---|---|")
B2 <- character(0)
for (nm in names(D)) {
  d <- D[[nm]]
  k <- KM_diff(d, "tte", "event", "treat", weight.name = "wv")
  ap <- k$at_points; idx <- ceiling(c(0.25, 0.5, 0.75) * length(ap))
  for (a in 0:1) {
    before <- sqrt(gold[[sprintf("%s/KM_diff_wv/sig2_surv%d", nm, a)]])[idx]
    after <- sqrt(k[[sprintf("sig2_surv%d", a)]])[idx]
    x <- d[d$treat == a, ]
    y <- survival::aeqSurv(survival::Surv(x$tte, x$event))
    rob <- summary(survival::survfit(y ~ 1, weights = x$wv, robust = TRUE), times = ap[idx], extend = TRUE)$std.err
    for (j in 1:3) say("| %s | %d | %.3g (%d%%) | %.5f | %.5f | %.3f | %.4f |", nm, a, ap[idx[j]], c(25, 50, 75)[j],
                       before[j], after[j], before[j] / rob[j], after[j] / rob[j])
  }
  kb <- KM_diff(d, "tte", "event", "treat", weight.name = "wv", draws.band = 200, show_resamples = FALSE)
  rb <- cumulative_rmst_bands(d, kb, "tte", "event", "treat", weight.name = "wv", plot = FALSE)
  n <- length(rb$rmst_time)
  g <- function(f) gold[[sprintf("%s/%s", nm, f)]]
  B2 <- c(B2, sprintf("| %s | %.4f -> %.4f | %.3f | %.3f, %.3f -> %.3f, %.3f | %.3f -> %.3f | %.3f, %.3f -> %.3f, %.3f |",
    nm, g("KM_diff_wv_band/c_alpha_band"), kb$c_alpha_band, rb$rmst_time[n],
    g("rmst_band_wv/rmst_time_lower")[n], g("rmst_band_wv/rmst_time_upper")[n], rb$rmst_time_lower[n], rb$rmst_time_upper[n],
    g("rmst_band_wv/c_alpha_band"), rb$c_alpha_band,
    g("rmst_band_wv/rmst_time_sb_lower")[n], g("rmst_band_wv/rmst_time_sb_upper")[n], rb$rmst_time_sb_lower[n], rb$rmst_time_sb_upper[n]))
}
say("")
say("B2. Bands: KM_diff(weights wv, draws.band = 200) critical value; cumulative_rmst_bands at the last grid time")
say("| dataset | KM band c_alpha 2b -> 2c | RMST (unchanged) | RMST pointwise limits 2b -> 2c | RMST band c_alpha 2b -> 2c | RMST band limits 2b -> 2c |")
say("|---|---|---|---|---|---|")
out <- c(out, B2)
writeLines(out, file.path(S2, "report_numbers.txt")); cat(out, sep = "\n")
