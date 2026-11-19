#!/usr/bin/env Rscript
# =============================================================================
# report_numbers.R -- tables B and C of REPORT_km_fixes_stage2b.md
#   B: weighted SE (KM_diff, weights wv, draws = 0) at the 25th / 50th / 75th percentile of the time
#      grid, before (golden_s2.rds) -> after (current code), and the ratio to survfit(robust = TRUE).
#   C: B8 -- RMST at the last grid time under "trapezoid" and "rectangle", survfit's restricted-mean
#      difference, and the band critical value under each rule (KM_diff(draws.band = 1000),
#      cumulative_rmst_bands() defaults otherwise).
# Writes dev/audit/stage2b/report_numbers.txt. Datasets are those of golden.R (same code, same seeds).
# =============================================================================
if (!file.exists("DESCRIPTION")) stop("Run from the repo root.")
suppressMessages(pkgload::load_all(".", quiet = TRUE, export_all = FALSE, helpers = FALSE))
S2 <- file.path("dev", "audit", "stage2b")
# dataset code from golden.R (function definitions only)
ex <- parse(file.path(S2, "golden.R"))
for (e in ex) if (is.call(e) && identical(e[[1]], as.name("<-")) && as.character(e[[2]]) %in% c("sim", "add_w", "add_wv", "make_data")) eval(e)
D <- make_data()
gold <- readRDS(file.path(S2, "golden_s2.rds"))$values
out <- character(0); say <- function(...) out <<- c(out, sprintf(...))

say("B. Weighted SE, KM_diff(weight.name = \"wv\", draws = 0): before (golden_s2) -> after; ratio to survfit(weights, robust = TRUE) before -> after")
say("| dataset | arm | t (grid pct) | SE before | SE after | ratio before | ratio after |")
say("|---|---|---|---|---|---|---|")
for (nm in names(D)) {
  d <- D[[nm]]
  k <- suppressWarnings(KM_diff(d, "tte", "event", "treat", weight.name = "wv"))
  ap <- k$at_points
  idx <- ceiling(c(0.25, 0.5, 0.75) * length(ap))
  for (a in 0:1) {
    before <- sqrt(gold[[sprintf("%s/KM_diff_wv/sig2_surv%d", nm, a)]])[idx]
    after <- sqrt(k[[sprintf("sig2_surv%d", a)]])[idx]
    x <- d[d$treat == a, ]
    fr <- survival::survfit(survival::Surv(tte, event) ~ 1, data = x, weights = wv, robust = TRUE)
    rob <- summary(fr, times = ap[idx], extend = TRUE)$std.err
    for (j in 1:3) say("| %s | %d | %.3g (%d%%) | %.5f | %.5f | %.3f | %.3f |", nm, a, ap[idx[j]], c(25, 50, 75)[j],
                       before[j], after[j], before[j] / rob[j], after[j] / rob[j])
  }
}

say("")
say("C. B8: RMST difference at the last grid time; band critical value under each rule")
say("| dataset | RMST trapezoid | RMST rectangle | survfit rmean diff | rectangle - survfit | c_alpha trapezoid | c_alpha rectangle |")
say("|---|---|---|---|---|---|---|")
v <- survival::veteran; v <- data.frame(tte = v$time, event = v$status, treat = v$trt - 1)
g <- survival::gbsg; g <- data.frame(tte = g$rfstime / 30.4375, event = g$status, treat = g$hormon)
C3 <- list(veteran = v, gbsg = g, sim_cont = D$sim_cont[, c("tte", "event", "treat")])
for (nm in names(C3)) {
  d <- C3[[nm]]
  k <- KM_diff(d, "tte", "event", "treat", draws.band = 1000)
  rt <- cumulative_rmst_bands(d, k, "tte", "event", "treat", plot = FALSE)
  rr <- cumulative_rmst_bands(d, k, "tte", "event", "treat", plot = FALSE, mean_rule = "rectangle")
  tau <- max(k$at_points)
  rm <- summary(survival::survfit(survival::Surv(tte, event) ~ treat, data = d), rmean = tau)$table[, "rmean"]
  sf <- unname(rm[2] - rm[1])
  n <- length(rt$rmst_time)
  say("| %s (tau = %.4g) | %.6f | %.6f | %.6f | %.2g | %.4f | %.4f |", nm, tau, rt$rmst_time[n], rr$rmst_time[n], sf,
      rr$rmst_time[n] - sf, rt$c_alpha_band, rr$c_alpha_band)
}
writeLines(out, file.path(S2, "report_numbers.txt"))
cat(out, sep = "\n")
