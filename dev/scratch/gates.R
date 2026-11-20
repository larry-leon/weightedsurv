#!/usr/bin/env Rscript
# gates.R -- Gates G7, G8, G9 (first part) of TASK_vignette_revision.md, from the objects of the knitted vignettes
# (dev/scratch/vignette_after/<file>/env.rds, written by knit_vignettes.R) and of the knit at ac76f82
# (dev/scratch/vignette_before/).   Rscript dev/scratch/gates.R <lib>
args <- commandArgs(trailingOnly = TRUE)
.libPaths(c(normalizePath(args[1]), .libPaths()))
suppressPackageStartupMessages({library(survival); library(weightedsurv)})
inr <- function(x, lo, hi) is.finite(x) && x >= lo && x <= hi
say <- function(gate, ok, txt) cat(sprintf("%s %s: %s\n", if (ok) "PASS" else "FAIL", gate, txt))
fh00 <- function(dc) cox_rhogamma(dc, scheme = "fh", scheme_params = list(rho = 0, gamma = 0))
all_ok <- TRUE
for (nm in c("examples", "martingale", "coxmeth")) {
  e <- readRDS(file.path("dev/scratch/vignette_after", nm, "env.rds"))
  d <- e$df_rotterdam
  dc <- if (nm == "coxmeth") e$dfcount_wtd else e$dfcount_rotterdam_wtd
  known <- df_counting(df = d, tte.name = "tte", event.name = "event", treat.name = "treat", arms = c("treat", "control"), weight.name = "sw.weights")
  v <- c(n = nrow(d), n0 = sum(d$treat == 0), n1 = sum(d$treat == 1), ess_ratio_treated = dc$ess$ratio[1], max_w = max(d$sw.weights),
         hr = fh00(dc)$hr_ci_asy$hr, z_known = known$z.score, z_estimated = dc$z.score)
  ok <- v[["n"]] == 829 && v[["n0"]] == 538 && v[["n1"]] == 291 && inr(v[["ess_ratio_treated"]], 0.92, 0.94) && inr(v[["max_w"]], 2.6, 2.8) &&
    inr(v[["hr"]], 0.655, 0.675) && inr(v[["z_known"]], 4.40, 4.50) && inr(v[["z_estimated"]], 4.82, 4.92) &&
    identical(dc$weight_variance, "estimated") && isTRUE(all.equal(d$sw.weights, e$pw$weights, tolerance = 0))
  all_ok <- all_ok && ok
  say(paste0("G7 [", nm, "]"), ok, sprintf("n = %d; %d / %d; treated ESS ratio %.4f; largest weight %.4f; weighted HR %.4f; z known %.4f, estimated %.4f; weight_variance of the vignette object \"%s\"",
                                           v[["n"]], v[["n0"]], v[["n1"]], v[["ess_ratio_treated"]], v[["max_w"]], v[["hr"]], v[["z_known"]], v[["z_estimated"]], dc$weight_variance))
}
e <- readRDS("dev/scratch/vignette_after/examples/env.rds")
# the vignette's own known-vs-estimated table and its known object
cat("   vignette table (examples):\n"); print(e$known_vs_estimated, row.names = FALSE)
cat(sprintf("   z of the vignette's known object: %.4f\n", e$dfcount_rotterdam_known$z.score))
ctc <- e$df_CTC; dcc <- e$dfcount_CTC_wtd; hc <- fh00(dcc)
ok <- sum(ctc$treat == 0) == 1207 && sum(ctc$treat == 1) == 246 && inr(dcc$ess$ratio[1], 0.29, 0.30) && inr(hc$hr_ci_asy$hr, 0.655, 0.665) &&
  identical(dcc$weight_variance, "known") && identical(hc$weight_variance, "known") && identical(e$band_CTC$fit_itt$weight_variance, "known") &&
  identical(e$rmst_CTC$weight_variance, "known") && identical(e$dfcount_CTC_cut$weight_variance, "known")
all_ok <- all_ok && ok
say("G8", ok, sprintf("%d / %d; treated ESS ratio %.4f; HR known %.4f (%.4f, %.4f); weight_variance: df_counting \"%s\", cox_rhogamma \"%s\", band \"%s\", RMST \"%s\", cut at 60 \"%s\"",
                      sum(ctc$treat == 0), sum(ctc$treat == 1), dcc$ess$ratio[1], hc$hr_ci_asy$hr, hc$hr_ci_asy$lower, hc$hr_ci_asy$upper,
                      dcc$weight_variance, hc$weight_variance, e$band_CTC$fit_itt$weight_variance, e$rmst_CTC$weight_variance, e$dfcount_CTC_cut$weight_variance))
h6 <- fh00(e$dfcount_CTC_cut)
cat(sprintf("   cross-trial: ESS treated %.1f (%.3f), control %.1f (%.3f); largest weight %.2f; z %.3f; cut at 60: HR %.2f (%.2f, %.2f), z %.2f; 10th treated event %.2f months\n",
            dcc$ess$ess[1], dcc$ess$ratio[1], dcc$ess$ess[2], dcc$ess$ratio[2], max(ctc$sw.weights), dcc$z.score, h6$hr_ci_asy$hr, h6$hr_ci_asy$lower, h6$hr_ci_asy$upper,
            e$dfcount_CTC_cut$z.score, sort(ctc$tte[ctc$treat == 1 & ctc$event == 1])[10]))
np <- e$dfcount_np_wtd$ess
cat(sprintf("   all node-positive: ESS treated %d -> %.1f (%.3f), control %d -> %.1f (%.3f); largest weight %.2f\n", np$n[1], np$ess[1], np$ratio[1], np$n[2], np$ess[2], np$ratio[2],
            max(e$df_rotterdam_np$sw.weights)))
ck <- fh00(e$dfcount_rotterdam_known)$hr_ci_asy; ce <- fh00(e$dfcount_rotterdam_wtd)$hr_ci_asy; cu <- fh00(e$dfcount_rotterdam_unwtd)$hr_ci_asy
cat(sprintf("   restricted: HR known %.2f (%.2f, %.2f); estimated (%.2f, %.2f); unweighted %.2f (%.2f, %.2f); control ESS ratio %.2f; worst imbalance %.2f -> %.2f\n",
            ck$hr, ck$lower, ck$upper, ce$lower, ce$upper, cu$hr, cu$lower, cu$upper, e$dfcount_rotterdam_wtd$ess$ratio[2], max(abs(e$balance$unweighted)), max(abs(e$balance$weighted))))
# G9 (first part): unweighted df_counting() on the unrestricted node-positive data, before (object of the knit at ac76f82) and after
b <- readRDS("dev/scratch/vignette_before/examples/env.rds")
before <- b$dfcount_rotterdam_unwtd
dnp <- e$df_rotterdam_np
after <- get_dfcounting(df = dnp, tte.name = "tte", event.name = "event", treat.name = "treat", arms = c("treat", "control"), by.risk = 24)
strip <- function(x) { x$cox_results$cox_fit <- NULL; x }   # the coxph object holds the call and an environment
# columns of the data as prepared (the weight columns and hormon2, added later in the vignette at ac76f82, are left out)
cols <- setdiff(names(b$df_rotterdam), c("sw.weights", "sw.weights_trunc", "hormon2"))
ok <- identical(b$df_rotterdam[, cols], dnp[, cols]) && identical(strip(before), strip(after)) &&
  isTRUE(all.equal(coef(before$cox_results$cox_fit), coef(after$cox_results$cox_fit), tolerance = 0)) &&
  isTRUE(all.equal(before$cox_results$cox_fit$var, after$cox_results$cox_fit$var, tolerance = 0))
all_ok <- all_ok && ok
say("G9a", ok, sprintf("node-positive data identical() before/after (%d rows, %d columns compared); unweighted df_counting() result identical() in all %d elements (cox_fit compared by coefficients and variance, tolerance 0); z = %.6f",
                       nrow(dnp), length(cols), length(names(after)), after$z.score))
# also: GBSG unweighted object of the vignette, before against after
ok2 <- identical(strip(b$dfcount_gbsg), strip(e$dfcount_gbsg)); all_ok <- all_ok && ok2
say("G9a (GBSG object of the vignette)", ok2, "dfcount_gbsg identical() before/after")
quit(status = if (all_ok) 0 else 1)
