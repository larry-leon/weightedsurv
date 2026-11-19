#!/usr/bin/env Rscript
# km_m2b_check.R -- M2b requirements (TASK_km_fixes_stage2c.md, Step 4) and Table C for the report.
#   (1) pkg equals ij: max relative difference < 1e-8 over replicates in which survfit merged no times;
#   (2) first 500 replicates of each arm-cell: mean(pkg_resampling SE) / mean(pkg SE) in [0.97, 1.03] at t = 5, 10, 20.
# Prints Table C (t = 10, SE/SD and coverage) beside the M2 (Stage 2b) values. Exit 0 when (1) and (2) hold.
rep <- readRDS("dev/audit/stage2c/km_m2_replicates.rds")$replicates
res <- read.csv("dev/audit/stage2c/km_m2_results.csv"); old <- read.csv("dev/audit/m2/km_m2_results.csv")
nm <- rep$n_merged == 0 & rep$ij > 0
r1 <- max(abs(rep$pkg[nm] / rep$ij[nm] - 1))
cat(sprintf("(1) max |pkg / ij - 1| over %d arm-cell x replicate x time values without merged times: %.3g\n", sum(nm), r1))
r1all <- max(abs(rep$pkg / rep$robust - 1)[rep$robust > 0 & rep$n_merged == 0])
cat(sprintf("    max |pkg / robust - 1| there: %.3g\n", r1all))
f <- rep[rep$r <= 500, ]
rt <- aggregate(cbind(pkg_resampling, pkg) ~ armcell + t, data = f, FUN = mean)
rt$ratio <- rt$pkg_resampling / rt$pkg
rt <- rt[order(match(rt$armcell, unique(rep$armcell)), rt$t), ]
cat("(2) mean(pkg_resampling) / mean(pkg), first 500 replicates:\n")
print(reshape(rt[, c("armcell", "t", "ratio")], idvar = "armcell", timevar = "t", direction = "wide"), row.names = FALSE, digits = 4)
ok2 <- all(rt$ratio >= 0.97 & rt$ratio <= 1.03)
cat(sprintf("    range %.4f to %.4f\n", min(rt$ratio), max(rt$ratio)))
cat("\nTable C, t = 5, 10, 20: SE/SD (coverage); M2 = Stage-2b code, M2b = Stage-2c code\n")
g <- function(d, e, t) { x <- d[d$estimator == e & d$t == t, ]; setNames(sprintf("%.3f (%.3f)", x$se_over_sd, x$coverage), x$armcell) }
for (t in c(5, 10, 20)) {
  cat(sprintf("\n t = %d\n| arm-cell | pkg M2 | pkg_resampling M2 | pkg M2b | pkg_resampling M2b | robust |\n|---|---|---|---|---|---|\n", t))
  a <- g(old, "pkg", t); b <- g(old, "pkg_resampling", t); c1 <- g(res, "pkg", t); d1 <- g(res, "pkg_resampling", t); r <- g(res, "robust", t)
  for (ac in names(a)) cat(sprintf("| %s | %s | %s | %s | %s | %s |\n", ac, a[ac], b[ac], c1[ac], d1[ac], r[ac]))
}
ok <- r1 < 1e-8 && ok2
cat(if (ok) "\nM2b: PASS\n" else "\nM2b: FAIL\n"); quit(status = if (ok) 0 else 1)
