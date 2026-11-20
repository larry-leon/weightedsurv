# oc_g0.R -- Gate G0 checks: both versions load; two identity checks on one cell-5 dataset under HEAD.
# Usage: Rscript dev/audit/oc/oc_g0.R <path to 0.1.0 worktree>
suppressPackageStartupMessages(library(survival))
source("dev/audit/oc/oc_common.R")
cran_path <- commandArgs(TRUE)[1]
ok <- tryCatch({ suppressMessages(pkgload::load_all(cran_path, quiet = TRUE)); TRUE }, error = function(e) { message(e); FALSE })
cat("load_all 0.1.0:", ok, " version:", as.character(utils::packageVersion("weightedsurv")), "\n")
pkgload::unload("weightedsurv")
ok2 <- tryCatch({ suppressMessages(pkgload::load_all(".", quiet = TRUE)); TRUE }, error = function(e) { message(e); FALSE })
cat("load_all HEAD :", ok2, " version:", as.character(utils::packageVersion("weightedsurv")), "\n")
stopifnot(ok, ok2)
d <- make_data(5, 1); df <- d$df
fit <- KM_diff(df, tte.name = "time", event.name = "status", treat.name = "treat", weight.name = "w",
               alpha = 0.05, draws = 0, draws.band = 2000, seedstart = 8316951 + 1, show_resamples = FALSE)
j <- last_le(fit$at_points, 18)
sf <- survfit(Surv(time, status) ~ 1, data = df[df$treat == 1, ], weights = w, robust = TRUE)
se_sf <- summary(sf, times = 18)$std.err
se_pk <- sqrt(fit$sig2_surv1[j])
rel <- abs(se_pk - se_sf) / se_sf
cat(sprintf("identity 1: SE surv1(18) package %.10f survfit(robust) %.10f relative difference %.3g\n", se_pk, se_sf, rel))
sel <- fit$surv1 > 0.1 & fit$surv1 < 0.9
ratio <- apply(fit$surv1_star, 1, var)[sel] / fit$sig2_surv1[sel]
cat(sprintf("identity 2: var over 2000 draws of surv1_star / sig2_surv1, where 0.1 < surv1 < 0.9: min %.4f max %.4f (%d grid points)\n",
            min(ratio), max(ratio), sum(sel)))
if (rel > 1e-6 || min(ratio) < 0.85 || max(ratio) > 1.15) { cat("G0 FAIL\n"); quit(status = 1) }
cat("G0 identity checks PASS\n")
