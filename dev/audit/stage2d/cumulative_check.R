#!/usr/bin/env Rscript
# cumulative_check.R -- the cumulative variance V(t) of wlr_cumulative(<df_counting result>, return_cumulative = TRUE),
# FH(0,0), against coxph(robust = TRUE, ties = "breslow") on the same data censored at t (time = pmin(time, t),
# status = status * (time <= t)): z(t)^2 should equal that fit's robust score test. Writes cumulative_check.txt.
if (!file.exists("DESCRIPTION")) stop("Run from the repo root.")
suppressMessages(pkgload::load_all(".", quiet = TRUE)); suppressPackageStartupMessages(library(survival))
source("tests/testthat/helper-score-resid.R")
out <- character(0)
for (ties in c(FALSE, TRUE)) {
  d <- ws_sim_ipw(300, 4, ties)
  dc <- ws_quiet(df_counting(d, "time", "status", "treat", weight.name = "w"))
  cum <- wlr_cumulative(dc, "fh", list(rho = 0, gamma = 0), return_cumulative = TRUE)
  tt <- as.numeric(stats::quantile(dc$time[dc$delta == 1], c(0.1, 0.25, 0.5, 0.75, 0.9)))
  for (t in tt) {
    dd <- d; dd$status <- dd$status * (dd$time <= t); dd$time <- pmin(dd$time, t)
    cp <- coxph(Surv(time, status) ~ treat, data = dd, weights = w, robust = TRUE, ties = "breslow")
    j <- max(which(cum$time <= t))
    out <- c(out, sprintf("ties = %-5s t = %6.2f: events 0/1 up to t %d/%d; cumulative z^2 %.8f, coxph robust score on data censored at t %.8f, relative difference %.2g",
                          ties, t, sum(dd$status == 1 & dd$treat == 0), sum(dd$status == 1 & dd$treat == 1), cum$z.score[j]^2, cp$rscore, abs(cum$z.score[j]^2 / cp$rscore - 1)))
  }
}
writeLines(out, "dev/audit/stage2d/cumulative_check.txt"); cat(out, sep = "\n")
