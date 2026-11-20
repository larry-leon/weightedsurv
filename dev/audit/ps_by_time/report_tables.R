#!/usr/bin/env Rscript
# report_tables.R -- the comparative table of the new cells for REPORT_ps_weight_results_by_time.md (section 3),
# from vignettes/articles/sim_results/ps_weight_results_by_time.csv.
#   Rscript --vanilla dev/audit/ps_by_time/report_tables.R > dev/audit/ps_by_time/report_tables.md
R <- read.csv("vignettes/articles/sim_results/ps_weight_results_by_time.csv", stringsAsFactors = FALSE)
f <- function(x, d) ifelse(is.na(x), "", formatC(x, format = "f", digits = d))
cv <- function(p, s) ifelse(is.na(p), "", sprintf("%s (%s)", f(p, 3), f(s, 3)))
row <- function(x, lab, d) c(x$cell, lab, f(x$true_value, d), f(x$bias, d + 1), f(x$empirical_sd, d + 1), f(x$mean_se_known, d + 1),
                             f(x$mean_se_estimated, d + 1), cv(x$coverage_known, x$coverage_known_mcse),
                             cv(x$coverage_estimated, x$coverage_estimated_mcse), f(x$mean_ci_length_known, d), f(x$mean_ci_length_estimated, d))
out <- NULL
for (cl in c("Es9-400", "Es9-400-null")) {
  for (q in c("delta", "R")) for (t in c(1, 6, 12, 24, 36, 48)) {
    x <- R[R$cell == cl & R$estimand == q & R$time %in% t, ]
    out <- rbind(out, row(x, sprintf("%s, month %d", if (q == "delta") "difference" else "RMST", t), 3))
  }
  for (q in c("loghr_fh00", "loghr_fh01")) {
    x <- R[R$cell == cl & R$estimand == q, ]
    out <- rbind(out, row(x, sprintf("log HR %s", if (q == "loghr_fh00") "FH(0,0)" else "FH(0,1)"), 3))
  }
  for (q in c("delta_band_3_42", "R_band_3_42")) {
    x <- R[R$cell == cl & R$estimand == q, ]
    out <- rbind(out, c(cl, sprintf("band [3, 42], %s (width)", if (q == "delta_band_3_42") "difference" else "RMST"), "", "", "", "", "",
                        cv(x$coverage_known, x$coverage_known_mcse), cv(x$coverage_estimated, x$coverage_estimated_mcse),
                        f(x$mean_ci_length_known, 3), f(x$mean_ci_length_estimated, 3)))
  }
}
h <- c("Cell", "Estimand", "True", "Bias", "Emp SD", "SE known", "SE est", "Cov known (MCSE)", "Cov est (MCSE)", "Length known", "Length est")
cat(paste0("| ", paste(h, collapse = " | "), " |"), paste0("|", paste(rep("---", length(h)), collapse = "|"), "|"),
    apply(out, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")), sep = "\n")
