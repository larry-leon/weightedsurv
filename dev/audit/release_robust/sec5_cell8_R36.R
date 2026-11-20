## Section 5: OC cell 8 (E-400), RMST, month 36, recomputed from the stored per-replicate files
## against the fixed truth R(36) (the CSV's true_value), as M7 does. No rerun. Run from the repo root.
## Prints the stored CSV values, the recomputed ones, and the M7 cell 31 values for comparison.
csv <- "vignettes/articles/sim_results/ps_weight_results_by_time.csv"
X <- read.csv(csv)
row <- which(X$source == "OC" & X$cell_id == 8 & X$estimand == "R" & X$time %in% 36)
stopifnot(length(row) == 1)
tru <- X$true_value[row]
P <- readRDS("dev/audit/oc/raw/proposed_cell08.rds")
M <- readRDS("dev/audit/stage2e/raw/m7_cell31.rds")
cat("same data (chk identical):", identical(unname(P$chk), unname(M$chk)),
    "; max |est cell 8 - est cell 31|:", max(abs(P$res[, "proposed.R.est.36"] - M$res[, "known.R.est.36"])), "\n")
P <- P$res
ev <- P[, "evaluable.36"] == 1
cov <- function(e, rows = rep(TRUE, nrow(P))) {
  lo <- P[rows, paste0(e, ".R.lo.36")]; hi <- P[rows, paste0(e, ".R.hi.36")]
  c_ <- (is.finite(lo) & is.finite(hi) & lo <= tru & tru <= hi)[ev[rows]]
  c(mean(c_), sqrt(mean(c_) * (1 - mean(c_)) / length(c_)))
}
e <- P[ev & is.finite(P[, "proposed.R.est.36"]), "proposed.R.est.36"]
first <- seq_len(X$boot_replicates[row])
new <- c(mean_estimate = mean(e), bias = mean(e - tru), bias_mcse = sd(e - tru) / sqrt(length(e)),
         coverage_known = cov("proposed")[1], coverage_known_mcse = cov("proposed")[2],
         coverage_estimated = cov("proposed_ew")[1], coverage_estimated_mcse = cov("proposed_ew")[2],
         coverage_known_bootreps = cov("proposed_rect", first)[1],
         coverage_estimated_bootreps = cov("proposed_ew", first)[1],
         empirical_sd = sd(e))
m7 <- X[X$source == "M7" & X$cell_id == 31 & X$estimand == "R" & X$time %in% 36, ]
out <- data.frame(column = names(new), stored = unlist(X[row, names(new)]), recomputed = unname(new),
                  m7_cell31 = if (nrow(m7) == 1) unlist(m7[, names(new)]) else NA, row.names = NULL)
print(format(out, digits = 15), row.names = FALSE)
cat("mean of the per-replicate truths R(t_j), t_j <= 36:", format(mean(P[ev, "R.truth.36"]), digits = 15),
    " (stored bias = mean(est - R(t_j)):", format(mean(P[ev, "proposed.R.est.36"] - P[ev, "R.truth.36"]), digits = 15), ")\n")
saveRDS(new, commandArgs(TRUE)[1])
