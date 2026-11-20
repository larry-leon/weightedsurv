# oc_g1.R -- Gate G1 item 3 (harness, 200 replicates of cells 1 and 5 under both versions) and the
# 200-replicate part of Part D gate 5a (cell 8). Reads dev/audit/oc/raw_g1/.
source("dev/audit/oc/oc_common.R")
rd <- function(v, k) readRDS(sprintf("dev/audit/oc/raw_g1/%s_cell%02d.rds", v, k))
fails <- character(0)
for (k in c(1, 5)) {
  a <- rd("cran", k); b <- rd("proposed", k)
  same <- identical(unname(a$chk), unname(b$chk))
  cat(sprintf("G1.3 cell %d: checksums identical between versions in all %d replicates: %s\n", k, length(a$chk), same))
  if (!same) fails <- c(fails, "checksums")
}
a <- rd("cran", 1)$res; b <- rd("proposed", 1)$res
nf <- c(cran = sum(a[, grep("\\.fail$", colnames(a))] == 1, na.rm = TRUE), proposed = sum(b[, grep("\\.fail$", colnames(b))] == 1, na.rm = TRUE))
cat(sprintf("G1.3 cell 1: calls that errored: cran %d, proposed %d\n", nf[1], nf[2]))
if (any(nf > 0)) fails <- c(fails, "cell 1 errors")
tr <- delta_true(18, 1)
bias <- mean(b[, "proposed.d.est.18"]) - tr
cov <- mean(b[, "proposed.d.lo.18"] <= tr & tr <= b[, "proposed.d.hi.18"])
cat(sprintf("G1.3 cell 1, proposed, delta(18): truth %.5f bias %.5f (limit |bias| 0.01); pointwise coverage %.3f (limits 0.91-0.99)\n", tr, bias, cov))
if (abs(bias) > 0.01) fails <- c(fails, "bias")
if (cov < 0.91 || cov > 0.99) fails <- c(fails, "coverage")
e <- rd("proposed", 8)$res
ratio <- mean(e[, "proposed_ew.d.se.18"]) / mean(e[, "proposed.d.se.18"])
sesd <- mean(e[, "proposed_ew.d.se.18"]) / sd(e[, "proposed.d.est.18"])
cat(sprintf("G1.5a (200 replicates of cell 8), delta(18): mean proposed_ew SE / mean proposed SE = %.3f (limits 0.75-0.98); proposed_ew SE/SD = %.3f (limits 0.80-1.15); proposed SE/SD = %.3f\n",
            ratio, sesd, mean(e[, "proposed.d.se.18"]) / sd(e[, "proposed.d.est.18"])))
g5a <- ratio >= 0.75 && ratio <= 0.98 && sesd >= 0.80 && sesd <= 1.15
cat("Part D gate 5a (200 replicates):", if (g5a) "PASS" else "FAIL", "\n")
if (length(fails)) { cat("G1.3 FAIL:", paste(fails, collapse = ", "), "\n"); quit(status = 1) }
cat("G1 item 3: PASS\n")
