#!/usr/bin/env Rscript
# ps9_summarise.R -- TASK_ps_weight_results_by_time.md, sections 2 and 3.
#
# 1. Summarises the gap run (raw/m7_cell41.rds, raw/m7_cell42.rds) with the rules of dev/audit/stage2e/oc_summarise.R:
#    the definitions and the per-cell loop of that script are evaluated here on this folder's raw/ (its file-writing
#    code is not evaluated, and the script is not edited). Output: ps9_summary.csv (long, the columns of m7_summary.csv).
# 2. Builds vignettes/articles/sim_results/ps_weight_results_by_time.csv (one row per cell x estimand x time point,
#    the columns of the inventory CSV) from ps9_summary.csv, dev/audit/stage2e/m7_summary.csv and
#    dev/audit/oc/oc_summary.csv; the ESS ratios of the OC reference cells are read from dev/audit/oc/raw/.
# 3. Gate checks G2, G4, G5 (and requirement 3 of M7 for the new cells), written to ps9_gates.txt.
#
#   Rscript --vanilla dev/audit/ps_by_time/ps9_summarise.R

.here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
S2E <- normalizePath(file.path(.here, "..", "stage2e")); OCD <- normalizePath(file.path(.here, "..", "oc"))
REPO <- normalizePath(file.path(.here, "..", "..", ".."))
source(file.path(S2E, "oc_common.R"))
PA <- -0.44; PB <- 0.51
new <- data.frame(k = c(41, 42), cell = c("Es9-400", "Es9-400-null"), sc = c("ph", "null"), type = "stabilized",
                  reps = 2000L, repro = FALSE, strat = FALSE, sub = FALSE, band = TRUE, rmst = TRUE,
                  fh01 = TRUE, both = FALSE, coxph = FALSE, stringsAsFactors = FALSE)
new$loghr_truth <- ifelse(new$sc == "ph", -0.4, 0)
m7 <- rbind(m7, new[, names(m7)])

# ---- 1. the summary rules of stage2e/oc_summarise.R ----
raw <- file.path(.here, "raw"); Q <- qnorm(0.975); NT <- length(tt)
ex <- as.list(parse(file.path(S2E, "oc_summarise.R"), keep.source = FALSE))
lhs <- function(e) if (is.call(e) && identical(e[[1]], as.name("<-"))) deparse(e[[2]]) else ""
i0 <- which(vapply(ex, lhs, "") == "rd"); i1 <- which(vapply(ex, lhs, "") == "S") - 1
stopifnot(length(i0) == 1, length(i1) == 1, i1 > i0)
for (e in ex[i0:i1]) eval(e, envir = globalenv())       # rd ... the per-cell loop and the bootstrap block (absent here)
stopifnot(identical(sort(kdone), c(41, 42)))
S <- do.call(rbind, rows)
write.csv(S, file.path(.here, "ps9_summary.csv"), row.names = FALSE, na = "NA")
cat(sprintf("ps9_summary.csv: %d rows (cells %s)\n", nrow(S), paste(kdone, collapse = ", ")))

# ---- 2. the article's CSV ----
rdcsv <- function(f) read.csv(f, stringsAsFactors = FALSE)
L7 <- rdcsv(file.path(S2E, "m7_summary.csv")); LO <- rdcsv(file.path(OCD, "oc_summary.csv"))
LN <- rdcsv(file.path(.here, "ps9_summary.csv"))
oc_ess <- function(k) { r <- readRDS(file.path(OCD, "raw", sprintf("proposed_cell%02d.rds", k))); unname(colMeans(r$res[, c("diag.ess0", "diag.ess1")], na.rm = TRUE)) }

cells_out <- data.frame(
  source  = c("sim", "sim", "M7", "M7", "M7", "M7", "M7", "M7", "OC", "OC"),
  cell_id = c(41, 42, 31, 32, 35, 37, 39, 38, 5, 10),
  regime  = c(rep("Near 0.9", 2), rep("About 0.69", 2), rep("Strong confounding, 0.33-0.39", 2), rep("14% treated, 0.35", 2),
              rep("Reference: true weights", 2)),
  N = c(400, 400, 400, 400, 400, 400, 1450, 1450, 400, 400),
  sc = c("ph", "null", "ph", "null", "ph", "null", "ph", "null", "ph", "ph"),
  propensity_model = c(rep(sprintf("logit e = %.2f + %.2fX", PA, PB), 2), rep("logit e = -0.5 + X", 2), rep("logit e = -0.5 + 2X", 2),
                       rep("logit e = -2.3 + 1.3X", 2), "logit e = -0.5 + X", "logit e = -0.5 + 2X"),
  weights = c(rep("estimated, stabilized", 4), rep("estimated, ipw", 2), rep("estimated, stabilized", 2), rep("true, ipw", 2)),
  stringsAsFactors = FALSE)
cells_out$cell <- c("Es9-400", "Es9-400-null", "Es-400", "Es-400-null", "Fe-400", "Fe-400-null", "V-1450-ph", "V-1450-null", "D-400", "F-400")
scen <- c(ph = "PH, HR = exp(-0.4) = 0.670", null = "null, HR = 1")

COLS <- c("source", "cell_id", "cell", "regime", "N", "scenario", "censoring", "propensity_model", "weights", "ess_ratio_control",
          "ess_ratio_treated", "estimand", "time", "true_value", "mean_estimate", "bias", "bias_mcse", "empirical_sd", "mean_se_known",
          "mean_se_estimated", "coverage_known", "coverage_known_mcse", "coverage_estimated", "coverage_estimated_mcse",
          "mean_ci_length_known", "mean_ci_length_estimated", "replicates", "replicates_evaluable", "n_fail_known", "n_fail_estimated",
          "bias_over_estimated_rows_if_different", "boot_replicates", "coverage_known_bootreps", "coverage_estimated_bootreps",
          "mean_ci_length_known_bootreps", "mean_ci_length_estimated_bootreps", "coverage_boot", "coverage_boot_mcse", "mean_ci_length_boot",
          "coverage_boot_pct", "mean_ci_length_boot_pct", "boot_se_over_sd", "note")
mcse <- function(p, n) sqrt(p * (1 - p) / n)
pick <- function(L, k, q, e, t = NA) {
  i <- L$cell_id == k & L$quantity == q & L$estimator == e & (if (is.na(t)) is.na(L$time) else (!is.na(L$time) & L$time == t))
  if (sum(i) != 1) stop(sprintf("cell %s %s %s %s: %d rows", k, q, e, t, sum(i)))
  L[i, ]
}
na_row <- function() stats::setNames(as.list(rep(NA, length(COLS))), COLS)

build_cell <- function(co) {
  k <- co$cell_id; src <- co$source
  L <- switch(src, sim = LN, M7 = L7, OC = LO)
  kn <- if (src == "OC") "proposed" else "known"; es <- if (src == "OC") NA else "estimated"
  ess <- if (src == "OC") oc_ess(k) else c(pick(L, k, "ess_ratio_control", "known")$mean_value, pick(L, k, "ess_ratio_treated", "known")$mean_value)
  note <- switch(src, sim = "new run (dev/audit/ps_by_time/ps9_sim.R)", M7 = "copied from dev/audit/stage2e/m7_summary.csv",
                 OC = "copied from dev/audit/oc/oc_summary.csv ('proposed'); true weights, 'known' only; ESS ratios from dev/audit/oc/raw/")
  out <- list()
  one <- function(q, t) {
    r <- na_row(); K <- pick(L, k, q, kn, t); E <- if (!is.na(es)) pick(L, k, q, es, t) else NULL
    r[c("source", "cell_id", "cell", "regime", "N", "scenario", "censoring", "propensity_model", "weights")] <-
      list(src, k, co$cell, co$regime, co$N, scen[[co$sc]], "uniform(0, 72) months", co$propensity_model, co$weights)
    r$ess_ratio_control <- ess[1]; r$ess_ratio_treated <- ess[2]; r$estimand <- q; r$time <- t
    band <- grepl("_band_", q)
    if (!band) {
      r$true_value <- K$truth; r$mean_estimate <- K$truth + K$bias; r$bias <- K$bias; r$bias_mcse <- K$bias_mcse
      r$empirical_sd <- K$emp_sd; r$mean_se_known <- K$mean_se
    }
    r$coverage_known <- K$coverage; r$coverage_known_mcse <- mcse(K$coverage, K$n_eval); r$mean_ci_length_known <- K$mean_length
    r$replicates <- K$n_rep; r$replicates_evaluable <- K$n_eval; r$n_fail_known <- K$n_fail
    if (!is.null(E)) {
      if (!band) {
        r$mean_se_estimated <- E$mean_se
        if (!identical(E$bias, K$bias)) r$bias_over_estimated_rows_if_different <- E$bias
      }
      r$coverage_estimated <- E$coverage; r$coverage_estimated_mcse <- mcse(E$coverage, E$n_eval)
      r$mean_ci_length_estimated <- E$mean_length; r$n_fail_estimated <- E$n_fail
    }
    r$note <- note
    r
  }
  for (q in c("S0", "S1", "delta", "R")) for (t in tt) out[[length(out) + 1]] <- one(q, t)
  for (q in c("loghr_fh00", "loghr_fh01", "delta_band_3_42", "R_band_3_42")) {
    if (q == "loghr_fh01" && k %in% c(38, 39)) next                        # not run in cells 38, 39
    out[[length(out) + 1]] <- one(q, NA)
  }
  out
}
W <- do.call(rbind, lapply(unlist(lapply(split(cells_out, seq_len(nrow(cells_out))), build_cell), recursive = FALSE), as.data.frame))

# bootstrap benchmark (first 500 replicates): cells 38, 39 (M7) and 8 (OC), on the rows delta 18, R 36, loghr_fh00
bootq <- c(delta18_first500 = "delta", R36_first500 = "R", loghr_first500 = "loghr_fh00")
boott <- c(delta18_first500 = 18, R36_first500 = 36, loghr_first500 = NA)
fill_boot <- function(r, L, k, q, kn, es) {
  K <- pick(L, k, q, kn); E <- pick(L, k, q, es); B <- pick(L, k, q, "boot"); P <- pick(L, k, q, "boot_pct")
  r$boot_replicates <- B$n_rep; r$coverage_known_bootreps <- K$coverage; r$coverage_estimated_bootreps <- E$coverage
  r$mean_ci_length_known_bootreps <- K$mean_length; r$mean_ci_length_estimated_bootreps <- E$mean_length
  r$coverage_boot <- B$coverage; r$coverage_boot_mcse <- mcse(B$coverage, B$n_eval); r$mean_ci_length_boot <- B$mean_length
  r$coverage_boot_pct <- P$coverage; r$mean_ci_length_boot_pct <- P$mean_length; r$boot_se_over_sd <- B$se_sd
  r
}
for (k in c(38, 39)) for (q in names(bootq)) {
  i <- which(W$source == "M7" & W$cell_id == k & W$estimand == bootq[[q]] & (if (is.na(boott[[q]])) is.na(W$time) else W$time %in% boott[[q]]))
  stopifnot(length(i) == 1)
  W[i, ] <- as.data.frame(fill_boot(as.list(W[i, ]), L7, k, q, "known", "estimated"))
}
# cell 8: rows of the three benchmark estimands only; main columns from the 2,000-replicate OC rows
# ('proposed' = "known", 'proposed_ew' = script-side "estimated"), boot columns from the first-500 rows
# (the bootstrap-replicate "known" of R36 is 'proposed_rect', the only 'known' row of R36_first500 stored).
ess8 <- oc_ess(8)
for (q in names(bootq)) {
  qq <- bootq[[q]]; t <- boott[[q]]
  K <- pick(LO, 8, qq, "proposed", t); E <- pick(LO, 8, qq, "proposed_ew", t)
  r <- na_row()
  r[c("source", "cell_id", "cell", "regime", "N", "scenario", "censoring", "propensity_model", "weights")] <-
    list("OC", 8, "E-400", "About 0.69", 400, scen[["ph"]], "uniform(0, 72) months", "logit e = -0.5 + X", "estimated, ipw")
  r$ess_ratio_control <- ess8[1]; r$ess_ratio_treated <- ess8[2]; r$estimand <- qq; r$time <- t
  r$true_value <- K$truth; r$mean_estimate <- K$truth + K$bias; r$bias <- K$bias; r$bias_mcse <- K$bias_mcse; r$empirical_sd <- K$emp_sd
  r$mean_se_known <- K$mean_se; r$mean_se_estimated <- E$mean_se
  r$coverage_known <- K$coverage; r$coverage_known_mcse <- mcse(K$coverage, K$n_eval)
  r$coverage_estimated <- E$coverage; r$coverage_estimated_mcse <- mcse(E$coverage, E$n_eval)
  r$mean_ci_length_known <- K$mean_length; r$mean_ci_length_estimated <- E$mean_length
  r$replicates <- K$n_rep; r$replicates_evaluable <- K$n_eval; r$n_fail_known <- K$n_fail; r$n_fail_estimated <- E$n_fail
  if (!identical(E$bias, K$bias)) r$bias_over_estimated_rows_if_different <- E$bias
  r <- fill_boot(r, LO, 8, q, if (q == "R36_first500") "proposed_rect" else "proposed", "proposed_ew")
  r$note <- "bootstrap benchmark only; copied from dev/audit/oc/oc_summary.csv: known = 'proposed', estimated = 'proposed_ew' (script-side, equal to the package option on shared replicates); ESS ratios from dev/audit/oc/raw/"
  W <- rbind(W, as.data.frame(r))
}
W <- W[, COLS]
OUTCSV <- file.path(REPO, "vignettes", "articles", "sim_results", "ps_weight_results_by_time.csv")
write.csv(W, OUTCSV, row.names = FALSE, na = "")
cat(sprintf("%s: %d rows\n", OUTCSV, nrow(W)))

# ---- 3. gates ----
G <- character(0); put <- function(...) { G <<- c(G, sprintf(...)); cat(sprintf(...), "\n") }
R <- read.csv(OUTCSV, stringsAsFactors = FALSE)                 # the stored file, read back
# G2: calibration, realised
p31 <- 0.3977; put("G2 calibration (realised, means over 2,000 replicates of the ESS ratio under the fitted stabilized weights):")
g2ok <- TRUE
for (k in c(41, 42)) {
  r <- readRDS(file.path(raw, sprintf("m7_cell%02d.rds", k)))
  e0 <- mean(r$res[, "ess.ratio0"], na.rm = TRUE); e1 <- mean(r$res[, "ess.ratio1"], na.rm = TRUE)
  ok <- e1 >= 0.88 && e1 <= 0.94 && e0 >= 0.95; g2ok <- g2ok && ok
  put("  cell %d: treated %.4f (range of replicates %.3f-%.3f), control %.4f (%.3f-%.3f): %s", k, e1, min(r$res[, "ess.ratio1"], na.rm = TRUE),
      max(r$res[, "ess.ratio1"], na.rm = TRUE), e0, min(r$res[, "ess.ratio0"], na.rm = TRUE), max(r$res[, "ess.ratio0"], na.rm = TRUE), if (ok) "in range" else "OUT OF RANGE")
}
put("G2 %s", if (g2ok) "PASS (ESS part; P(treated) in calibration.txt)" else "FAIL")
# G3 counts
for (k in c(41, 42)) {
  r <- readRDS(file.path(raw, sprintf("m7_cell%02d.rds", k))); M <- r$res
  fc <- grep("\\.fail$", colnames(M), value = TRUE); anyfail <- rowSums(M[, fc, drop = FALSE] == 1, na.rm = TRUE) > 0
  put("G3 cell %d: %d replicates; replicates with any failed call: %d (%.2f%%); ps_weights() failures %d; %.1f s on %d cores; package %s at %s",
      k, nrow(M), sum(anyfail), 100 * mean(anyfail), sum(M[, "ps.fail"] == 1), r$secs, r$cores, r$pkg_version, substr(r$pkg_sha, 1, 7))
}
# G4: truth
T41 <- LN[LN$cell_id == 41 & !is.na(LN$truth), ]; T31 <- L7[L7$cell_id == 31 & !is.na(L7$truth) & L7$estimator %in% c("known", "estimated"), ]
key <- function(D) paste(D$quantity, D$time, D$estimator)
T31 <- T31[key(T31) %in% key(T41), ]; m <- match(key(T41), key(T31))
d4 <- abs(T41$truth - T31$truth[m])
put("G4 cell 41 vs cell 31: %d (quantity, time, option) truths compared (unmatched %d); identical %d; largest |difference| %.3g",
    nrow(T41), sum(is.na(m)), sum(T41$truth == T31$truth[m], na.rm = TRUE), max(d4, na.rm = TRUE))
T42 <- LN[LN$cell_id == 42 & LN$quantity %in% c("delta", "R", "loghr_fh00", "loghr_fh01"), ]
put("G4 cell 42: difference, RMST difference and log HR truths: %d rows, all exactly 0: %s", nrow(T42), all(T42$truth == 0))
g4ok <- !any(is.na(m)) && max(d4) < 1e-12 && all(T42$truth == 0)
put("G4 %s", if (g4ok) "PASS" else "FAIL")
# G5: copied rows against their source, read back from the stored file
fields <- list(true_value = "truth", bias = "bias", bias_mcse = "bias_mcse", empirical_sd = "emp_sd", mean_se_known = c("mean_se", "K"),
               mean_se_estimated = c("mean_se", "E"), coverage_known = c("coverage", "K"), coverage_estimated = c("coverage", "E"),
               mean_ci_length_known = c("mean_length", "K"), mean_ci_length_estimated = c("mean_length", "E"),
               replicates = "n_rep", replicates_evaluable = "n_eval", n_fail_known = c("n_fail", "K"), n_fail_estimated = c("n_fail", "E"))
nrow5 <- 0; nval5 <- 0; bad5 <- character(0)
for (i in which(R$source %in% c("M7", "OC"))) {
  x <- R[i, ]; L <- if (x$source == "M7") L7 else LO
  kn <- if (x$source == "OC") "proposed" else "known"; es <- if (x$source == "OC") (if (x$cell_id == 8) "proposed_ew" else NA) else "estimated"
  t <- if (is.na(x$time)) NA else x$time
  K <- pick(L, x$cell_id, x$estimand, kn, t); E <- if (!is.na(es)) pick(L, x$cell_id, x$estimand, es, t) else NULL
  band <- grepl("_band_", x$estimand)
  for (f in names(fields)) {
    sf <- fields[[f]]; which_ <- if (length(sf) == 2) sf[2] else "K"; sf <- sf[1]
    if (band && !(sf %in% c("coverage", "mean_length", "n_rep", "n_eval", "n_fail"))) next
    if (which_ == "E" && is.null(E)) { if (!is.na(x[[f]])) bad5 <- c(bad5, sprintf("row %d %s not empty", i, f)); next }
    sv <- (if (which_ == "E") E else K)[[sf]]
    nval5 <- nval5 + 1
    if (!identical(as.numeric(x[[f]]), as.numeric(sv))) bad5 <- c(bad5, sprintf("row %d (%s %s %s %s) %s: %s vs %s", i, x$source, x$cell_id, x$estimand, t, f, x[[f]], sv))
  }
  nrow5 <- nrow5 + 1
}
# boot columns of the copied rows
bootmap <- list(coverage_known_bootreps = c("K", "coverage"), coverage_estimated_bootreps = c("E", "coverage"),
                mean_ci_length_known_bootreps = c("K", "mean_length"), mean_ci_length_estimated_bootreps = c("E", "mean_length"),
                coverage_boot = c("B", "coverage"), mean_ci_length_boot = c("B", "mean_length"), coverage_boot_pct = c("P", "coverage"),
                mean_ci_length_boot_pct = c("P", "mean_length"), boot_se_over_sd = c("B", "se_sd"), boot_replicates = c("B", "n_rep"))
nb5 <- 0
for (i in which(!is.na(R$boot_replicates))) {
  x <- R[i, ]; L <- if (x$source == "M7") L7 else LO
  q <- names(bootq)[bootq == x$estimand & (if (is.na(x$time)) is.na(boott) else boott %in% x$time)]
  src <- list(K = pick(L, x$cell_id, q, if (x$source == "M7") "known" else if (q == "R36_first500") "proposed_rect" else "proposed"),
              E = pick(L, x$cell_id, q, if (x$source == "M7") "estimated" else "proposed_ew"), B = pick(L, x$cell_id, q, "boot"), P = pick(L, x$cell_id, q, "boot_pct"))
  for (f in names(bootmap)) { nval5 <- nval5 + 1; sv <- src[[bootmap[[f]][1]]][[bootmap[[f]][2]]]
    if (!identical(as.numeric(x[[f]]), as.numeric(sv))) bad5 <- c(bad5, sprintf("boot row %d %s: %s vs %s", i, f, x[[f]], sv)) }
  nb5 <- nb5 + 1
}
# ESS of M7 rows
for (i in which(R$source == "M7")) {
  x <- R[i, ]; nval5 <- nval5 + 2
  if (!identical(x$ess_ratio_control, pick(L7, x$cell_id, "ess_ratio_control", "known")$mean_value)) bad5 <- c(bad5, sprintf("row %d ess control", i))
  if (!identical(x$ess_ratio_treated, pick(L7, x$cell_id, "ess_ratio_treated", "known")$mean_value)) bad5 <- c(bad5, sprintf("row %d ess treated", i))
}
put("G5 copied rows compared with their source summary (stored CSV read back): %d rows (M7 %d, OC %d; %d of them carry bootstrap columns), %d values; differing %d",
    nrow5, sum(R$source == "M7"), sum(R$source == "OC"), nb5, nval5, length(bad5))
if (length(bad5)) put("  %s", head(bad5, 20))
# consistency with the inventory CSV (read-only record, tolerance 1e-12 relative) on the shared statistical columns
INV <- read.csv(file.path(.here, "ps_weight_simulation_table.csv"), stringsAsFactors = FALSE)
cmpc <- setdiff(COLS[14:42], c("bias_over_estimated_rows_if_different"))
ki <- function(D) paste(D$source, D$cell_id, D$estimand, D$time)
mi <- match(ki(R[R$source != "sim", ]), ki(INV)); Rc <- R[R$source != "sim", ]
dif <- 0; ncmp <- 0; nb510 <- 0
bootcols <- COLS[32:42]
for (cc in c(cmpc, "ess_ratio_control", "ess_ratio_treated")) {
  a <- as.numeric(Rc[[cc]]); b <- as.numeric(INV[[cc]][mi])
  if (cc %in% bootcols) {                                  # true-weight bootstrap of cells 5, 10: in the inventory, not in this file
    x510 <- Rc$cell_id %in% c(5, 10) & Rc$source == "OC"
    nb510 <- nb510 + sum(x510 & is.na(a) & !is.na(b)); a[x510] <- NA; b[x510] <- NA
  }
  s <- !(is.na(a) & is.na(b)); ncmp <- ncmp + sum(s)
  dif <- dif + sum(is.na(a[s]) != is.na(b[s]) | abs(a[s] - b[s]) > 1e-12 * pmax(1, abs(b[s])), na.rm = TRUE)
}
put("G5 cross-check with the inventory CSV: %d rows matched (unmatched %d), %d values compared, %d differ beyond 1e-12 relative; %d inventory values not carried (the true-weight bootstrap columns of OC cells 5 and 10, outside the benchmark cells 8, 38, 39)",
    sum(!is.na(mi)), sum(is.na(mi)), ncmp, dif, nb510)
# realised P(treated): the data of the new cells and of cell 31 regenerated (data only, no package call)
pt <- function(k) mean(vapply(1:2000, function(r) mean(m7_data9(k, r)$treat), numeric(1)))
m7_data9 <- function(k, r) {
  if (k == 31) return(make_data(8, r)$df)
  set.seed(1e6 * k + r); X <- rnorm(400); H <- -log(pnorm(0.6 * X + 0.8 * rnorm(400))); e <- plogis(PA + PB * X)
  data.frame(treat = rbinom(400, 1, e))
}
p <- c(`31` = pt(31), `41` = pt(41), `42` = pt(42))
put("G2 realised P(treated), mean over 2,000 replicates: cell 31 %.4f, cell 41 %.4f, cell 42 %.4f; largest difference from cell 31 %.4f (limit 0.02)",
    p[1], p[2], p[3], max(abs(p[2:3] - p[1])))
put("G5 %s", if (length(bad5) == 0) "PASS" else "FAIL")
# requirement 3 of M7 in the new cells
x3 <- req$r3; put("M7 requirement 3 in cells 41, 42 (estimates identical under the two options): %d of %d replicates failing; %s", sum(x3$nbad), sum(x3$nrep), paste(x3$note, collapse = "; "))
writeLines(G, file.path(.here, "ps9_gates.txt"))
