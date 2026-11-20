#!/usr/bin/env Rscript
# 02_recompute.R -- TASK_rmst_readoff_truth.md, section 2, steps 2-3; gates G3-G6.
# Recomputes every RMST row (estimand "R") of vignettes/articles/sim_results/ps_weight_results_by_time.csv against
# the read-off truth R(t_j) of 01_recover_tj.R (tj.rds), from the stored per-replicate estimates and limits, with
# pw() of the summary scripts (dev/audit/oc/oc_summarise.R:32-42 = dev/audit/stage2e/oc_summarise.R:36-46, copied).
# Adds the column mean_truth_readoff after true_value (empty on non-RMST rows). Non-RMST lines are carried over as
# text with one empty field inserted. No package call; no simulation.
#
#   Rscript --vanilla dev/audit/rmst_readoff/02_recompute.R [--write]
#   without --write: computes and checks, writes nothing to the CSV

.here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
AUD <- normalizePath(file.path(.here, "..")); REPO <- normalizePath(file.path(AUD, "..", ".."))
WRITE <- "--write" %in% commandArgs(TRUE)
source(file.path(AUD, "stage2e", "oc_common.R"))
CSV <- file.path(REPO, "vignettes", "articles", "sim_results", "ps_weight_results_by_time.csv")
TJ <- readRDS(file.path(.here, "tj.rds"))
Q <- qnorm(0.975)

# pw(): verbatim from dev/audit/stage2e/oc_summarise.R:36-46 (identical to dev/audit/oc/oc_summarise.R:32-42)
pw <- function(est, lo, hi, se, truth, ev, fail, truth_nom = truth[1]) {
  fail <- !is.na(fail) & fail == 1
  finI <- is.finite(lo) & is.finite(hi)
  cover <- finI & lo <= truth & truth <= hi
  fe <- ev & is.finite(est)
  err <- (est - truth)[fe]
  list(n_rep = length(est), n_eval = sum(ev), n_fail = sum(fail), n_nonfinite = sum(ev & !fail & !finI),
       truth = truth_nom, bias = mean(err), bias_mcse = sd(err) / sqrt(length(err)), emp_sd = sd(est[fe]),
       mean_se = mean(se[ev & is.finite(se)]), se_sd = mean(se[ev & is.finite(se)]) / sd(est[fe]),
       coverage = if (is.na(truth_nom)) NA else mean(cover[ev]), mean_length = mean((hi - lo)[ev & finI]))
}
mcse <- function(p, n) sqrt(p * (1 - p) / n)      # as dev/audit/ps_by_time/ps9_summarise.R:64

rawfile <- function(src, k) switch(src,
  OC = file.path(AUD, "oc", "raw", sprintf("proposed_cell%02d.rds", k)),
  M7 = file.path(AUD, "stage2e", "raw", sprintf("m7_cell%02d.rds", k)),
  sim = file.path(AUD, "ps_by_time", "raw", sprintf("m7_cell%02d.rds", k)))
bootfile <- function(src, k) if (src == "OC") file.path(AUD, "oc", "raw", sprintf("boot_cell%02d.rds", k)) else
  file.path(AUD, "stage2e", "raw", sprintf("boot_cell%02d.rds", k))
RAW <- list(); getM <- function(src, k) { key <- paste(src, k); if (is.null(RAW[[key]])) RAW[[key]] <<- readRDS(rawfile(src, k)); RAW[[key]] }

LOG <- character(0); say <- function(...) { s <- sprintf(...); LOG <<- c(LOG, s); cat(s, "\n") }
say("02_recompute.R  %s  R %s  write=%s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), getRversion(), WRITE)

# ---- the stored CSV: raw lines and parsed values ----
lines <- readLines(CSV)
splitq <- function(s) strsplit(s, ',(?=(?:[^"]*"[^"]*")*[^"]*$)', perl = TRUE)[[1]]
hdr <- splitq(lines[1]); body <- lines[-1]
OLD <- read.csv(CSV, stringsAsFactors = FALSE)
stopifnot(nrow(OLD) == length(body), identical(gsub('"', "", hdr), names(OLD)))
iR <- which(OLD$estimand == "R")
say("CSV: %d rows, %d RMST rows (estimand R), %d other rows", nrow(OLD), length(iR), nrow(OLD) - length(iR))

# ---- recompute ----
NEW <- OLD
NEW$mean_truth_readoff <- NA_real_
cmp <- list(); bootinfo <- list()
for (i in iR) {
  x <- OLD[i, ]; src <- x$source; k <- x$cell_id; t <- x$time
  st <- getM(src, k); M <- st$res
  kn <- if (src == "OC") "proposed" else "known"
  es <- if (src == "OC") (if (k == 8) "proposed_ew" else NA) else "estimated"
  rec <- TJ[[sprintf("%s_%02d", src, k)]]; stopifnot(!is.null(rec), identical(unname(rec$chk_ok), rep(TRUE, nrow(M))))
  truth <- rec$Rtj[, as.character(t)]
  ev <- M[, paste0("evaluable.", t)] == 1
  P <- function(e, f) M[, sprintf("%s.R.%s.%s", e, f, t)]
  run <- function(e) pw(P(e, "est"), P(e, "lo"), P(e, "hi"), P(e, "se"), truth, ev, M[, paste0(e, ".rmst.fail")], x$true_value)
  K <- run(kn); E <- if (!is.na(es)) run(es) else NULL
  est <- P(kn, "est"); fe <- ev & is.finite(est)
  NEW$mean_truth_readoff[i] <- mean(truth[fe])
  NEW$mean_estimate[i] <- mean(est[fe])
  NEW$bias[i] <- K$bias; NEW$bias_mcse[i] <- K$bias_mcse
  NEW$empirical_sd[i] <- K$emp_sd; NEW$mean_se_known[i] <- K$mean_se
  NEW$coverage_known[i] <- K$coverage; NEW$coverage_known_mcse[i] <- mcse(K$coverage, K$n_eval)
  NEW$mean_ci_length_known[i] <- K$mean_length
  if (!is.null(E)) {
    NEW$mean_se_estimated[i] <- E$mean_se
    NEW$coverage_estimated[i] <- E$coverage; NEW$coverage_estimated_mcse[i] <- mcse(E$coverage, E$n_eval)
    NEW$mean_ci_length_estimated[i] <- E$mean_length
    NEW$bias_over_estimated_rows_if_different[i] <- if (!identical(E$bias, K$bias)) E$bias else NA
  }
  # also against the nominal R(t) (for the comparative table)
  Kn <- pw(P(kn, "est"), P(kn, "lo"), P(kn, "hi"), P(kn, "se"), rep(x$true_value, nrow(M)), ev, M[, paste0(kn, ".rmst.fail")], x$true_value)
  En <- if (!is.na(es)) pw(P(es, "est"), P(es, "lo"), P(es, "hi"), P(es, "se"), rep(x$true_value, nrow(M)), ev, M[, paste0(es, ".rmst.fail")], x$true_value) else NULL
  cmp[[length(cmp) + 1]] <- data.frame(row = i, source = src, cell_id = k, cell = x$cell, time = t, true_value = x$true_value,
    mean_truth_readoff = NEW$mean_truth_readoff[i], gap = x$true_value - NEW$mean_truth_readoff[i],
    n_eval = K$n_eval, n_fe = sum(fe), n_fail_K = K$n_fail, n_nonfinite_K = K$n_nonfinite,
    mean_est_all_eval = mean(est[ev], na.rm = TRUE),
    bias_stored = x$bias, bias_Rt = Kn$bias, bias_Rtj = K$bias,
    covk_stored = x$coverage_known, covk_Rt = Kn$coverage, covk_Rtj = K$coverage,
    cove_stored = x$coverage_estimated, cove_Rt = if (is.null(En)) NA else En$coverage, cove_Rtj = if (is.null(E)) NA else E$coverage,
    E_bias_identical = if (is.null(E)) NA else identical(E$bias, K$bias),
    stringsAsFactors = FALSE)

  # bootstrap benchmark columns (cells 8, 38, 39, month 36)
  if (!is.na(x$boot_replicates)) {
    bt <- readRDS(bootfile(src, k)); Bm <- bt$res; nb <- nrow(Bm)
    stopifnot(identical(unname(bt$chk), unname(st$chk[seq_len(nb)])))
    s5 <- seq_len(nb)
    evb <- if (src == "OC") rep(TRUE, nb) else ev[s5]        # oc_summarise.R:148; stage2e/oc_summarise.R:169,172
    knb <- if (src == "OC") "proposed_rect" else "known"     # ps9_summarise.R:146
    pb <- function(e, tr) pw(P(e, "est")[s5], P(e, "lo")[s5], P(e, "hi")[s5], P(e, "se")[s5], tr, evb, M[s5, paste0(e, ".rmst.fail")], x$true_value)
    Kb <- pb(knb, truth[s5]); Eb <- pb(es, truth[s5])
    Kb0 <- pb(knb, rep(x$true_value, nb)); Eb0 <- pb(es, rep(x$true_value, nb))
    NEW$coverage_known_bootreps[i] <- Kb$coverage; NEW$coverage_estimated_bootreps[i] <- Eb$coverage
    NEW$mean_ci_length_known_bootreps[i] <- Kb$mean_length; NEW$mean_ci_length_estimated_bootreps[i] <- Eb$mean_length
    # the bootstrap's own intervals (est.R36 from survfit(..., rmean = 36): integral of each curve over [0, 36])
    est <- Bm[, "est.R36"]; s <- Bm[, "sd.R36"]; nof <- rep(0, nb)
    bn_t <- pw(est, est - Q * s, est + Q * s, s, rep(x$true_value, nb), evb, nof, x$true_value)
    bp_t <- pw(est, Bm[, "pct_lo.R36"], Bm[, "pct_hi.R36"], s, rep(x$true_value, nb), evb, nof, x$true_value)
    bn_j <- pw(est, est - Q * s, est + Q * s, s, truth[s5], evb, nof, x$true_value)
    bp_j <- pw(est, Bm[, "pct_lo.R36"], Bm[, "pct_hi.R36"], s, truth[s5], evb, nof, x$true_value)
    bootinfo[[length(bootinfo) + 1]] <- data.frame(row = i, cell_id = k, nb = nb, chk_same = TRUE,
      covk_stored = x$coverage_known_bootreps, covk_Rt = Kb0$coverage, covk_Rtj = Kb$coverage,
      cove_stored = x$coverage_estimated_bootreps, cove_Rt = Eb0$coverage, cove_Rtj = Eb$coverage,
      lenk_same = isTRUE(abs(Kb$mean_length - x$mean_ci_length_known_bootreps) <= 1e-12),
      lene_same = isTRUE(abs(Eb$mean_length - x$mean_ci_length_estimated_bootreps) <= 1e-12),
      boot_stored = x$coverage_boot, boot_Rt = bn_t$coverage, boot_Rtj = bn_j$coverage,
      bootpct_stored = x$coverage_boot_pct, bootpct_Rt = bp_t$coverage, bootpct_Rtj = bp_j$coverage,
      bootlen_same = isTRUE(abs(bn_t$mean_length - x$mean_ci_length_boot) <= 1e-12), stringsAsFactors = FALSE)
  }
}
CMP <- do.call(rbind, cmp); BOOT <- do.call(rbind, bootinfo)
say("recomputed %d RMST rows; replicates with evaluable month but no finite estimate: %d; failed RMST calls among evaluable: %d; non-finite intervals: %d",
    nrow(CMP), sum(CMP$n_eval - CMP$n_fe), sum(CMP$n_fail_K), sum(CMP$n_nonfinite_K))
say("'estimated' bias identical to 'known' bias (estimates identical under the two options) in %d of %d rows with an 'estimated' option",
    sum(CMP$E_bias_identical, na.rm = TRUE), sum(!is.na(CMP$E_bias_identical)))
say("bootstrap benchmark rows: per-replicate checksums of boot_cell*.rds equal those of the main files in all %d rows", nrow(BOOT))
for (j in seq_len(nrow(BOOT))) with(BOOT[j, ], say(
  "  cell %d: known_bootreps stored %.4f, R(t) %.4f, R(t_j) %.4f; estimated_bootreps stored %.4f, R(t) %.4f, R(t_j) %.4f; lengths unchanged %s/%s; boot normal stored %.4f, R(36) %.4f, R(t_j) %.4f; boot percentile stored %.4f, R(36) %.4f, R(t_j) %.4f; boot length reproduced %s",
  cell_id, covk_stored, covk_Rt, covk_Rtj, cove_stored, cove_Rt, cove_Rtj, lenk_same, lene_same, boot_stored, boot_Rt, boot_Rtj,
  bootpct_stored, bootpct_Rt, bootpct_Rtj, bootlen_same))
say("  coverage_boot, coverage_boot_mcse, coverage_boot_pct, boot_se_over_sd: left unchanged (est.R36 of oc_boot.R is survfit(..., rmean = 36), integrated to 36 itself, not read at t_j; stored values reproduced against R(36): %s)",
    all(abs(BOOT$boot_stored - BOOT$boot_Rt) <= 1e-12 & abs(BOOT$bootpct_stored - BOOT$bootpct_Rt) <= 1e-12))

# ---- gates ----
tol <- 1e-12; near <- function(a, b) (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & abs(a - b) <= tol)
# G3: OC cells 5, 10 against the stored OC summary (dev/audit/oc/oc_summary.csv, read only) and the stored CSV
LO <- read.csv(file.path(AUD, "oc", "oc_summary.csv"), stringsAsFactors = FALSE)
g3 <- 0; g3n <- 0; g3bad <- character(0)
for (i in iR[OLD$source[iR] == "OC" & OLD$cell_id[iR] %in% c(5, 10)]) {
  s <- LO[LO$cell_id == OLD$cell_id[i] & LO$quantity == "R" & LO$estimator == "proposed" & !is.na(LO$time) & LO$time == OLD$time[i], ]
  stopifnot(nrow(s) == 1)
  pairs <- list(bias = c(NEW$bias[i], s$bias), bias_from_means = c(NEW$mean_estimate[i] - NEW$mean_truth_readoff[i], s$bias),
                bias_mcse = c(NEW$bias_mcse[i], s$bias_mcse), coverage_known = c(NEW$coverage_known[i], s$coverage),
                coverage_known_mcse = c(NEW$coverage_known_mcse[i], mcse(s$coverage, s$n_eval)),
                bias_csv = c(NEW$bias[i], OLD$bias[i]), coverage_csv = c(NEW$coverage_known[i], OLD$coverage_known[i]),
                coverage_mcse_csv = c(NEW$coverage_known_mcse[i], OLD$coverage_known_mcse[i]))
  for (nm in names(pairs)) { g3n <- g3n + 1; d <- abs(diff(pairs[[nm]])); g3 <- max(g3, d)
    if (!(d <= tol)) g3bad <- c(g3bad, sprintf("cell %d month %d %s: %.17g vs %.17g", OLD$cell_id[i], OLD$time[i], nm, pairs[[nm]][1], pairs[[nm]][2])) }
}
G3 <- length(g3bad) == 0
say("G3 OC cells 5, 10 (22 rows): %d comparisons (bias, mean_estimate - mean_truth_readoff, bias MCSE, coverage, coverage MCSE; against oc_summary.csv and the stored CSV); largest |difference| %.3g: %s",
    g3n, g3, if (G3) "PASS" else "FAIL"); if (!G3) say("  %s", g3bad)
# G4: OC cell 8 against M7 cell 31, RMST month 36
a <- which(NEW$source == "OC" & NEW$cell_id == 8 & NEW$estimand == "R"); b <- which(NEW$cell_id == 31 & NEW$estimand == "R" & NEW$time == 36)
f4 <- c("mean_estimate", "bias", "bias_mcse", "coverage_known", "coverage_estimated")
d4 <- vapply(f4, function(f) abs(NEW[[f]][a] - NEW[[f]][b]), 0)
G4 <- all(d4 <= tol)
say("G4 OC cell 8 against M7 cell 31, RMST month 36: %s: %s", paste(sprintf("%s %.3g", f4, d4), collapse = "; "), if (G4) "PASS" else "FAIL")
say("   values: %s", paste(sprintf("%s %.15g / %.15g", f4, unlist(NEW[a, f4]), unlist(NEW[b, f4])), collapse = "; "))
# G5: means of the stored estimates; unchanged columns
g5n <- 0; g5bad <- character(0)
for (i in iR) {
  x <- OLD[i, ]; M <- getM(x$source, x$cell_id)$res; kn <- if (x$source == "OC") "proposed" else "known"
  e <- M[, sprintf("%s.R.est.%s", kn, x$time)]; ev <- M[, paste0("evaluable.", x$time)] == 1
  chk <- list(mean_estimate = c(NEW$mean_estimate[i], mean(e[ev & !is.na(e)])))
  for (f in c("empirical_sd", "mean_se_known", "mean_se_estimated", "mean_ci_length_known", "mean_ci_length_estimated",
              "replicates", "replicates_evaluable", "n_fail_known", "n_fail_estimated", "true_value",
              "mean_ci_length_known_bootreps", "mean_ci_length_estimated_bootreps"))
    chk[[f]] <- c(NEW[[f]][i], OLD[[f]][i])
  for (f in names(chk)) { g5n <- g5n + 1; if (!near(chk[[f]][1], chk[[f]][2])) g5bad <- c(g5bad, sprintf("row %d cell %d month %d %s: %.17g vs %.17g", i, x$cell_id, x$time, f, chk[[f]][1], chk[[f]][2])) }
}
G5 <- length(g5bad) == 0
say("G5 %d RMST rows, %d values (mean_estimate against the mean of the stored per-replicate estimates over evaluable replicates; empirical SD, mean SEs, mean lengths, counts, true_value, bootreps lengths against the stored CSV): %d differ beyond 1e-12: %s",
    length(iR), g5n, length(g5bad), if (G5) "PASS" else "FAIL"); if (!G5) say("  %s", head(g5bad, 30))

# ---- the new file: mean_truth_readoff after true_value ----
COLS <- append(names(OLD), "mean_truth_readoff", after = match("true_value", names(OLD)))
NEW <- NEW[, COLS]
tc <- textConnection("newtxt", "w", local = TRUE); write.csv(NEW, tc, row.names = FALSE, na = ""); close(tc)
pos <- match("true_value", names(OLD))                            # insert after field `pos`
ins <- function(s) { f <- splitq(s); stopifnot(length(f) == length(hdr)); paste(c(f[seq_len(pos)], "", f[-seq_len(pos)]), collapse = ",") }
# RMST rows: a recomputed value within 1e-12 of the stored one keeps the stored text (the stored MCSEs of copied rows were
# computed from coverages read back from a CSV, so they can differ in the 16th digit; such differences are not changes)
merge_rmst <- function(i) {
  fn <- splitq(newtxt[i + 1]); fo <- splitq(ins(body[i])); nn <- names(NEW)
  for (j in seq_along(fn)) {
    a <- suppressWarnings(as.numeric(gsub('"', "", fo[j]))); b <- suppressWarnings(as.numeric(gsub('"', "", fn[j])))
    if (fn[j] != fo[j] && is.finite(a) && is.finite(b) && abs(a - b) <= 1e-12 * max(1, abs(a))) fn[j] <- fo[j]
  }
  paste(fn, collapse = ",")
}
out <- c(newtxt[1], vapply(seq_along(body), function(i) if (i %in% iR) merge_rmst(i) else ins(body[i]), ""))
stopifnot(identical(splitq(out[1]), append(hdr, '"mean_truth_readoff"', after = pos)))
# G6: non-RMST lines byte-identical apart from the added empty field
drop_field <- function(s) { f <- splitq(s); paste(f[-(pos + 1)], collapse = ",") }
non <- setdiff(seq_along(body), iR)
same6 <- vapply(non, function(i) identical(drop_field(out[i + 1]), body[i]) && splitq(out[i + 1])[pos + 1] == "", TRUE)
# and write.csv's own rendering of these rows equals the old text (round trip), for the record
rt6 <- vapply(non, function(i) identical(newtxt[i + 1], out[i + 1]), TRUE)
G6 <- all(same6) && length(non) == nrow(OLD) - length(iR)
say("G6 non-RMST rows: %d of %d byte-identical to before apart from the added empty field (write.csv round trip also identical in %d of %d): %s",
    sum(same6), length(non), sum(rt6), length(non), if (G6) "PASS" else "FAIL")
# RMST rows: fields not recomputed keep their text
keepf <- setdiff(names(OLD), c("mean_estimate", "bias", "bias_mcse", "empirical_sd", "mean_se_known", "mean_se_estimated", "coverage_known",
                               "coverage_known_mcse", "coverage_estimated", "coverage_estimated_mcse", "mean_ci_length_known", "mean_ci_length_estimated",
                               "bias_over_estimated_rows_if_different", "coverage_known_bootreps", "coverage_estimated_bootreps",
                               "mean_ci_length_known_bootreps", "mean_ci_length_estimated_bootreps"))
kt <- vapply(iR, function(i) identical(splitq(drop_field(out[i + 1]))[match(keepf, names(OLD))], splitq(body[i])[match(keepf, names(OLD))]), TRUE)
say("RMST rows: the %d columns not recomputed keep their text in %d of %d rows", length(keepf), sum(kt), length(iR))
# changed text by column in RMST rows
chg <- vapply(names(OLD), function(f) sum(vapply(iR, function(i) splitq(drop_field(out[i + 1]))[match(f, names(OLD))] != splitq(body[i])[match(f, names(OLD))], TRUE)), 0)
say("RMST rows, columns whose text changed (rows): %s", paste(sprintf("%s %d", names(chg)[chg > 0], chg[chg > 0]), collapse = "; "))

saveRDS(list(CMP = CMP, BOOT = BOOT), file.path(.here, "recompute.rds"))
write.csv(CMP, file.path(.here, "rmst_rows_old_new.csv"), row.names = FALSE)
ok <- G3 && G4 && G5 && G6
if (WRITE && ok) { writeLines(out, CSV); say("written: %s (%d lines)", CSV, length(out)) } else say("CSV not written (write=%s, gates %s)", WRITE, if (ok) "pass" else "FAIL")
writeLines(LOG, file.path(.here, if (WRITE) "02_recompute.txt" else "02_recompute_dryrun.txt"))
if (!ok) quit(status = 1)
