#!/usr/bin/env Rscript
# oc_summarise.R -- Stage 2e, M7 (TASK_km_fixes_stage2e.md, section 6): builds m7_required.txt, m7_summary.csv and
# m7_tables.md from the per-replicate results in dev/audit/stage2e/raw/.
# Adapted copy of dev/audit/oc/oc_summarise.R. Kept from it: the summary definitions (pw(): coverage over the
# evaluable replicates, an errored call or a non-finite interval counts as not covering; band coverage; the
# log hazard ratio and log-rank rows). Changed: the cells are those of M7; the estimators are "known" and
# "estimated" (A, B, C in cell 40; boot and boot_pct in cells 38, 39); the four requirements are checked.
#
#   Rscript dev/audit/stage2e/oc_summarise.R [--raw=raw]
#
# m7_summary.csv is long: one row per cell x quantity x time x estimator. `time` is NA for bands, log HR, log-rank
# and diagnostic rows. Quantities: S0, S1, delta, R (pointwise, 11 time points); delta_band_3_42, delta_band_any,
# R_band_3_42, R_band_any; loghr_fh00, loghr_fh01; logrank_fh00, logrank_fh01, logrank_strat_fh00,
# logrank_strat_fh01 (column `reject`); ess_ratio_control, ess_ratio_treated, cor_lr_dhat18_estimated,
# cor_lr_dhat18_empirical (column `mean_value`); and, for the bootstrap benchmark (first 500 replicates of cells
# 38, 39): delta18_first500, R36_first500, loghr_first500.

.here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(.here, "oc_common.R"))
opt <- parse_args(list(raw = "raw"))
raw <- file.path(.here, opt$raw); OUT <- .here
Q <- qnorm(0.975); NT <- length(tt)
rd <- function(k) { f <- file.path(raw, sprintf("m7_cell%02d.rds", k)); if (file.exists(f)) readRDS(f) else NULL }
has <- function(M, nm) nm %in% colnames(M)
g11 <- function(M, prefix) M[, paste0(prefix, ".", tt), drop = FALSE]

rows <- list()
add <- function(...) rows[[length(rows) + 1]] <<- data.frame(..., stringsAsFactors = FALSE)
blank <- list(n_rep = NA, n_eval = NA, n_fail = NA, n_nonfinite = NA, truth = NA, bias = NA, bias_mcse = NA,
              emp_sd = NA, mean_se = NA, se_sd = NA, coverage = NA, mean_length = NA, reject = NA, mean_calpha = NA, mean_value = NA)
mk <- function(k, quantity, time, estimator, v) {
  b <- blank; b[names(v)] <- v
  do.call(add, c(list(cell_id = k, cell = m7$cell[m7$k == k], quantity = quantity, time = time, estimator = estimator), b))
}
# pointwise metrics over replicates (the evaluation's definition)
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
pw_block <- function(k, quantity, estimator, M, est, lo, hi, se, truth_nom, EV, fail) {
  for (i in seq_len(NT)) {
    mk(k, quantity, tt[i], estimator, pw(est[, i], lo[, i], hi[, i], se[, i], rep(truth_nom[i], nrow(M)), EV[, i] == 1, fail, truth_nom[i]))
  }
}
band_block <- function(k, quantity, estimator, M, prefix) {
  ca <- M[, paste0(prefix, ".calpha")]; fail <- M[, paste0(prefix, ".fail")] == 1
  bad <- !fail & (!is.finite(ca) | M[, paste0(prefix, ".nonfin")] > 0)
  okc <- !fail & is.finite(ca)
  mk(k, paste0(quantity, "_band_3_42"), NA, estimator,
     list(n_rep = nrow(M), n_eval = nrow(M), n_fail = sum(fail), n_nonfinite = sum(bad),
          coverage = mean(okc & M[, paste0(prefix, ".cov_in")] == 1),
          mean_length = mean(M[okc, paste0(prefix, ".width")], na.rm = TRUE), mean_calpha = mean(ca[okc])))
  mk(k, paste0(quantity, "_band_any"), NA, estimator,
     list(n_rep = nrow(M), n_eval = nrow(M), n_fail = sum(fail), n_nonfinite = sum(bad),
          coverage = mean(okc & M[, paste0(prefix, ".cov_any")] == 1), mean_calpha = mean(ca[okc])))
}
cox_block <- function(k, g, estimator, M, prefix) {
  p <- function(f) M[, paste0(prefix, ".cox", g, ".", f)]
  tr <- m7$loghr_truth[m7$k == k]
  v <- pw(p("b"), p("lo"), p("hi"), p("se"), rep(tr, nrow(M)), rep(TRUE, nrow(M)), p("fail"), tr)
  if (is.na(tr)) v$bias <- v$bias_mcse <- NA
  mk(k, sprintf("loghr_fh0%d", g), NA, estimator, v)
  z <- p("z")
  if (any(is.finite(z))) mk(k, sprintf("logrank_fh0%d", g), NA, estimator,
                            list(n_rep = nrow(M), n_eval = sum(is.finite(z)), n_fail = sum(p("fail") == 1, na.rm = TRUE),
                                 reject = mean(abs(z[is.finite(z)]) > 1.96)))
}

# ---- requirements ----
req <- list(r1 = NULL, r2 = NULL, r3 = NULL, r4 = NULL)
acc <- function(id, nbad, maxd, nrep, note = "") req[[id]] <<- rbind(req[[id]], data.frame(nbad = nbad, maxd = maxd, nrep = nrep, note = note, stringsAsFactors = FALSE))
timing <- list(); psfail <- list(); kdone <- integer(0)

for (k in m7$k) {
  r <- rd(k); if (is.null(r)) next
  kdone <- c(kdone, k); M <- r$res; cl <- m7[m7$k == k, ]; n <- nrow(M)
  timing[[length(timing) + 1]] <- data.frame(cell_id = k, cell = cl$cell, reps = n, secs = r$secs, cores = r$cores, sha = substr(r$pkg_sha, 1, 7))
  psfail[[length(psfail) + 1]] <- data.frame(cell_id = k, ps_fail = sum(M[, "ps.fail"] == 1))
  EV <- g11(M, "evaluable")
  sc <- cl$sc
  ests <- if (cl$sub) c("A", "B", "C") else if (cl$strat) character(0) else c("known", "estimated")
  for (e in ests) {
    fail <- M[, paste0(e, ".km.fail")]
    S0 <- exp(-Lam(tt, 0, sc)); S1 <- exp(-Lam(tt, 1, sc))
    for (a in 0:1) {
      est <- g11(M, sprintf("%s.S%d.est", e, a)); se <- g11(M, sprintf("%s.S%d.se", e, a))
      pw_block(k, paste0("S", a), e, M, est, est - Q * se, est + Q * se, se, if (a == 0) S0 else S1, EV, fail)
    }
    pw_block(k, "delta", e, M, g11(M, paste0(e, ".d.est")), g11(M, paste0(e, ".d.lo")), g11(M, paste0(e, ".d.hi")), g11(M, paste0(e, ".d.se")),
             delta_sc(tt, sc), EV, fail)
    if (cl$band) band_block(k, "delta", e, M, paste0(e, ".dband"))
    if (cl$rmst) {
      pw_block(k, "R", e, M, g11(M, paste0(e, ".R.est")), g11(M, paste0(e, ".R.lo")), g11(M, paste0(e, ".R.hi")), g11(M, paste0(e, ".R.se")),
               R_sc(tt, sc), EV, M[, paste0(e, ".rmst.fail")])
      band_block(k, "R", e, M, paste0(e, ".rband"))
    }
    for (g in if (cl$fh01) c(0, 1) else 0) cox_block(k, g, e, M, e)
  }
  if (cl$strat) for (g in c(0, 1)) for (e in c("known", "estimated")) {
    z <- M[, sprintf("%s.slr%d.z", e, g)]
    mk(k, sprintf("logrank_strat_fh0%d", g), NA, e, list(n_rep = n, n_eval = sum(is.finite(z)), n_fail = sum(M[, sprintf("%s.slr%d.fail", e, g)] == 1, na.rm = TRUE),
                                                         reject = mean(abs(z[is.finite(z)]) > 1.96)))
  }
  # effective-sample-size ratio by arm
  en <- if (cl$sub) list(A = c("A.ess.ratio0", "A.ess.ratio1"), B = c("ess.ratio0", "ess.ratio1"), C = c("ess.ratio0", "ess.ratio1"))
        else list(known = c("ess.ratio0", "ess.ratio1"), estimated = c("ess.ratio0", "ess.ratio1"))
  for (e in names(en)) {
    mk(k, "ess_ratio_control", NA, e, list(n_rep = n, n_eval = sum(is.finite(M[, en[[e]][1]])), mean_value = mean(M[, en[[e]][1]], na.rm = TRUE)))
    mk(k, "ess_ratio_treated", NA, e, list(n_rep = n, n_eval = sum(is.finite(M[, en[[e]][2]])), mean_value = mean(M[, en[[e]][2]], na.rm = TRUE)))
  }
  # correlation of the log-rank score with the difference at 18 (cells 31, 32)
  if (cl$both) {
    ok <- is.finite(M[, "wd.score"]) & is.finite(M[, "wd.dhat"]) & is.finite(M[, "wd.cor"])
    mk(k, "cor_lr_dhat18_estimated", NA, "estimated", list(n_rep = n, n_eval = sum(ok), mean_value = mean(M[ok, "wd.cor"])))
    mk(k, "cor_lr_dhat18_empirical", NA, "estimated", list(n_rep = n, n_eval = sum(ok), mean_value = stats::cor(M[ok, "wd.score"], M[ok, "wd.dhat"])))
  }

  # requirement 3, every cell: identical estimates under the two options
  c3 <- grep("^req3\\.", colnames(M), value = TRUE)
  X3 <- M[, c3, drop = FALSE]
  designed <- colSums(!is.na(X3)) > 0                                  # fields evaluated in this cell
  acc("r3", sum(rowSums(X3 == 0, na.rm = TRUE) > 0), sum(X3 == 0, na.rm = TRUE), n,
      sprintf("cell %d: %d comparisons not made (a call failed)", k, sum(is.na(X3[, designed, drop = FALSE]))))
  # requirement 1
  if (cl$repro) {
    f8 <- c("req1.km", "req1.band", "req1.rmst")
    X <- M[, f8, drop = FALSE]
    cx <- NULL
    for (g in 0:1) {
      al <- M[, paste0("req1.aligned", g)] == 1 & !is.na(M[, paste0("req1.aligned", g)])
      v <- pmax(M[, paste0("req1.cox", g, ".se")], M[, paste0("req1.cox", g, ".z")])
      v[!al] <- 0                                                      # only where the script's alignment confirmation passes
      cx <- cbind(cx, v)
    }
    X <- cbind(X, cx)
    bad <- rowSums(is.na(X) | X > 1e-8) > 0
    acc("r1", sum(bad), max(X, na.rm = TRUE), n, sprintf("requirement 1, cell %d: %d of %d replicates failing; largest difference of the Kaplan-Meier, band and RMST quantities %.3g, of the Cox / log-rank quantities %.3g; the script's alignment confirmation passed in %d (FH(0,0)) and %d (FH(0,1)) of %d replicates", k,
                                                         sum(bad), n, max(M[, f8], na.rm = TRUE), max(cx, na.rm = TRUE),
                                                         sum(M[, "req1.aligned0"] == 1, na.rm = TRUE), sum(M[, "req1.aligned1"] == 1, na.rm = TRUE), n))
  }
  # requirement 2
  if (cl$both) {
    X <- M[, grep("^req2\\.", colnames(M)), drop = FALSE]
    acc("r2", sum(rowSums(is.na(X) | X > 1e-8) > 0), max(X, na.rm = TRUE), n)
  }
  # requirement 4
  if (cl$coxph) {
    bad <- is.na(M[, "req4.z"]) | M[, "req4.z"] > 1e-8 | is.na(M[, "req4.se"]) | M[, "req4.se"] > 1e-6
    acc("r4", sum(bad), max(M[, "req4.z"], na.rm = TRUE), n, sprintf("largest difference of the SE (limit 1e-6): %.3g", max(M[, "req4.se"], na.rm = TRUE)))
  }
}

# ---- bootstrap benchmark (cells 38, 39; first 500 replicates) ----
boot_done <- all(file.exists(file.path(raw, sprintf("boot_cell%02d.rds", c(38, 39)))))
if (boot_done) for (k in c(38, 39)) {
  bt <- readRDS(file.path(raw, sprintf("boot_cell%02d.rds", k))); Bm <- bt$res; nb <- nrow(Bm)
  r <- rd(k); if (is.null(r)) next
  M <- r$res[seq_len(nb), , drop = FALSE]; sc <- m7$sc[m7$k == k]
  stopifnot(identical(unname(bt$chk), unname(r$chk[seq_len(nb)])))     # the same data in the two scripts
  tru <- c(d18 = delta_sc(18, sc), R36 = R_sc(36, sc), loghr = m7$loghr_truth[m7$k == k])
  qn <- c(d18 = "delta18_first500", R36 = "R36_first500", loghr = "loghr_first500")
  ev18 <- M[, "evaluable.18"] == 1; ev36 <- M[, "evaluable.36"] == 1
  for (q in names(qn)) {
    est <- Bm[, paste0("est.", q)]; s <- Bm[, paste0("sd.", q)]; lo <- Bm[, paste0("pct_lo.", q)]; hi <- Bm[, paste0("pct_hi.", q)]
    ev <- if (q == "d18") ev18 else if (q == "R36") ev36 else rep(TRUE, nb)
    nofail <- rep(0, nb)
    mk(k, qn[q], NA, "boot", pw(est, est - Q * s, est + Q * s, s, rep(tru[q], nb), ev, nofail, tru[q]))
    mk(k, qn[q], NA, "boot_pct", pw(est, lo, hi, s, rep(tru[q], nb), ev, nofail, tru[q]))
    for (e in c("known", "estimated")) {
      if (q == "d18") v <- pw(M[, paste0(e, ".d.est.18")], M[, paste0(e, ".d.lo.18")], M[, paste0(e, ".d.hi.18")], M[, paste0(e, ".d.se.18")], rep(tru[q], nb), ev, M[, paste0(e, ".km.fail")], tru[q])
      if (q == "R36") v <- pw(M[, paste0(e, ".R.est.36")], M[, paste0(e, ".R.lo.36")], M[, paste0(e, ".R.hi.36")], M[, paste0(e, ".R.se.36")], rep(tru[q], nb), ev, M[, paste0(e, ".rmst.fail")], tru[q])
      if (q == "loghr") v <- pw(M[, paste0(e, ".cox0.b")], M[, paste0(e, ".cox0.lo")], M[, paste0(e, ".cox0.hi")], M[, paste0(e, ".cox0.se")], rep(tru[q], nb), ev, M[, paste0(e, ".cox0.fail")], tru[q])
      mk(k, qn[q], NA, e, v)
    }
  }
}

S <- do.call(rbind, rows)
write.csv(S, file.path(OUT, "m7_summary.csv"), row.names = FALSE, na = "NA")

# ---- m7_required.txt ----
line <- function(id, label, lim) {
  x <- req[[id]]
  if (is.null(x)) return(sprintf("FAIL %s: no results", label))
  nb <- sum(x$nbad); mx <- max(x$maxd)
  sprintf("%s %s: %d of %d replicates failing; largest difference %.3g%s", if (nb == 0) "OK" else "FAIL", label, nb, sum(x$nrep), mx, lim)
}
x3 <- req$r3
l3 <- if (is.null(x3)) "FAIL 3: no results" else sprintf("%s 3 (every cell; dhat, rmst_time, bhat and lr identical under \"known\" and \"estimated\"): %d of %d replicates failing; %d comparisons differ; %d cells",
                                                       if (sum(x3$nbad) == 0) "OK" else "FAIL", sum(x3$nbad), sum(x3$nrep), sum(x3$maxd), nrow(x3))
x4 <- req$r4
reqlines <- c(
  paste0(line("r1", "1 (cells 7, 8, 9, 15, 19, replicates 1-500; package \"estimated\" against the evaluation's script: sig2_surv0, sig2_surv1, sig2_dhat at the time points, band critical value and limits, RMST variance, critical value and band limits, sig_bhat_asy and |z| for FH(0,0) and FH(0,1))", " (limit 1e-8)"),
         if (!is.null(req$r1) && sum(req$r1$nbad) > 0) "; the failing replicates are in the Cox / log-rank part only, see m7_req1_diagnosis.txt" else ""),
  line("r2", "2 (cells 31, 32; \"stabilized\" against \"ipw\" from the same model, both options: sig2_dhat, band critical value, RMST variance)", " (limit 1e-8)"),
  l3,
  if (is.null(x4)) "FAIL 4: no results" else sprintf("%s 4 (cell 34, tied times; \"estimated\" FH(0,0) against the coxph residual routes, Breslow): %d of %d replicates failing; largest difference of z^2 %.3g (limit 1e-8); %s",
                                                    if (sum(x4$nbad) == 0) "OK" else "FAIL", sum(x4$nbad), sum(x4$nrep), max(x4$maxd), x4$note[1]))
writeLines(reqlines, file.path(OUT, "m7_required.txt"))
cat(reqlines, sep = "\n")
notes <- c(req$r1$note, req$r3$note[grepl(": [1-9][0-9]* comparisons not made", req$r3$note)])

# ---- m7_tables.md ----
val <- function(k, q, e, f, t = NA) {
  i <- S$cell_id == k & S$quantity == q & S$estimator == e & (if (is.na(t)) is.na(S$time) else (!is.na(S$time) & S$time == t))
  if (!any(i)) NA_real_ else S[[f]][i][1]
}
f3 <- function(x, d = 3) ifelse(is.na(x), "—", formatC(x, digits = d, format = "f"))
md <- function(head, M) c(paste0("| ", paste(head, collapse = " | "), " |"), paste0("|", paste(rep("---", length(head)), collapse = "|"), "|"),
                          apply(M, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
out <- character(0); put <- function(...) out <<- c(out, ...)
put("### M7 requirements (m7_required.txt)", "", "```", reqlines, "```", "")
if (length(notes)) put("Notes: ", paste0("- ", notes), "")
put("### Table M7 — operating characteristics with the package code (nominal 0.95; two-sided 5% log-rank test)", "",
    "Rows: cell x option (analyses A, B, C in cell 40; the two bootstrap rows in cells 38 and 39, first 500 replicates). ESS ratio = mean of ess / n by arm (control / treated).",
    "Length ratio = mean interval length under \"estimated\" over that under \"known\" (shown on the \"estimated\" row). `—` = not evaluated.", "")
rowsT <- NULL
for (k in kdone) {
  cl <- m7[m7$k == k, ]
  ests <- if (cl$sub) c("A", "B", "C") else c("known", "estimated")
  for (e in ests) {
    lr <- function(q, t = NA) if (e == "estimated") val(k, q, "estimated", "mean_length", t) / val(k, q, "known", "mean_length", t) else if (e == "B") val(k, q, "B", "mean_length", t) / val(k, q, "C", "mean_length", t) else NA
    rej <- function(g) { v <- val(k, sprintf(if (cl$strat) "logrank_strat_fh0%d" else "logrank_fh0%d", g), e, "reject"); if (cl$sc == "null") v else NA }
    rowsT <- rbind(rowsT, c(k, cl$cell, sprintf("%s%s", e, if (cl$repro) " (500)" else ""),
      sprintf("%s / %s", f3(val(k, "ess_ratio_control", e, "mean_value"), 2), f3(val(k, "ess_ratio_treated", e, "mean_value"), 2)),
      f3(val(k, "delta", e, "coverage", 18)), f3(val(k, "delta_band_3_42", e, "coverage")), f3(val(k, "R", e, "coverage", 36)),
      f3(val(k, "loghr_fh00", e, "coverage")), f3(val(k, "loghr_fh01", e, "coverage")), f3(val(k, "loghr_fh00", e, "se_sd")),
      f3(rej(0)), f3(rej(1)), f3(lr("delta", 18)), f3(lr("R", 36)), f3(lr("loghr_fh00"))))
  }
  if (boot_done && k %in% c(38, 39)) for (e in c("boot", "boot_pct")) {
    rowsT <- rbind(rowsT, c(k, cl$cell, sprintf("%s (500)", e), "—", f3(val(k, "delta18_first500", e, "coverage")), "—", f3(val(k, "R36_first500", e, "coverage")),
                            f3(val(k, "loghr_first500", e, "coverage")), "—", f3(val(k, "loghr_first500", e, "se_sd")), "—", "—", "—", "—", "—"))
  }
}
put(md(c("k", "Cell", "Option", "ESS ratio c / t", "delta(18)", "band [3, 42]", "RMST(36)", "log HR FH(0,0)", "log HR FH(0,1)", "SE/SD FH(0,0)",
         "reject FH(0,0)", "reject FH(0,1)", "length ratio delta(18)", "length ratio RMST(36)", "length ratio log HR"), rowsT), "")
if (boot_done) {
  put("### Table M7-B — bootstrap benchmark, the same first 500 replicates (coverage; SE/SD in brackets)", "")
  m <- NULL
  for (k in c(38, 39)) for (q in c("delta18_first500", "R36_first500", "loghr_first500")) {
    cp <- function(e) sprintf("%s (%s)", f3(val(k, q, e, "coverage")), f3(val(k, q, e, "se_sd"), 2))
    m <- rbind(m, c(k, m7$cell[m7$k == k], c(delta18_first500 = "delta(18)", R36_first500 = "R(36)", loghr_first500 = "log HR")[q],
                    cp("known"), cp("estimated"), cp("boot"), cp("boot_pct")))
  }
  put(md(c("k", "Cell", "Quantity", "known", "estimated", "boot", "boot_pct"), m), "")
}
put("### Correlation of the log-rank score with the survival difference at month 18 (FH(0,0), \"estimated\")", "")
m <- NULL
for (k in intersect(kdone, c(31, 32))) m <- rbind(m, c(k, m7$cell[m7$k == k], f3(val(k, "cor_lr_dhat18_estimated", "estimated", "mean_value")),
                                                        f3(val(k, "cor_lr_dhat18_empirical", "estimated", "mean_value"))))
if (!is.null(m)) put(md(c("k", "Cell", "mean estimated correlation (wlr_dhat_estimates)", "empirical correlation over replicates"), m), "")
put("### Failures and timing", "")
tm <- do.call(rbind, timing); pf <- do.call(rbind, psfail)
fl <- vapply(tm$cell_id, function(k) { x <- S[S$cell_id == k & !is.na(S$n_fail), ]; if (nrow(x)) max(x$n_fail, na.rm = TRUE) else 0 }, numeric(1))
nf <- vapply(tm$cell_id, function(k) { x <- S[S$cell_id == k & !is.na(S$n_nonfinite), ]; if (nrow(x)) max(x$n_nonfinite, na.rm = TRUE) else 0 }, numeric(1))
put(md(c("k", "Cell", "replicates", "ps_weights() errors", "largest count of failed calls (any quantity)", "largest count of non-finite intervals", "seconds", "cores", "package SHA"),
       cbind(tm$cell_id, tm$cell, tm$reps, pf$ps_fail, fl, nf, sprintf("%.0f", tm$secs), tm$cores, tm$sha)), "")
writeLines(out, file.path(OUT, "m7_tables.md"), useBytes = TRUE)
cat(sprintf("oc_summarise.R (M7): %d rows in m7_summary.csv; cells %s; bootstrap %s; tables in m7_tables.md\n",
            nrow(S), paste(kdone, collapse = ","), if (boot_done) "included" else "absent"))
