#!/usr/bin/env Rscript
# oc_summarise.R -- builds oc_summary.csv, oc_profiles.pdf and oc_tables.md (the tables of
# REPORT_km_oc.md) from the per-replicate results in dev/audit/oc/raw/ (TASK_km_oc_evaluation.md, 10-12).
#
#   Rscript dev/audit/oc/oc_summarise.R [--raw=raw] [--outdir=<dir>]
#
# oc_summary.csv is long: one row per cell x quantity x time x estimator. `time` is NA for bands, log HR
# and log-rank rows. Quantities: S0, S1, delta, R (pointwise, 11 time points); delta_band_3_42,
# delta_band_any, R_band_3_42, R_band_any; loghr_fh00, loghr_fh01; logrank_fh00, logrank_fh01 (column
# `reject`); and, for Part E (first 500 replicates of cells 5, 8, 10): delta18_first500, R36_first500,
# loghr_first500.

.here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(.here, "oc_common.R"))
opt <- parse_args(list(raw = "raw", outdir = .here))   # --raw / --outdir: for trial runs only
raw <- file.path(.here, opt$raw); OUT <- opt$outdir
Q <- qnorm(0.975); NT <- length(tt)
rd <- function(v, k) { f <- file.path(raw, sprintf("%s_cell%02d.rds", v, k)); if (file.exists(f)) readRDS(f) else NULL }
has <- function(M, nm) nm %in% colnames(M)
g11 <- function(M, prefix) M[, paste0(prefix, ".", tt), drop = FALSE]

rows <- list()
add <- function(...) rows[[length(rows) + 1]] <<- data.frame(..., stringsAsFactors = FALSE)
blank <- list(n_rep = NA, n_eval = NA, n_fail = NA, n_nonfinite = NA, truth = NA, bias = NA, bias_mcse = NA,
              emp_sd = NA, mean_se = NA, se_sd = NA, coverage = NA, mean_length = NA, reject = NA, mean_calpha = NA)
mk <- function(k, quantity, time, estimator, v) {
  b <- blank; b[names(v)] <- v
  do.call(add, c(list(cell_id = k, cell = cells$cell[k], quantity = quantity, time = time, estimator = estimator), b))
}
# pointwise metrics over replicates: evaluable replicates only; an errored call or a non-finite interval
# counts as not covering; bias, SD, SE and length average finite values
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
pw_block <- function(k, quantity, estimator, M, est, lo, hi, se, truth_nom, EV, fail, truth_mat = NULL) {
  for (i in seq_len(NT)) {
    tr <- if (is.null(truth_mat)) rep(truth_nom[i], nrow(M)) else truth_mat[, i]
    mk(k, quantity, tt[i], estimator, pw(est[, i], lo[, i], hi[, i], se[, i], tr, EV[, i] == 1, fail, truth_nom[i]))
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
cox_block <- function(k, g, estimator, M, prefix, sel = rep(TRUE, nrow(M))) {
  p <- function(f) M[sel, paste0(prefix, ".cox", g, ".", f)]
  tr <- cells$loghr_truth[k]
  v <- pw(p("b"), p("lo"), p("hi"), p("se"), rep(tr, sum(sel)), rep(TRUE, sum(sel)), p("fail"), tr)
  if (is.na(tr)) v$bias <- v$bias_mcse <- NA
  mk(k, sprintf("loghr_fh0%d", g), NA, estimator, v)
  z <- p("z")
  if (any(is.finite(z))) mk(k, sprintf("logrank_fh0%d", g), NA, estimator,
                            list(n_rep = sum(sel), n_eval = sum(is.finite(z)), n_fail = sum(p("fail") == 1, na.rm = TRUE),
                                 reject = mean(abs(z[is.finite(z)]) > 1.96)))
}

diag <- list(); chk_same <- logical(0); pkgdiff <- list(); timing <- list(); maxdiff_est <- list(); failtab <- list()
for (k in cells$cell_id) {
  rc <- rd("cran", k); rp <- rd("proposed", k)
  if (is.null(rc) || is.null(rp)) next
  C <- rc$res; P <- rp$res
  chk_same[as.character(k)] <- identical(unname(rc$chk), unname(rp$chk))
  timing[[length(timing) + 1]] <- data.frame(cell_id = k, cran_secs = rc$secs, proposed_secs = rp$secs, cores = rp$cores)
  EV <- g11(P, "evaluable")
  stopifnot(identical(unname(EV), unname(g11(C, "evaluable"))))
  S0t <- S_true(tt, 0, k); S1t <- S_true(tt, 1, k); dt_ <- S1t - S0t; Rt <- R_true(tt, k)

  ests <- list(cran = C, proposed = P)
  if (cells$ew[k]) ests$proposed_ew <- P
  if (cells$d20[k]) ests$proposed_d20 <- P
  for (e in names(ests)) {
    M <- ests[[e]]; fail <- M[, paste0(e, ".km.fail")]
    for (a in 0:1) {
      est <- g11(M, sprintf("%s.S%d.est", e, a)); se <- g11(M, sprintf("%s.S%d.se", e, a))
      pw_block(k, paste0("S", a), e, M, est, est - Q * se, est + Q * se, se, if (a == 0) S0t else S1t, EV, fail)
    }
    pw_block(k, "delta", e, M, g11(M, paste0(e, ".d.est")), g11(M, paste0(e, ".d.lo")), g11(M, paste0(e, ".d.hi")),
             g11(M, paste0(e, ".d.se")), dt_, EV, fail)
    band_block(k, "delta", e, M, paste0(e, ".dband"))
  }
  maxdiff_est[[length(maxdiff_est) + 1]] <- data.frame(cell_id = k,
    maxdiff_delta = max(abs(g11(C, "cran.d.est") - g11(P, "proposed.d.est")), na.rm = TRUE),
    maxdiff_R = if (cells$do_rmst[k]) max(abs(g11(C, "cran.R.est") - g11(P, "proposed.R.est")), na.rm = TRUE) else NA)

  if (cells$do_rmst[k]) {
    ests <- list(cran = C, proposed = P, proposed_rect = P)
    if (cells$ew[k]) ests$proposed_ew <- P
    for (e in names(ests)) {
      M <- ests[[e]]
      pw_block(k, "R", e, M, g11(M, paste0(e, ".R.est")), g11(M, paste0(e, ".R.lo")), g11(M, paste0(e, ".R.hi")),
               g11(M, paste0(e, ".R.se")), Rt, EV, M[, paste0(e, ".rmst.fail")], truth_mat = g11(M, "R.truth"))
      band_block(k, "R", e, M, paste0(e, ".rband"))
    }
  }
  if (cells$do_cox[k]) {
    for (g in 0:1) {
      cox_block(k, g, "cran", C, "cran")
      if (g == 0) {
        cox_block(k, g, "cran_mean1", C, "cran_mean1")
        cox_block(k, g, "cran_star", C, "cran_star", sel = !is.na(C[, "cran_star.cox0.fail"]))
      }
      cox_block(k, g, "proposed", P, "proposed")
      if (cells$ew[k]) cox_block(k, g, "proposed_ew", P, "proposed_ew")
      f <- c("b", "se", "lo", "hi", "z")
      pkgdiff[[length(pkgdiff) + 1]] <- data.frame(cell_id = k, g = g,
        maxdiff = max(abs(C[, paste0("cran.cox", g, ".", f)] - P[, paste0("pkg.cox", g, ".", f)]), na.rm = TRUE),
        rootcheck = max(P[, paste0("proposed.cox", g, ".rootcheck")], na.rm = TRUE))
    }
  }
  # errored calls; the proposed_ew Cox "failures" are replicates that did not pass the alignment
  # confirmation of section 7 (cx$fit$time against df$time[order(df$time)]) and are counted apart
  fc <- function(M, ew) { cc <- grep("\\.fail$", colnames(M), value = TRUE); cc <- cc[grepl("^proposed_ew\\.cox", cc) == ew]
                          if (length(cc)) sum(M[, cc] == 1, na.rm = TRUE) else 0 }
  rcv <- if (cells$do_cox[k]) pmax(P[, "proposed.cox0.rootcheck"], P[, "proposed.cox1.rootcheck"], na.rm = TRUE) else 0
  failtab[[length(failtab) + 1]] <- data.frame(cell_id = k, cran = fc(C, FALSE), proposed = fc(P, FALSE),
                                               ew_align = fc(P, TRUE) / 2, root_gt = sum(rcv > 1e-6))
  diag[[length(diag) + 1]] <- data.frame(cell_id = k, cell = cells$cell[k],
    ess0 = mean(P[, "diag.ess0"]), ess1 = mean(P[, "diag.ess1"]),
    share0 = mean(P[, "diag.share180"], na.rm = TRUE), share1 = mean(P[, "diag.share181"], na.rm = TRUE),
    share0_max = max(P[, "diag.share180"], na.rm = TRUE), share1_max = max(P[, "diag.share181"], na.rm = TRUE),
    eval42 = sum(EV[, "evaluable.42"]), eval48 = sum(EV[, "evaluable.48"]))
}

# ---- Part E ----
boot_done <- all(file.exists(file.path(raw, sprintf("boot_cell%02d.rds", c(5, 8, 10)))))
if (boot_done) for (k in c(5, 8, 10)) {
  bt <- readRDS(file.path(raw, sprintf("boot_cell%02d.rds", k))); Bm <- bt$res; nb <- nrow(Bm)
  P <- rd("proposed", k)$res[seq_len(nb), , drop = FALSE]
  stopifnot(identical(unname(bt$chk), unname(rd("proposed", k)$chk[seq_len(nb)])))
  tru <- c(d18 = delta_true(18, k), R36 = R_true(36, k), loghr = cells$loghr_truth[k])
  qn <- c(d18 = "delta18_first500", R36 = "R36_first500", loghr = "loghr_first500")
  ev <- rep(TRUE, nb); nofail <- rep(0, nb)
  for (q in names(qn)) {
    est <- Bm[, paste0("est.", q)]; s <- Bm[, paste0("sd.", q)]
    mk(k, qn[q], NA, "boot", pw(est, est - Q * s, est + Q * s, s, rep(tru[q], nb), ev, nofail, tru[q]))
    lo <- Bm[, paste0("pct_lo.", q)]; hi <- Bm[, paste0("pct_hi.", q)]
    mk(k, qn[q], NA, "boot_pct", pw(est, lo, hi, s, rep(tru[q], nb), ev, nofail, tru[q]))
  }
  for (e in c("proposed", if (cells$ew[k]) "proposed_ew")) {
    mk(k, qn["d18"], NA, e, pw(P[, paste0(e, ".d.est.18")], P[, paste0(e, ".d.lo.18")], P[, paste0(e, ".d.hi.18")],
                              P[, paste0(e, ".d.se.18")], rep(tru["d18"], nb), ev, P[, paste0(e, ".km.fail")], tru["d18"]))
    mk(k, qn["loghr"], NA, e, pw(P[, paste0(e, ".cox0.b")], P[, paste0(e, ".cox0.lo")], P[, paste0(e, ".cox0.hi")],
                                P[, paste0(e, ".cox0.se")], rep(tru["loghr"], nb), ev, P[, paste0(e, ".cox0.fail")], tru["loghr"]))
  }
  for (e in c("proposed_rect", if (cells$ew[k]) "proposed_ew")) {
    mk(k, qn["R36"], NA, e, pw(P[, paste0(e, ".R.est.36")], P[, paste0(e, ".R.lo.36")], P[, paste0(e, ".R.hi.36")],
                              P[, paste0(e, ".R.se.36")], P[, "R.truth.36"], ev, P[, paste0(e, ".rmst.fail")], tru["R36"]))
  }
}

S <- do.call(rbind, rows)
S <- S[order(S$cell_id), ]
num <- vapply(S, is.numeric, logical(1))
write.csv(S, file.path(OUT, "oc_summary.csv"), row.names = FALSE, na = "NA")
diag <- do.call(rbind, diag); pkgdiff <- do.call(rbind, pkgdiff); timing <- do.call(rbind, timing)
maxdiff_est <- do.call(rbind, maxdiff_est); failtab <- do.call(rbind, failtab)

# ---- lookups ----
val <- function(k, quantity, estimator, col, time = NA) {
  s <- S$cell_id == k & S$quantity == quantity & S$estimator == estimator & (if (is.na(time)) is.na(S$time) else (!is.na(S$time) & S$time == time))
  if (any(s)) S[[col]][s][1] else NA
}
prof <- function(k, quantity, estimator, col) vapply(tt, function(t) val(k, quantity, estimator, col, t), numeric(1))
f3 <- function(x) ifelse(is.na(x), "—", sprintf("%.3f", x))
f4 <- function(x) ifelse(is.na(x), "—", sprintf("%.4f", x))
f2 <- function(x) ifelse(is.na(x), "—", sprintf("%.2f", x))
md <- function(header, mat) {
  c(paste0("| ", paste(header, collapse = " | "), " |"), paste0("|", paste(rep("---", length(header)), collapse = "|"), "|"),
    apply(mat, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
}
kdone <- sort(unique(S$cell_id))
out <- character(0); put <- function(...) out <<- c(out, ...)
tcols <- paste0("t=", tt)

# ---- Table H ----
put("### Table H — headline coverage (nominal 0.95)", "",
    "Each entry: `cran` / `proposed` / `proposed_ew`. `—` = not evaluated (`proposed_ew` exists only in the E cells).", "")
tri <- function(k, q, time = NA, ests = c("cran", "proposed", "proposed_ew"))
  paste(f3(vapply(ests, function(e) val(k, q, e, "coverage", time), numeric(1))), collapse = " / ")
put(md(c("k", "Cell", "delta(18)", "band for delta, [3, 42]", "R(36)", "log HR, FH(0,0)"),
       t(sapply(kdone, function(k) c(k, cells$cell[k], tri(k, "delta", 18), tri(k, "delta_band_3_42"), tri(k, "R", 36), tri(k, "loghr_fh00"))))), "")

# ---- profile tables ----
profile_table <- function(title, quantity, col, ests, fmt = f3, band = TRUE, bandcol = col, note = NULL) {
  put(paste0("### ", title), "", note, if (!is.null(note)) "")
  m <- NULL
  for (k in kdone) for (e in ests) {
    p <- prof(k, quantity, e, col)
    if (all(is.na(p))) next
    r <- c(k, cells$cell[k], e, fmt(p))
    if (band) r <- c(r, fmt(val(k, paste0(quantity, "_band_3_42"), e, bandcol)),
                     fmt(if (bandcol == "coverage") val(k, paste0(quantity, "_band_any"), e, bandcol) else NA))
    m <- rbind(m, r)
  }
  put(md(c("k", "Cell", "Estimator", tcols, if (band) c("band [3, 42]", "band anywhere")), m), "")
}
E5 <- c("cran", "proposed", "proposed_ew", "proposed_d20")
profile_table("Table P1 — coverage of delta(t)", "delta", "coverage", E5)
profile_table("Table P2 — mean interval length of delta(t) (band: mean width on [3, 42])", "delta", "mean_length", E5, bandcol = "mean_length")
put("### Table P3 — bias of the delta(t) estimate", "",
    "One row of bias per cell (`proposed`; the `cran` estimate is the same up to the last column), then the truth. Last column: the largest absolute difference between the two versions' estimates over all replicates and time points.", "")
m <- NULL
for (k in kdone) {
  m <- rbind(m, c(k, cells$cell[k], "bias", f4(prof(k, "delta", "proposed", "bias")), sprintf("%.2g", maxdiff_est$maxdiff_delta[maxdiff_est$cell_id == k])),
             c("", "", "truth", f4(prof(k, "delta", "proposed", "truth")), ""))
}
put(md(c("k", "Cell", "", tcols, "max abs(cran - proposed)"), m), "")
mcse_d <- range(S$bias_mcse[S$quantity == "delta" & S$estimator == "proposed"], na.rm = TRUE)
put(sprintf("Monte Carlo SE of the bias of delta(t): %.4f to %.4f.", mcse_d[1], mcse_d[2]), "")
R4 <- c("cran", "proposed", "proposed_rect", "proposed_ew")
profile_table("Table P4 — coverage of R(t)", "R", "coverage", R4)
profile_table("Table P5 — mean interval length of R(t) (band: mean width on [3, 42])", "R", "mean_length", R4, bandcol = "mean_length")
put("### Table P6 — bias of R(t)", "",
    "Bias = mean of (estimate at the last grid point t_j <= t) minus R(t_j). Truth row: R(t).", "")
m <- NULL
for (k in kdone) if (cells$do_rmst[k]) {
  for (e in c("cran", "proposed", "proposed_rect")) m <- rbind(m, c(k, cells$cell[k], e, f4(prof(k, "R", e, "bias"))))
  m <- rbind(m, c("", "", "truth", f4(prof(k, "R", "proposed", "truth"))))
}
put(md(c("k", "Cell", "Estimator", tcols), m), "")
mcse_r <- range(S$bias_mcse[S$quantity == "R" & S$estimator == "proposed"], na.rm = TRUE)
put(sprintf("Monte Carlo SE of the bias of R(t): %.4f to %.4f (largest at month 48).", mcse_r[1], mcse_r[2]), "")

# ---- Table C ----
trip <- function(k, q, e) {
  v <- c(val(k, q, e, "se_sd"), val(k, q, e, "coverage"), val(k, q, e, "mean_length"))
  if (all(is.na(v))) "—" else paste(f3(v), collapse = " · ")
}
biasm <- function(k, q, e) { b <- val(k, q, e, "bias"); if (is.na(b)) "—" else sprintf("%.4f (%.4f)", b, val(k, q, e, "bias_mcse")) }
put("### Table C — log HR", "",
    "Each estimator entry: SE/SD · coverage · mean interval length. Bias (Monte Carlo SE) is that of the package estimate `bhat`, which `cran`, `cran_mean1`, `proposed` and `proposed_ew` share; `cran_star` is centred at the package's de-biased estimate and uses the first 500 replicates. Coverage `—` where there is no single true hazard ratio (cells 3, 17, 18, 19); cell 12 not evaluated.", "",
    "**FH(0,0)**", "")
kc <- kdone[cells$do_cox[kdone]]
put(md(c("k", "Cell", "bias (MC SE)", "cran", "cran_mean1", "cran_star (500)", "proposed", "proposed_ew"),
       t(sapply(kc, function(k) c(k, cells$cell[k], biasm(k, "loghr_fh00", "cran"),
                                  sapply(c("cran", "cran_mean1", "cran_star", "proposed", "proposed_ew"), function(e) trip(k, "loghr_fh00", e)))))), "",
    "**FH(0,1)**", "")
put(md(c("k", "Cell", "bias (MC SE)", "cran", "proposed", "proposed_ew"),
       t(sapply(kc, function(k) c(k, cells$cell[k], biasm(k, "loghr_fh01", "cran"),
                                  sapply(c("cran", "proposed", "proposed_ew"), function(e) trip(k, "loghr_fh01", e)))))), "")

# ---- Table N ----
put("### Table N — null cells: rejection rate of the two-sided 5% log-rank test", "")
kn <- intersect(13:16, kdone)
put(md(c("k", "Cell", "FH(0,0) cran", "FH(0,0) cran_mean1", "FH(0,0) proposed", "FH(0,0) proposed_ew", "FH(0,1) cran", "FH(0,1) proposed", "FH(0,1) proposed_ew"),
       t(sapply(kn, function(k) c(k, cells$cell[k],
                                  f3(sapply(c("cran", "cran_mean1", "proposed", "proposed_ew"), function(e) val(k, "logrank_fh00", e, "reject"))),
                                  f3(sapply(c("cran", "proposed", "proposed_ew"), function(e) val(k, "logrank_fh01", e, "reject"))))))), "")

# ---- Table B ----
if (boot_done) {
  put("### Table B — bootstrap benchmark (first 500 replicates, 200 resamples)", "",
      "Each entry: coverage · mean interval length. `proposed` for R(36) is `proposed_rect`.", "")
  m <- NULL
  for (k in c(5, 8, 10)) for (q in c("delta18_first500", "R36_first500", "loghr_first500")) {
    cp <- function(e) { v <- c(val(k, q, e, "coverage"), val(k, q, e, "mean_length")); if (all(is.na(v))) "—" else paste(f3(v), collapse = " · ") }
    pe <- if (q == "R36_first500") "proposed_rect" else "proposed"
    m <- rbind(m, c(k, cells$cell[k], c(delta18_first500 = "delta(18)", R36_first500 = "R(36)", loghr_first500 = "log HR")[q],
                    cp(pe), cp("proposed_ew"), cp("boot"), cp("boot_pct")))
  }
  put(md(c("k", "Cell", "Quantity", "proposed", "proposed_ew", "boot", "boot_pct"), m), "")
}

# ---- Table D ----
put("### Table D — diagnostics and checks", "",
    "ESS ratio = (sum w)^2 / (n_a sum w^2), mean over replicates. Influence share = largest single share of the influence sum at month 18, mean (maximum) over replicates. Non-finite: pointwise intervals of S_0, S_1, delta, R that are non-finite in an evaluable replicate (count over replicates x time points) + band replicates (delta and RMST) with a non-finite critical value or limit. c_alpha: mean critical value of the band for delta.", "")
nonfin <- function(k, e) {
  s <- S$cell_id == k & S$estimator == e
  paste0(sum(S$n_nonfinite[s & S$quantity %in% c("S0", "S1", "delta", "R")]), " + ",
         sum(S$n_nonfinite[s & S$quantity %in% c("delta_band_3_42", "R_band_3_42")]))
}
put(md(c("k", "Cell", "ESS ratio, control", "ESS ratio, treated", "influence share, control", "influence share, treated",
         "evaluable at 42", "evaluable at 48", "non-finite, cran", "non-finite, proposed", "c_alpha, cran", "c_alpha, proposed"),
       t(sapply(kdone, function(k) { d <- diag[diag$cell_id == k, ]
         c(k, cells$cell[k], f3(d$ess0), f3(d$ess1), sprintf("%.3f (%.3f)", d$share0, d$share0_max), sprintf("%.3f (%.3f)", d$share1, d$share1_max),
           d$eval42, d$eval48, nonfin(k, "cran"), nonfin(k, "proposed"),
           f3(val(k, "delta_band_3_42", "cran", "mean_calpha")), f3(val(k, "delta_band_3_42", "proposed", "mean_calpha"))) }))), "")
put(sprintf("- Data checksums identical between the two processes in every replicate: %s (%d cells x 2,000 replicates).",
            if (all(chk_same)) "yes" else paste("NO, cells", paste(names(chk_same)[!chk_same], collapse = ", ")), length(chk_same)),
    sprintf("- Package Cox output (bhat, sig_bhat_asy, log limits of hr_ci_asy, df_counting z.score), HEAD against `cran`: largest absolute difference over all cells, replicates and both schemes %.3g.", max(pkgdiff$maxdiff)),
    sprintf("- `cox_robust()` root check, abs(U) / sqrt(sum psi^2) at bhat: above 1e-6 in %d of %d replicates (largest %.3g); at most %.3g in the others. In those replicates the package's time fix (`survival::aeqSurv`) merged two nearly equal times into a tie, which `cox_robust()` does not handle (cells: %s).",
            sum(failtab$root_gt), 2000 * sum(cells$do_cox), max(pkgdiff$rootcheck), 2e-10, paste(sprintf("%d: %d", failtab$cell_id[failtab$root_gt > 0], failtab$root_gt[failtab$root_gt > 0]), collapse = ", ")),
    sprintf("- Calls that errored: cran %d, proposed %d (all cells, all package and script calls).", sum(failtab$cran), sum(failtab$proposed)),
    sprintf("- `proposed_ew`, log HR and log-rank: replicates that did not pass the alignment confirmation (`cx$fit$time` equal to `df$time[order(df$time)]`): %s. They are the replicates with merged times in the E cells; they are counted as failures and as not covering, for both schemes.",
            paste(sprintf("cell %d: %d", failtab$cell_id[cells$ew[failtab$cell_id]], failtab$ew_align[cells$ew[failtab$cell_id]]), collapse = ", ")),
    "- Non-finite pointwise intervals: one replicate of cell 7 at month 1 (no grid point at or before month 1), every estimator.",
    sprintf("- Wall-clock (sum over cells, %d cores per version, the two versions running at the same time): cran %.1f min, proposed %.1f min.",
            timing$cores[1], sum(timing$cran_secs) / 60, sum(timing$proposed_secs) / 60), "")

# ---- Appendix tables ----
put("### Appendix A1 — S_0(t) and S_1(t): coverage, mean interval length, SE/SD", "")
for (a in 0:1) for (cc in list(c("coverage", "coverage"), c("mean_length", "mean interval length"), c("se_sd", "SE/SD"))) {
  profile_table(sprintf("A1: S_%d(t), %s", a, cc[2]), paste0("S", a), cc[1], c("cran", "proposed", "proposed_ew"), band = FALSE)
}
profile_table("Appendix A2 — SE/SD of delta(t)", "delta", "se_sd", E5, band = FALSE)
profile_table("Appendix A3 — SE/SD of R(t)", "R", "se_sd", R4, band = FALSE)
put("### Appendix A4 — evaluable replicates by time point (of 2,000)", "")
put(md(c("k", "Cell", tcols), t(sapply(kdone, function(k) c(k, cells$cell[k], prof(k, "delta", "proposed", "n_eval"))))), "")
writeLines(out, file.path(OUT, "oc_tables.md"), useBytes = TRUE)

# ---- figures (base graphics) ----
mcsd <- sqrt(0.95 * 0.05 / 2000)
cols <- c(cran = "#B4502D", proposed = "#1F5FA8", proposed_ew = "#2E8B57", proposed_d20 = "#8E5BB5", proposed_rect = "#C89A1B", boot = "#444444")
ltys <- c(cran = 1, proposed = 1, proposed_ew = 1, proposed_d20 = 2, proposed_rect = 2, boot = 1)
pchs <- c(cran = 1, proposed = 16, proposed_ew = 15, proposed_d20 = 2, proposed_rect = 4, boot = 17)
panel_page <- function(quantity, col, ests, main, ylab, ref = NULL, rel = FALSE) {
  par(mfrow = c(4, 5), mar = c(2.6, 2.8, 1.8, 0.6), oma = c(1.5, 1.2, 2.4, 0), mgp = c(1.6, 0.5, 0), cex = 0.62, tcl = -0.3)
  getp <- function(k, e) { p <- prof(k, quantity, e, col); if (rel) p / prof(k, quantity, "proposed", col) else p }
  allv <- unlist(lapply(kdone, function(k) lapply(ests, function(e) getp(k, e))))
  allv <- allv[is.finite(allv)]
  yl <- if (rel) range(c(allv, 1)) else c(max(0.5, floor(min(allv) * 20) / 20), 1)
  for (k in cells$cell_id) {
    plot(range(tt), yl, type = "n", xlab = "month", ylab = ylab, main = sprintf("%d  %s", k, cells$cell[k]), xaxt = "n")
    axis(1, at = c(1, 6, 12, 18, 24, 30, 36, 42, 48))
    if (!is.null(ref)) { abline(h = ref, col = "grey40"); abline(h = ref + c(-2, 2) * mcsd, lty = 2, col = "grey55") } else abline(h = 1, col = "grey40")
    any_drawn <- FALSE
    for (e in ests) {
      p <- getp(k, e); if (all(is.na(p))) next
      if (rel && e == "proposed") next
      lines(tt, pmax(p, yl[1]), col = cols[e], lty = ltys[e], lwd = 1.3); points(tt, pmax(p, yl[1]), col = cols[e], pch = pchs[e], cex = 0.7)
      any_drawn <- TRUE
    }
    if (!any_drawn && !(rel)) text(mean(range(tt)), mean(yl), "not evaluated", col = "grey40")
  }
  plot.new()
  le <- if (rel) setdiff(ests, "proposed") else ests
  legend("center", legend = le, col = cols[le], lty = ltys[le], pch = pchs[le], lwd = 1.3, bty = "n", cex = 1.1)
  mtext(main, outer = TRUE, cex = 0.95, font = 2, line = 0.6)
}
pdf(file.path(OUT, "oc_profiles.pdf"), width = 13, height = 9.5)
panel_page("delta", "coverage", E5, sprintf("Page 1. Coverage of the pointwise 95%% interval for delta(t). Line at 0.95; dashed: 0.95 +/- 2 Monte Carlo SD (%.4f). Values below the axis minimum are drawn at the minimum.", mcsd),
           "coverage", ref = 0.95)
panel_page("R", "coverage", R4, "Page 2. Coverage of the pointwise 95% interval for R(t). Line at 0.95; dashed: 0.95 +/- 2 Monte Carlo SD. Values below the axis minimum are drawn at the minimum.",
           "coverage", ref = 0.95)
panel_page("delta", "mean_length", E5, "Page 3. Mean interval length of delta(t), each estimator relative to proposed (= 1).", "length / proposed", rel = TRUE)
par(mfrow = c(3, 1), mar = c(5.2, 4, 2.2, 1), oma = c(0, 0, 2.2, 0), cex = 0.8, mgp = c(2.2, 0.7, 0))
dot_panel <- function(q, title, bootq = NULL) {
  ests <- c("cran", "proposed", "proposed_ew")
  v <- sapply(ests, function(e) sapply(cells$cell_id, function(k) val(k, q, e, "coverage")))
  bv <- if (!is.null(bootq) && boot_done) sapply(cells$cell_id, function(k) val(k, bootq, "boot", "coverage")) else rep(NA, 19)
  yl <- c(max(0.4, floor(min(c(v, bv), na.rm = TRUE) * 20) / 20), 1)
  plot(c(1, 19), yl, type = "n", xaxt = "n", xlab = "", ylab = "coverage", main = title)
  axis(1, at = 1:19, labels = paste(1:19, cells$cell), las = 2, cex.axis = 0.75)
  abline(h = 0.95, col = "grey40"); abline(h = 0.95 + c(-2, 2) * mcsd, lty = 2, col = "grey55"); abline(v = 1:19, col = "grey92")
  off <- c(cran = -0.18, proposed = 0, proposed_ew = 0.18)
  for (e in ests) points(1:19 + off[e], pmax(v[, e], yl[1]), pch = pchs[e], col = cols[e], cex = 1.2)
  points(1:19 + 0.34, pmax(bv, yl[1]), pch = 17, col = cols["boot"], cex = 1.1)
  legend("bottomleft", legend = c("cran", "proposed", "proposed_ew", if (!is.null(bootq) && boot_done) "boot (first 500 replicates)"),
         pch = c(1, 16, 15, if (!is.null(bootq) && boot_done) 17), col = c(cols[ests], if (!is.null(bootq) && boot_done) cols["boot"]), bty = "n", horiz = TRUE, cex = 0.9)
}
dot_panel("delta_band_3_42", "Simultaneous band for delta, coverage on [3, 42]")
dot_panel("R_band_3_42", "Simultaneous band for R (RMST difference), coverage on [3, 42]")
dot_panel("loghr_fh00", "Log hazard ratio, FH(0,0): coverage of the 95% interval", bootq = "loghr_first500")
mtext("Page 4. Coverage by cell. Open circles cran, filled proposed, squares proposed_ew, triangles boot. Line at 0.95; dashed: 0.95 +/- 2 Monte Carlo SD (2,000 replicates).", outer = TRUE, cex = 0.8, font = 2, line = 0.5)
invisible(dev.off())
cat(sprintf("oc_summarise.R: %d rows in oc_summary.csv; cells %s; Part E %s; tables in oc_tables.md; figures in oc_profiles.pdf\n",
            nrow(S), paste(range(kdone), collapse = "-"), if (boot_done) "included" else "absent"))
