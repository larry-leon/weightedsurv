#!/usr/bin/env Rscript
# =============================================================================
# report_numbers.R -- tables B and E of REPORT_km_fixes_stage2e.md. Writes dev/audit/stage2e/report_numbers.txt.
#   B1. Riders, on the golden datasets (same datasets and calls as golden.R): the zero-weight cases before
#       (golden_s2d.rds) -> after (current code) beside the positive-weight subset; cox_results$score,
#       logrank_results$chisq and check_results() before -> after on the weighted datasets; ess by arm.
#   B2. The option, on one propensity-weighted dataset without tied times and one with: "known" -> "estimated"
#       for the log-rank z, the cox_rhogamma() SE, the SE of the survival difference at the median time and the
#       band critical value, with the coxph-route value beside each FH(0,0) number.
#   E.  All weights equal to one constant: z.score and the cox_rhogamma() SE with no weights, weights all 1 and
#       weights all 2.
#   Rscript dev/audit/stage2e/report_numbers.R [<lib with the START_SHA build>]   (for the "before" check_results() output)
# =============================================================================
if (!file.exists("DESCRIPTION")) stop("Run from the repo root.")
args <- commandArgs(trailingOnly = TRUE)
suppressMessages(pkgload::load_all(".", quiet = TRUE, export_all = FALSE, helpers = FALSE))
suppressPackageStartupMessages(library(survival))
grDevices::pdf(NULL)
S2 <- file.path("dev", "audit", "stage2e")
ex <- parse(file.path(S2, "golden.R"))
for (e in ex) if (is.call(e) && identical(e[[1]], as.name("<-")) && as.character(e[[2]]) %in% c("sim", "add_w", "add_wv", "make_data", "STRATA_OK", "ZW")) eval(e)
source(file.path("tests", "testthat", "helper-ew-ref.R"))
D <- make_data()
gold <- readRDS(file.path(S2, "golden_s2d.rds"))$values
out <- character(0); say <- function(...) out <<- c(out, sprintf(...))
q <- function(expr) { r <- NULL; invisible(utils::capture.output(r <- suppressWarnings(suppressMessages(expr)))); r }
P00 <- list(rho = 0, gamma = 0)

say("B1a. Zero weights (Rd): before (golden_s2d.rds) -> after (current code) | the positive-weight rows analysed with the same weights (no weights for w01)")
say("| dataset | weights | zero-weight rows | log-rank z | stratified z | cox_rhogamma SE FH(0,0) | SE of the survival difference at the median grid point | band critical value (draws.band = 200) | resampling SE (draws = 200) |")
say("|---|---|---|---|---|---|---|---|---|")
for (nm in ZW) for (wt in c("w01", "wz")) {
  d <- D[[nm]]; sub <- d[d[[wt]] > 0, ]; wsub <- if (wt == "w01") NULL else wt
  g <- function(p) gold[[paste0(nm, "/zw_", p)]]
  f <- q(df_counting(d, "tte", "event", "treat", weight.name = wt)); p <- q(df_counting(sub, "tte", "event", "treat", weight.name = wsub))
  cf <- q(cox_rhogamma(f, scheme = "fh", scheme_params = P00, draws = 200)); cp <- q(cox_rhogamma(p, scheme = "fh", scheme_params = P00, draws = 200))
  kf <- q(KM_diff(d, "tte", "event", "treat", weight.name = wt, draws.band = 200)); kp <- q(KM_diff(sub, "tte", "event", "treat", weight.name = wsub, at_points = kf$at_points, draws.band = 200))
  j <- ceiling(length(kf$at_points) / 2)
  strat <- "–"
  if (nm %in% STRATA_OK) {
    fs <- q(df_counting(d, "tte", "event", "treat", weight.name = wt, strata.name = "st")); ps <- q(df_counting(sub, "tte", "event", "treat", weight.name = wsub, strata.name = "st"))
    strat <- sprintf("%.6f -> %.6f | %.6f", g(paste0("dcs_", wt, "/z.score_stratified")), fs$z.score_stratified, ps$z.score_stratified)
  }
  say("| %s | %s | %d of %d | %.6f -> %.6f | %.6f | %s | %.6f -> %.6f | %.6f | %.6f -> %.6f | %.6f | %.4f -> %.4f | %.4f | %.5f -> %.5f | %.5f |",
      nm, wt, sum(d[[wt]] == 0), nrow(d), g(paste0("dc_", wt, "/z.score")), f$z.score, p$z.score, strat,
      g(paste0("cox_d0_", wt, "/fit/sig_bhat_asy")), cf$fit$sig_bhat_asy, cp$fit$sig_bhat_asy,
      sqrt(g(paste0("km_", wt, "/sig2_dhat"))[j]), sqrt(kf$sig2_dhat[j]), sqrt(kp$sig2_dhat[j]),
      g(paste0("km_", wt, "/c_alpha_band")), kf$c_alpha_band, kp$c_alpha_band,
      g(paste0("cox_d200_", wt, "/fit/sig_bhat_star")), cf$fit$sig_bhat_star, cp$fit$sig_bhat_star)
}
say("")
say("B1b. Weighted datasets (Rb): before (golden_s2d.rds) -> after. z^2 is df_counting()$z.score^2 (unchanged by Rb); coxph is coxph(weights = , robust = TRUE): $rscore with ties = \"breslow\" beside logrank chisq, summary()$robscore (default ties) beside cox score.")
say("| dataset | weights | z.score^2 | logrank_results$chisq: before -> after (coxph rscore, Breslow) | cox_results$score: before -> after (coxph robscore) | check_results() largest relative difference: before -> after | ess / n control | ess / n treated |")
say("|---|---|---|---|---|---|---|---|")
rd <- function(v) (max(v) - min(v)) / mean(c(max(v), min(v)))
for (nm in names(D)) for (wt in c("wv", if (nm == "sim_weighted") "w")) {
  d <- D[[nm]]; d$.w <- d[[wt]]
  call <- if (wt == "wv") "df_counting_wv" else "df_counting_weighted"
  g <- function(p) gold[[paste0(nm, "/", call, "/", p)]]
  dc <- q(df_counting(d, "tte", "event", "treat", weight.name = wt))
  cb <- coxph(Surv(tte, event) ~ treat, data = d, weights = .w, robust = TRUE, ties = "breslow")
  ce <- summary(coxph(Surv(tte, event) ~ treat, data = d, weights = .w, robust = TRUE))
  before3 <- c(g("lr")^2 / g("sig2_lr"), g("logrank_results/chisq"), g("z.score")^2)
  after3 <- unlist(check_results(dc, verbose = FALSE))
  say("| %s | %s | %.4f | %.4f -> %.4f (%.4f) | %.4f -> %.4f (%.4f) | %.1f%% -> %.2g%% | %.3f | %.3f |", nm, wt, dc$z.score^2,
      g("logrank_results/chisq"), dc$logrank_results$chisq, cb$rscore, g("cox_results/score"), dc$cox_results$score, ce$robscore[1],
      100 * rd(before3), 100 * rd(after3), dc$ess$ratio[2], dc$ess$ratio[1])
}
say("")
say("B1c. check_results() on sim_weighted (weights w), printed output")
d <- D$sim_weighted
if (length(args) && dir.exists(args[1])) {
  code <- sprintf('.libPaths(c("%s", .libPaths())); suppressPackageStartupMessages(library(weightedsurv)); d <- readRDS("%s"); r <- withCallingHandlers(check_results(suppressWarnings(df_counting(d, "tte", "event", "treat", weight.name = "w"))), warning = function(w) { cat("Warning:", conditionMessage(w), "\\n"); invokeRestart("muffleWarning") })',
                  normalizePath(args[1]), tf <- tempfile(fileext = ".rds"))
  saveRDS(d, tf)
  say("before (START_SHA build):"); out <- c(out, paste0("    ", system2(file.path(R.home("bin"), "Rscript"), c("-e", shQuote(code)), stdout = TRUE, stderr = TRUE)))
}
say("after (current code):")
ws <- character(0)
o <- utils::capture.output(withCallingHandlers(check_results(q(df_counting(d, "tte", "event", "treat", weight.name = "w"))),
                                               warning = function(w) { ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning") }))
out <- c(out, paste0("    ", o), if (length(ws)) paste0("    Warning: ", ws) else "    (no warning)")
dce <- q(df_counting(d, "tte", "event", "treat", weight.name = "w")); say("ess on sim_weighted (weights w):"); out <- c(out, paste0("    ", utils::capture.output(print(dce$ess, row.names = FALSE))))
say("")

say("B2. The option: \"known\" -> \"estimated\" (the coxph-route or n x m reference value under \"estimated\" in brackets). Data: n = 400, propensity model treat ~ x1 + x2 (set.seed(11), as in Gate G2 checks 6 and 7).")
say("| tied times | weights | ess / n control, treated | log-rank z FH(0,0) | stratified z FH(0,0) | cox_rhogamma SE FH(0,0) | log-rank z FH(0,1) | cox_rhogamma SE FH(0,1) | SE of the survival difference at the median time | band critical value (draws.band = 1000) | RMST difference at the last time: SE |")
say("|---|---|---|---|---|---|---|---|---|---|---|")
for (ties in c(FALSE, TRUE)) for (type in c("ipw", "stabilized")) {
  d <- ws_sim_ps(400, 11, ties = ties); pw <- ps_weights(treat ~ x1 + x2, data = d, type = type); d$w <- pw$weights
  adj <- function(r) { r <- as.numeric(r); as.numeric(r + pw$IF %*% colSums((r / d$w) * pw$dw)) }
  sres <- function(f) { r <- as.numeric(residuals(f, type = "score", weighted = TRUE)); sum(r)^2 / sum(adj(r)^2) }
  chi <- sres(coxph(Surv(time, status) ~ treat, data = d, weights = w, init = 0, iter.max = 0, robust = TRUE, ties = "breslow"))
  chis <- sres(coxph(Surv(time, status) ~ treat + strata(s), data = d, weights = w, init = 0, iter.max = 0, robust = TRUE, ties = "breslow"))
  se_b <- sqrt(sum(adj(residuals(coxph(Surv(time, status) ~ treat, data = d, weights = w, robust = TRUE, ties = "breslow"), type = "dfbeta", weighted = TRUE))^2))
  P01 <- list(rho = 0, gamma = 1)
  dfc <- function(wv, sp = P00, st = NULL) q(df_counting(d, "time", "status", "treat", weight.name = "w", strata.name = st, scheme = "fh", scheme_params = sp,
                                                      weight_variance = wv, weight_model = if (wv == "estimated") pw else NULL))
  dk <- dfc("known"); de <- dfc("estimated"); dks <- dfc("known", st = "s"); des <- dfc("estimated", st = "s"); dk1 <- dfc("known", P01); de1 <- dfc("estimated", P01)
  ck <- q(cox_rhogamma(dk, scheme = "fh", scheme_params = P00)); ce <- q(cox_rhogamma(de, scheme = "fh", scheme_params = P00))
  ck1 <- q(cox_rhogamma(dk1, scheme = "fh", scheme_params = P01)); ce1 <- q(cox_rhogamma(de1, scheme = "fh", scheme_params = P01))
  km <- function(wv) q(KM_diff(d, "time", "status", "treat", weight.name = "w", draws.band = 1000, weight_variance = wv, weight_model = if (wv == "estimated") pw else NULL))
  kk <- km("known"); ke <- km("estimated"); j <- findInterval(median(d$time), kk$at_points)
  tf <- aeqSurv(Surv(d$time, d$status)); g <- sort(unique(tf[, "time"]))
  k1 <- ws_km_ew_ref(tf[, "time"], tf[, "status"], d$w, pw$IF, pw$dw, which(d$treat == 1), g); k0 <- ws_km_ew_ref(tf[, "time"], tf[, "status"], d$w, pw$IF, pw$dw, which(d$treat == 0), g)
  vref <- ws_km_ew_diff_ref(k1, k0)[findInterval(median(d$time), g)]
  rk <- q(cumulative_rmst_bands(d, kk, "time", "status", "treat", weight.name = "w", plot = FALSE)); re <- q(cumulative_rmst_bands(d, ke, "time", "status", "treat", weight.name = "w", plot = FALSE, weight_variance = "estimated", weight_model = pw))
  n <- length(rk$rmst_time)
  say("| %s | %s | %.3f, %.3f | %.4f -> %.4f (%.4f) | %.4f -> %.4f (%.4f) | %.5f -> %.5f (%.5f) | %.4f -> %.4f | %.5f -> %.5f | %.5f -> %.5f (%.5f) | %.3f -> %.3f | %.4f -> %.4f |",
      if (ties) "yes" else "no", type, dk$ess$ratio[2], dk$ess$ratio[1], dk$z.score, de$z.score, sign(de$lr) * sqrt(chi), dks$z.score_stratified, des$z.score_stratified, sign(des$lr_stratified) * sqrt(chis),
      ck$fit$sig_bhat_asy, ce$fit$sig_bhat_asy, se_b, dk1$z.score, de1$z.score, ck1$fit$sig_bhat_asy, ce1$fit$sig_bhat_asy,
      sqrt(kk$sig2_dhat[j]), sqrt(ke$sig2_dhat[j]), sqrt(vref), kk$c_alpha_band, ke$c_alpha_band, sqrt(rk$sig2_rmst_time[n]), sqrt(re$sig2_rmst_time[n]))
  if (type == "ipw") {
    cr <- unlist(check_results(de, verbose = FALSE))
    say("  check_results() on the \"estimated\" result (%s): zlr_sq %.6f, logrank_chisq %.6f, zCox_sq %.6f; estimates identical under the two options: lr %s, bhat %s, dhat %s, rmst %s",
        if (ties) "tied times" else "no tied times", cr[1], cr[2], cr[3], identical(dk$lr, de$lr), identical(ck$fit$bhat, ce$fit$bhat), identical(kk$dhat, ke$dhat), identical(rk$rmst_time, re$rmst_time))
  }
}
say("")

say("E. All weights equal to one constant: the equal-weight formulas (gbsg, golden dataset)")
say("| weights | log-rank z FH(0,0) | sig2_lr | cox_rhogamma SE FH(0,0) | SE of S1 at the median grid point | logrank_results$method, chisq | cox_results$score | ess (control, treated) |")
say("|---|---|---|---|---|---|---|---|")
d <- D$gbsg; d$one <- 1; d$two <- 2
for (wt in list(NULL, "one", "two")) {
  dc <- q(df_counting(d, "tte", "event", "treat", weight.name = wt)); cx <- q(cox_rhogamma(dc, scheme = "fh", scheme_params = P00))
  j <- ceiling(length(dc$at_points) / 2)
  say("| %s | %.10f | %.6f | %.10f | %.10f | %s, %.6f | %.6f | %.1f, %.1f |", if (is.null(wt)) "none" else if (wt == "one") "all 1" else "all 2", dc$z.score, dc$sig2_lr, cx$fit$sig_bhat_asy,
      sqrt(dc$sig2_surv1[j]), dc$logrank_results$method, dc$logrank_results$chisq, dc$cox_results$score, dc$ess$ess[2], dc$ess$ess[1])
}
writeLines(out, file.path(S2, "report_numbers.txt"))
cat(out, sep = "\n")
