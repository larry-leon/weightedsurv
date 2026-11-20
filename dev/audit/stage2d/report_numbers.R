#!/usr/bin/env Rscript
# =============================================================================
# report_numbers.R -- tables B and E of REPORT_km_fixes_stage2d.md
#   B. Per golden dataset (weights wv): log-rank z for FH(0,0) and FH(0,1), the cox_rhogamma() SE for FH(0,0)
#      and FH(0,1), the stratified z (datasets with strata) and the resampling SE (draws = 200, FH(0,0)),
#      Stage 2c (golden_s2c.rds) -> Stage 2d (current code), with the coxph(weights = wv, robust = TRUE,
#      ties = "breslow") value beside the FH(0,0) ones. Same datasets and calls as golden.R.
#   E. On one weighted dataset: the outputs that are still computed without the weights or compared with an
#      unweighted reference.
# Writes dev/audit/stage2d/report_numbers.txt.
# =============================================================================
if (!file.exists("DESCRIPTION")) stop("Run from the repo root.")
suppressMessages(pkgload::load_all(".", quiet = TRUE, export_all = FALSE, helpers = FALSE))
suppressPackageStartupMessages(library(survival))
grDevices::pdf(NULL)
S2 <- file.path("dev", "audit", "stage2d")
ex <- parse(file.path(S2, "golden.R"))
for (e in ex) if (is.call(e) && identical(e[[1]], as.name("<-")) && as.character(e[[2]]) %in% c("sim", "add_w", "add_wv", "make_data", "STRATA_OK")) eval(e)
D <- make_data()
gold <- readRDS(file.path(S2, "golden_s2c.rds"))$values
out <- character(0); say <- function(...) out <<- c(out, sprintf(...))
q <- function(expr) suppressWarnings(suppressMessages(expr))
P00 <- list(rho = 0, gamma = 0); P01 <- list(rho = 0, gamma = 1)
say("B. Golden datasets, weights wv: Stage 2c (golden_s2c.rds) -> Stage 2d (current code)")
say("| dataset | n | log-rank z FH(0,0): 2c -> 2d (coxph robust) | log-rank z FH(0,1): 2c -> 2d | cox_rhogamma SE FH(0,0): 2c -> 2d (coxph robust) | SE FH(0,1): 2c -> 2d | stratified z FH(0,0): 2c -> 2d (coxph robust) | resampling SE FH(0,0), draws = 200: 2c -> 2d |")
say("|---|---|---|---|---|---|---|---|")
unch <- character(0)
for (nm in names(D)) {
  d <- D[[nm]]
  cp <- coxph(Surv(tte, event) ~ treat, data = d, weights = wv, robust = TRUE, ties = "breslow")
  dc0 <- q(df_counting(d, "tte", "event", "treat", weight.name = "wv"))
  dc1 <- q(df_counting(d, "tte", "event", "treat", weight.name = "wv", scheme = "fh", scheme_params = P01))
  c0 <- q(cox_rhogamma(dc0, scheme = "fh", scheme_params = P00, draws = 200))
  c1 <- q(cox_rhogamma(dc0, scheme = "fh", scheme_params = P01, draws = 0))
  g <- function(p) gold[[paste0(nm, "/", p)]]
  strat <- "–"
  if (nm %in% STRATA_OK) {
    dcs <- q(df_counting(d, "tte", "event", "treat", weight.name = "wv", strata.name = "st"))
    cps <- coxph(Surv(tte, event) ~ treat + strata(st), data = d, weights = wv, robust = TRUE, ties = "breslow")
    strat <- sprintf("%.3f -> %.3f (%.3f)", g("df_counting_strata_wv/z.score_stratified"), dcs$z.score_stratified, sign(dcs$lr_stratified) * sqrt(cps$rscore))
    unch <- c(unch, sprintf("%s stratified score %.10g -> %.10g", nm, g("df_counting_strata_wv/lr_stratified"), dcs$lr_stratified))
  }
  say("| %s | %d | %.3f -> %.3f (%.3f) | %.3f -> %.3f | %.4f -> %.4f (%.4f) | %.4f -> %.4f | %s | %.4f -> %.4f |", nm, nrow(d),
      g("df_counting_wv/z.score"), dc0$z.score, sign(dc0$lr) * sqrt(cp$rscore), g("df_counting_fh01_wv/z.score"), dc1$z.score,
      g("cox_fh00_d0_wv/fit/sig_bhat_asy"), c0$fit$sig_bhat_asy, sqrt(cp$var[1, 1]), g("cox_fh01_d0_wv/fit/sig_bhat_asy"), c1$fit$sig_bhat_asy,
      strat, g("cox_fh00_d200_wv/fit/sig_bhat_star"), c0$fit$sig_bhat_star)
  unch <- c(unch, sprintf("%s: lr %.10g -> %.10g; bhat FH(0,0) %.10g -> %.10g; bhat FH(0,1) %.10g -> %.10g", nm, g("df_counting_wv/lr"), dc0$lr,
                          g("cox_fh00_d0_wv/fit/bhat"), c0$fit$bhat, g("cox_fh01_d0_wv/fit/bhat"), c1$fit$bhat))
}
say(""); say("Point estimates (unchanged): score and bhat, 2c -> 2d"); out <- c(out, paste0("  ", unch))

# ---- E: outputs still computed without the weights, or compared with an unweighted reference ----
d <- D$sim_weighted
dc <- q(df_counting(d, "tte", "event", "treat", weight.name = "wv"))
cp <- coxph(Surv(tte, event) ~ treat, data = d, weights = wv, robust = TRUE, ties = "breslow")
sd0 <- survdiff(Surv(tte, event) ~ treat, data = d)
say(""); say("E. dataset sim_weighted (n = %d), weights wv", nrow(d))
say("df_counting()$z.score^2 = %.4f (coxph robust score test %.4f)", dc$z.score^2, cp$rscore)
say("df_counting()$logrank_results$chisq = %.4f; survdiff() without weights = %.4f; logrank_text: %s; zlogrank_text: %s", dc$logrank_results$chisq, sd0$chisq,
    dc$logrank_results$logrank_text, dc$zlogrank_text)
cr <- NULL; wmsg <- character(0)
invisible(utils::capture.output(cr <- withCallingHandlers(check_results(dc, verbose = TRUE), warning = function(w) { wmsg <<- c(wmsg, conditionMessage(w)); invokeRestart("muffleWarning") })))
say("check_results(): zlr_sq = %.4f, logrank_chisq = %.4f, zCox_sq = %.4f; warning: %s", cr$zlr_sq, cr$logrank_chisq, cr$zCox_sq, if (length(wmsg)) wmsg else "(none)")
sm <- summary(dc$cox_results$cox_fit)
say("df_counting()$cox_results$score = %.4f = summary(coxph)$sctest (model-based score test of the weighted fit; Efron ties); summary(coxph)$robscore = %.4f; coxph(ties = 'breslow')$rscore = %.4f",
    dc$cox_results$score, sm$robscore[1], cp$rscore)
say("df_counting()$cox_results: coxph(weights, robust = TRUE) with the default ties (efron): HR %.4f, robust SE %.4f; cox_rhogamma FH(0,0): HR %.4f, SE %.4f (Breslow)",
    dc$cox_results$hr, sqrt(dc$cox_results$cox_fit$var[1, 1]), exp(q(cox_rhogamma(dc, scheme = "fh", scheme_params = P00))$fit$bhat),
    q(cox_rhogamma(dc, scheme = "fh", scheme_params = P00))$fit$sig_bhat_asy)
pk <- plotKM.band_subgroups(d, "tte", "event", "treat", weight.name = "wv", draws.band = 0)
rp <- pk$rpoints
say("plotKM.band_subgroups() risk table at t = %s: printed counts %s (unweighted numbers of subjects with time >= t); weighted numbers at risk %s",
    paste(rp[1:4], collapse = "/"), paste(colSums(outer(d$tte, rp[1:4], ">=")), collapse = "/"),
    paste(round(colSums(outer(d$tte, rp[1:4], ">=") * d$wv), 1), collapse = "/"))
say("plotKM.band_subgroups() max_tau = %.4f = smaller of the arms' last event times, weights not used (a zero-weight event would count)", pk$fit_itt$max_tau)
say("formals: plotKM.band_subgroups draws = %s, draws.band = %s; KM_diff draws = %s, draws.band = %s", deparse(formals(plotKM.band_subgroups)$draws),
    deparse(formals(plotKM.band_subgroups)$draws.band), deparse(formals(KM_diff)$draws), deparse(formals(KM_diff)$draws.band))
writeLines(out, file.path(S2, "report_numbers.txt")); cat(out, sep = "\n")
