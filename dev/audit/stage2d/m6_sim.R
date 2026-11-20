#!/usr/bin/env Rscript
# m6_sim.R -- M6 of TASK_km_fixes_stage2d.md (section 6): Monte Carlo of the subject-weighted log-rank test and
# Cox estimate with the package code (HEAD), 11 cells x 2,000 replicates, N = 400, FH(0,0) and FH(0,1).
#
#   Rscript dev/audit/stage2d/m6_sim.R [reps = 2000] [cores = 12]
#
# Per replicate and scheme: df_counting() and cox_rhogamma(draws = 0) as in the evaluation (Appendix A5 of
# REPORT_km_oc.md); first 500 replicates, FH(0,0): cox_rhogamma(draws = 200) after set.seed(3e7 + 1e6 * k + r);
# cells 5 and 14: wlr_dhat_estimates(dc, tzero = 18); cells 22 and 23: stratified log-rank only.
# Required checks (per replicate, recorded): cox_robust() agreement in cells 2, 5, 8, 10, 14, 16 when there is no
# tied time after the time fix; coxph(robust = TRUE, ties = "breslow") agreement in cells 12 and 21, FH(0,0).
# Output: m6_raw/cellNN.rds (replicates x fields), written cell by cell; m6_summarise.R makes the tables.
# Stops and records if the projected wall-clock passes 30 min.
suppressPackageStartupMessages(library(survival))
a <- commandArgs(trailingOnly = TRUE)
NREP <- if (length(a) >= 1) as.integer(a[1]) else 2000L
CORES <- if (length(a) >= 2) as.integer(a[2]) else 12L
if (!file.exists("DESCRIPTION")) stop("Run from the repo root.")
S <- file.path("dev", "audit", "stage2d")
source(file.path(S, "m6_common.R"))
suppressMessages(pkgload::load_all(".", quiet = TRUE, export_all = FALSE, helpers = FALSE, attach_testthat = FALSE))
grDevices::pdf(NULL)
dir.create(file.path(S, "m6_raw"), showWarnings = FALSE)

call_dc <- function(df, wn, g, strata = NULL)
  df_counting(df, tte.name = "time", event.name = "status", treat.name = "treat", weight.name = wn, strata.name = strata,
              scheme = "fh", scheme_params = list(rho = 0, gamma = g))
call_cox <- function(dc, g, draws = 0)
  cox_rhogamma(dc, scheme = "fh", scheme_params = list(rho = 0, gamma = g), draws = draws)
relc <- function(x, y) abs(x / y - 1)
FIELDS <- c("fail", "bhat", "se", "lo", "hi", "z", "zcox", "ties", "chk_se", "chk_z", "chk_zcox", "cph_se", "cph_z",
            "star_fail", "star_se", "star_b", "star_lo", "star_hi",
            "l4_score", "l4_sig2", "l4_dhat", "l4_sig2_dhat", "l4_cov", "l4_cor", "km_sig2_dhat")

one_rep <- function(k, r) {
  cl <- m6_cells[m6_cells$k == k, ]
  d <- m6_data(k, r); df <- d$df; wn <- d$wname
  out <- NULL
  for (g in c(0, 1)) {
    v <- stats::setNames(rep(NA_real_, length(FIELDS)), FIELDS)
    res <- tryCatch({
      dc <- suppressWarnings(call_dc(df, wn, g, if (cl$strata) "s" else NULL))
      v["ties"] <- as.numeric(anyDuplicated(dc$time) > 0)
      if (cl$strata) {
        v["z"] <- dc$z.score_stratified
        if (g == 0) {
          cps <- coxph(Surv(time, status) ~ treat + strata(s), data = df, weights = w, robust = TRUE, ties = "breslow")
          v["cph_z"] <- relc(dc$z.score_stratified^2, cps$rscore)
        }
      } else {
        cx <- suppressWarnings(call_cox(dc, g))
        if (is.null(cx$fit) || !is.finite(cx$fit$bhat)) stop("root finding failed")
        v[c("bhat", "se", "lo", "hi", "z", "zcox")] <- c(cx$fit$bhat, cx$fit$sig_bhat_asy, log(cx$hr_ci_asy$lower), log(cx$hr_ci_asy$upper),
                                                       dc$z.score, cx$fit$z.score)
        if (k %in% c(2, 5, 8, 10, 14, 16)) {
          cr <- cox_robust(cx$fit); cr0 <- cox_robust(cx$fit, 0)
          v[c("chk_se", "chk_z", "chk_zcox")] <- c(relc(cx$fit$sig_bhat_asy, cr$se), relc(abs(dc$z.score), abs(cr0$z)), relc(abs(cx$fit$z.score), abs(cr0$z)))
        }
        if (g == 0 && !is.null(wn)) {
          cp <- coxph(Surv(time, status) ~ treat, data = df, weights = w, robust = TRUE, ties = "breslow")
          v[c("cph_se", "cph_z")] <- c(relc(cx$fit$sig_bhat_asy, sqrt(cp$var[1, 1])), relc(dc$z.score^2, cp$rscore))
        }
        if (g == 0 && r <= 500) {
          rs <- tryCatch({
            set.seed(3e7 + 1e6 * k + r)
            cxs <- suppressWarnings(call_cox(dc, 0, draws = 200))
            c(0, cxs$fit$sig_bhat_star, cxs$fit$bhat_debiased, log(cxs$hr_ci_star$lower), log(cxs$hr_ci_star$upper))
          }, error = function(e) c(1, NA, NA, NA, NA))
          v[c("star_fail", "star_se", "star_b", "star_lo", "star_hi")] <- rs
        }
        if (k %in% c(5, 14)) {
          l4 <- wlr_dhat_estimates(dc, scheme = "fh", scheme_params = list(rho = 0, gamma = g), tzero = 18)
          km <- KM_diff(df, "time", "status", "treat", weight.name = wn, at_points = dc$at_points, draws = 0)
          v[c("l4_score", "l4_sig2", "l4_dhat", "l4_sig2_dhat", "l4_cov", "l4_cor", "km_sig2_dhat")] <-
            c(l4$score, l4$sig2.score, l4$dhat, l4$sig2_dhat, l4$cov_wlr_dhat, l4$cor_wlr_dhat, km$sig2_dhat[max(which(dc$at_points <= 18))])
        }
      }
      0
    }, error = function(e) 1)
    v["fail"] <- res
    out <- c(out, stats::setNames(v, paste0("fh0", g, ".", FIELDS)))
  }
  out
}

t_all <- proc.time()[["elapsed"]]
cat(sprintf("m6_sim.R: weightedsurv %s at %s; reps = %d, cores = %d; started %s\n", as.character(utils::packageVersion("weightedsurv")),
            system("git rev-parse --short HEAD", intern = TRUE), NREP, CORES, format(Sys.time())))
done <- 0
for (k in m6_cells$k) {
  t0 <- proc.time()[["elapsed"]]
  L <- parallel::mclapply(seq_len(NREP), function(r) one_rep(k, r), mc.cores = CORES, mc.preschedule = TRUE)
  bad <- !vapply(L, is.numeric, TRUE)
  if (any(bad)) stop(sprintf("cell %d: %d replicates did not return (%s)", k, sum(bad), as.character(L[[which(bad)[1]]])))
  res <- do.call(rbind, L)
  secs <- proc.time()[["elapsed"]] - t0
  saveRDS(list(k = k, cell = m6_cells$cell[m6_cells$k == k], res = res, secs = secs, cores = CORES,
               sha = system("git rev-parse HEAD", intern = TRUE)), file.path(S, "m6_raw", sprintf("cell%02d.rds", k)))
  done <- done + 1
  tot <- proc.time()[["elapsed"]] - t_all
  proj <- tot / done * nrow(m6_cells)
  cat(sprintf("  cell %2d %-20s %6.1f s; failures FH(0,0) %d, FH(0,1) %d; elapsed %.1f min, projected total %.1f min\n", k,
              m6_cells$cell[m6_cells$k == k], secs, sum(res[, "fh00.fail"]), sum(res[, "fh01.fail"]), tot / 60, proj / 60))
  if (proj > 30 * 60) { cat("STOP: projected wall-clock passes 30 min; recorded and stopped.\n"); quit(status = 2) }
}
cat(sprintf("m6_sim.R: done in %.1f min\n", (proc.time()[["elapsed"]] - t_all) / 60))
