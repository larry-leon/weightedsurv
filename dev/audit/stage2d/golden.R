#!/usr/bin/env Rscript
# =============================================================================
# golden.R -- golden snapshot for Stage 2d (TASK_km_fixes_stage2d.md, Step 1)
# Copy of dev/audit/stage2c/golden.R, changed as follows (the datasets and the Stage-2c calls are unchanged):
#   * baseline is stage2d/golden_s2c.rds (the branch at the start of Stage 2d);
#   * calls added (Step 1), each with weights `wv`, with no weights and with a constant weight column `wc` (= 2):
#       df_counting() FH(0,0) and FH(0,1), and with strata.name on the datasets with strata;
#       wlr_cumulative(<df_counting result>, return_cumulative = TRUE) FH(0,0) and FH(0,1), and
#       wlr_cumulative() on a bare data frame of counts (cumulative and not);
#       wlr_dhat_estimates(<df_counting result>, tzero = median of the grid) FH(0,0) and FH(0,1);
#       cox_rhogamma() draws = 0 for FH(0,0), FH(0,1), MB, schemper; draws = 200 for FH(0,0);
#       plotKM.band_subgroups(draws.band = 200) with `draws` not supplied (null device);
#   * whitelists are the Stage-2d items L1, L2, L3, L4, B12; "final" = union over stage2d/packages_done.txt;
#   * "identical" is exact: numeric fields outside the whitelist must not differ at all (tolerance 0);
#   * contract checks: every FH(0,0) z.score (df_counting, unstratified and stratified; cox_rhogamma) and
#     sig_bhat_asy that moved is compared with coxph(weights = , robust = TRUE, ties = "breslow"):
#     z^2 against $rscore within 1e-8 relative, sig_bhat_asy against sqrt($var) within 1e-6 relative.
#
#   Rscript dev/audit/stage2d/golden.R make                 -> golden_s2c.rds (refuses to overwrite)
#   Rscript dev/audit/stage2d/golden.R compare <whitelist>  -> exit 0 when every field outside the
#                                                             whitelist is identical and the contract checks pass
#
# NA pattern, names, lengths and strings exact. Fields present only in the new run are allowed and listed.
# No call may raise a warning it did not raise in golden_s2c.rds.
# Local tree via pkgload::load_all. All data seeded.
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[1] else "make"
if (!file.exists("DESCRIPTION")) stop("Run from the repo root.")
suppressMessages(pkgload::load_all(".", quiet = TRUE, export_all = FALSE, helpers = FALSE))
grDevices::pdf(NULL)
S2 <- file.path("dev", "audit", "stage2d")
GOLD <- file.path(S2, "golden_s2c.rds")

# ---- datasets -------------------------------------------------------------------
sim <- function(n, seed, ties = FALSE, weighted = FALSE, nstrata = 1) {
  set.seed(seed)
  d <- data.frame(treat = rbinom(n, 1, 0.5))
  d$st <- if (nstrata > 1) sample(paste0("S", seq_len(nstrata)), n, replace = TRUE) else "S1"
  rate <- 0.08 * exp(-0.3 * d$treat) * if (nstrata > 1) c(S1 = 1, S2 = 1.5, S3 = 0.7)[d$st] else 1
  d$tte <- rexp(n, rate); cc <- runif(n, 0, 30)
  d$event <- as.numeric(d$tte <= cc); d$tte <- pmin(d$tte, cc)
  if (ties) d$tte <- ceiling(d$tte * 2) / 2
  d$w <- if (weighted) runif(n, 0.3, 3) else 1
  d
}
add_w <- function(d, seed) { set.seed(seed); d$w <- runif(nrow(d), 0.3, 3); d }
add_wv <- function(d, seed) {
  if (length(unique(d$w)) > 1) d$wv <- d$w else { set.seed(seed); d$wv <- runif(nrow(d), 0.3, 3) }
  d
}
make_data <- function() {
  v <- survival::veteran
  v <- data.frame(tte = v$time, event = v$status, treat = v$trt - 1, st = as.character(v$celltype))
  g <- survival::gbsg
  g <- data.frame(tte = g$rfstime / 30.4375, event = g$status, treat = g$hormon, st = paste0("G", g$grade), er = g$er)
  D <- list(
    veteran      = add_w(v, 101),
    gbsg         = add_w(g, 102),
    gbsg_er_le0  = add_w(g[g$er <= 0, ], 103),
    gbsg_er_gt0  = add_w(g[g$er > 0, ], 104),
    sim_cont     = sim(200, 201),
    sim_ties     = sim(200, 202, ties = TRUE),
    sim_weighted = sim(300, 203, weighted = TRUE),
    sim_strata   = sim(300, 204, nstrata = 3)
  )
  for (i in seq_along(D)) D[[i]] <- add_wv(D[[i]], 300 + i)
  for (i in seq_along(D)) D[[i]]$wc <- 2          # Stage 2d: a constant weight column
  D
}
STRATA_OK <- c("veteran", "gbsg", "gbsg_er_gt0", "sim_strata")
# datasets whose `w` is non-constant (df_counting_weighted passes non-constant weights there)
W_NONCONST <- c("veteran", "gbsg", "gbsg_er_le0", "gbsg_er_gt0", "sim_weighted")

# ---- calls ------------------------------------------------------------------------
WARN <- new.env()
rec <- function(key, expr) {
  msgs <- character(0)
  val <- withCallingHandlers(expr,
    warning = function(w) { msgs <<- c(msgs, conditionMessage(w)); invokeRestart("muffleWarning") },
    message = function(m) invokeRestart("muffleMessage"))
  assign(key, sort(unique(msgs)), envir = WARN)
  val
}
run_calls <- function(d, nm) {
  out <- list()
  q <- function(call, expr) rec(paste0(nm, "/", call), expr)
  tz <- stats::median(d$tte); tstar <- unname(stats::quantile(d$tte, 0.5))
  dc <- q("df_counting_default", df_counting(d, "tte", "event", "treat"))
  out$df_counting_default <- dc
  out$df_counting_weighted <- q("df_counting_weighted", df_counting(d, "tte", "event", "treat", weight.name = "w"))
  out$df_counting_MB <- q("df_counting_MB", df_counting(d, "tte", "event", "treat", scheme = "MB", scheme_params = list(mb_tstar = tstar)))
  out$df_counting_q25 <- q("df_counting_q25", df_counting(d, "tte", "event", "treat", qprob = 0.25))
  out$KM_diff_default <- q("KM_diff_default", KM_diff(d, "tte", "event", "treat"))
  kb <- q("KM_diff_band", KM_diff(d, "tte", "event", "treat", draws.band = 200, show_resamples = FALSE))
  out$KM_diff_band <- kb
  out$KM_diff_band_modtau <- q("KM_diff_band_modtau", KM_diff(d, "tte", "event", "treat", draws.band = 200, show_resamples = FALSE, modify_tau = TRUE))
  out$KM_diff_weighted_draws <- q("KM_diff_weighted_draws", KM_diff(d, "tte", "event", "treat", weight.name = "w", draws = 200))
  out$cox_fh00_d0 <- q("cox_fh00_d0", cox_rhogamma(dc, scheme = "fh", scheme_params = list(rho = 0, gamma = 0), draws = 0))
  out$cox_fh01_d0 <- q("cox_fh01_d0", cox_rhogamma(dc, scheme = "fh", scheme_params = list(rho = 0, gamma = 1), draws = 0))
  out$cox_fh00_d200 <- q("cox_fh00_d200", cox_rhogamma(dc, scheme = "fh", scheme_params = list(rho = 0, gamma = 0), draws = 200))
  out$cox_fh01_d200 <- q("cox_fh01_d200", cox_rhogamma(dc, scheme = "fh", scheme_params = list(rho = 0, gamma = 1), draws = 200))
  out$wlr_dhat <- q("wlr_dhat", wlr_dhat_estimates(dc, tzero = tz))
  out$rmst_band <- q("rmst_band", cumulative_rmst_bands(d, kb, "tte", "event", "treat", plot = FALSE))
  kd <- out$KM_diff_default
  out$km_quantile_arm0 <- q("km_quantile_arm0", km_quantile(kd$at_points, kd$surv0, sqrt(kd$sig2_surv0)))
  out$km_quantile_arm1 <- q("km_quantile_arm1", km_quantile(kd$at_points, kd$surv1, sqrt(kd$sig2_surv1)))
  # Stage 2b additions: non-constant weights, draws = 0
  out$df_counting_wv <- q("df_counting_wv", df_counting(d, "tte", "event", "treat", weight.name = "wv"))
  out$KM_diff_wv <- q("KM_diff_wv", KM_diff(d, "tte", "event", "treat", weight.name = "wv"))
  kbw <- q("KM_diff_wv_band", KM_diff(d, "tte", "event", "treat", weight.name = "wv", draws.band = 200, show_resamples = FALSE))
  out$KM_diff_wv_band <- kbw
  out$rmst_band_wv <- q("rmst_band_wv", cumulative_rmst_bands(d, kbw, "tte", "event", "treat", weight.name = "wv", plot = FALSE))
  # Stage 2d additions: log-rank and Cox, with weights wv, no weights (none) and a constant weight column (wc)
  P00 <- list(rho = 0, gamma = 0); P01 <- list(rho = 0, gamma = 1)
  out$df_counting_wc <- q("df_counting_wc", df_counting(d, "tte", "event", "treat", weight.name = "wc"))
  for (wt in c("none", "wv", "wc")) {
    wn <- if (wt == "none") NULL else wt
    k <- function(x) paste0(x, "_", wt)
    dc0 <- switch(wt, none = dc, wv = out$df_counting_wv, wc = out$df_counting_wc)
    out[[k("df_counting_fh01")]] <- q(k("df_counting_fh01"), df_counting(d, "tte", "event", "treat", weight.name = wn, scheme = "fh", scheme_params = P01))
    out[[k("wlr_cum_fh00")]] <- q(k("wlr_cum_fh00"), wlr_cumulative(dc0, "fh", P00, return_cumulative = TRUE))
    out[[k("wlr_cum_fh01")]] <- q(k("wlr_cum_fh01"), wlr_cumulative(dc0, "fh", P01, return_cumulative = TRUE))
    bare <- data.frame(at_points = dc0$at_points, nbar0 = dc0$nbar0, nbar1 = dc0$nbar1, ybar0 = dc0$ybar0, ybar1 = dc0$ybar1,
                       survP = dc0$survP, survG = dc0$survG)
    out[[k("wlr_cum_bare")]] <- q(k("wlr_cum_bare"), wlr_cumulative(bare, "fh", P00, return_cumulative = TRUE))
    out[[k("wlr_bare")]] <- q(k("wlr_bare"), wlr_cumulative(bare, "fh", P01, return_cumulative = FALSE))
    tzg <- stats::median(dc0$at_points)
    out[[k("wlr_dhat_grid")]] <- q(k("wlr_dhat_grid"), wlr_dhat_estimates(dc0, tzero = tzg))
    out[[k("wlr_dhat_fh01_grid")]] <- q(k("wlr_dhat_fh01_grid"), wlr_dhat_estimates(dc0, scheme = "fh", scheme_params = P01, tzero = tzg))
    out[[k("cox_fh00_d0")]] <- q(k("cox_fh00_d0"), cox_rhogamma(dc0, scheme = "fh", scheme_params = P00, draws = 0))
    out[[k("cox_fh01_d0")]] <- q(k("cox_fh01_d0"), cox_rhogamma(dc0, scheme = "fh", scheme_params = P01, draws = 0))
    out[[k("cox_MB_d0")]] <- q(k("cox_MB_d0"), cox_rhogamma(dc0, scheme = "MB", scheme_params = list(mb_tstar = tstar), draws = 0))
    out[[k("cox_schemper_d0")]] <- q(k("cox_schemper_d0"), cox_rhogamma(dc0, scheme = "schemper", scheme_params = list(), draws = 0))
    out[[k("cox_fh00_d200")]] <- q(k("cox_fh00_d200"), cox_rhogamma(dc0, scheme = "fh", scheme_params = P00, draws = 200))
    out[[k("pkb")]] <- q(k("pkb"), plotKM.band_subgroups(d, "tte", "event", "treat", weight.name = wn, draws.band = 200))
  }
  out
}

# ---- flatten: every numeric / character / logical leaf, plus the names of every node ----
flatten <- function(x, path, acc) {
  if (is.environment(x) || is.function(x) || is.language(x) || is.null(x)) return(acc)
  if (is.factor(x)) x <- as.character(x)
  if (is.list(x)) {
    nm <- names(x)
    acc[[paste0(path, "@names")]] <- if (is.null(nm)) character(0) else nm
    for (i in seq_along(x)) {
      key <- if (!is.null(nm) && nzchar(nm[i])) nm[i] else paste0("[[", i, "]]")
      acc <- flatten(x[[i]], paste0(path, "/", key), acc)
    }
    return(acc)
  }
  if (is.numeric(x) || is.character(x) || is.logical(x)) {
    v <- x; attributes(v) <- NULL
    acc[[path]] <- v
    if (!is.null(names(x))) acc[[paste0(path, "@names")]] <- names(x)
    if (!is.null(dim(x))) acc[[paste0(path, "@dim")]] <- dim(x)
  }
  acc
}
snapshot <- function() {
  D <- make_data()
  acc <- list(); errors <- character(0)
  for (nm in names(D)) {
    d <- D[[nm]]
    r <- tryCatch(run_calls(d, nm), error = function(e) { errors <<- c(errors, paste(nm, conditionMessage(e))); NULL })
    if (!is.null(r)) acc <- flatten(r, nm, acc)
    if (nm %in% STRATA_OK) {
      s <- tryCatch(rec(paste0(nm, "/df_counting_strata"), df_counting(d, "tte", "event", "treat", strata.name = "st")),
                    error = function(e) { errors <<- c(errors, paste(nm, "strata:", conditionMessage(e))); NULL })
      if (!is.null(s)) acc <- flatten(list(df_counting_strata = s), nm, acc)
      # Stage 2d additions: stratified, FH(0,0) with wv / wc, FH(0,1) with none / wv / wc
      for (cs in list(list("df_counting_strata_wv", "wv", 0), list("df_counting_strata_wc", "wc", 0),
                      list("df_counting_strata_fh01_none", NULL, 1), list("df_counting_strata_fh01_wv", "wv", 1),
                      list("df_counting_strata_fh01_wc", "wc", 1))) {
        s <- tryCatch(rec(paste0(nm, "/", cs[[1]]), df_counting(d, "tte", "event", "treat", weight.name = cs[[2]], strata.name = "st",
                                                              scheme = "fh", scheme_params = list(rho = 0, gamma = cs[[3]]))),
                      error = function(e) { errors <<- c(errors, paste(nm, cs[[1]], conditionMessage(e))); NULL })
        if (!is.null(s)) acc <- flatten(stats::setNames(list(s), cs[[1]]), nm, acc)
      }
    }
  }
  warns <- mget(sort(ls(WARN)), envir = WARN)
  list(values = acc, warnings = warns, errors = errors, data = D)
}

# ---- whitelists (regular expressions on flattened paths) -----------------------------
# A field may move only in a call that passes non-constant subject weights (suffix _wv, and
# df_counting_weighted on the datasets whose `w` is non-constant), and only if it is listed for a DONE item.
# L2 note: wald_debiased1 (= bhat_debiased / sig_bhat_asy), z.score_debiased (= debiased score / sqrt(sig2_score))
# and fit_resamples/Z_bstar (= bhat_center_star / sig_bhat_asy) are listed under L3 in the task document, but they
# are divided by an L2 quantity and so already move with L2; they are whitelisted from L2 on.
WNC <- paste0("(", paste(W_NONCONST, collapse = "|"), ")")
LR_F <- "(z\\.score|sig2_lr|sig2_lr_stratified|z\\.score_stratified|zlogrank_text|get_score/(sig2\\.score|z\\.score))"
WL <- list(
  none = character(0),
  L1 = c(paste0("^", WNC, "/df_counting_weighted/", LR_F, "$"),
         paste0("^[^/]+/df_counting(_fh01|_strata|_strata_fh01)?_wv/", LR_F, "$"),
         "^[^/]+/wlr_cum_fh0[01]_wv/(sig2\\.score|z\\.score)$"),
  L2 = c("^[^/]+/cox_(fh00|fh01|MB|schemper)_d(0|200)_wv/(fit/(sig2_score|sig_bhat_asy|z\\.score|wald)|z\\.score|hr_ci_asy/(sig_bhat|lower|upper)|zlogrank_text|cox_text_asy)$",
         "^[^/]+/cox_fh00_d200_wv/(fit/wald_debiased1|z\\.score_debiased|fit_resamples/Z_bstar)$"),
  L3 = c("^[^/]+/cox_fh00_d200_wv/(fit_resamples/(score_star|score_star_null|bhat_center_star|Z_bstar|sig_bhat_star)|fit/(sig_bhat_star|bhat_debiased|wald_debiased1|wald_debiased2)|hr_ci_star/(beta|sig_bhat|hr|lower|upper)|sig2_score_star|z\\.score_debiased|cox_text_star)$"),
  L4 = c("^[^/]+/wlr_dhat(_fh01)?_grid_wv/(sig2\\.score|z\\.score|sig2_dhat|cov_wlr_dhat|cor_wlr_dhat)$"),
  B12 = c("^[^/]+/pkb_wv/(fit_itt/(sig2_surv0|sig2_surv1|sig2_dhat|lower|upper|Zdhat_star|c_alpha_band|sb_lower|sb_upper)|lower|upper)$")
)
whitelist_for <- function(w) {
  if (w == "final") {
    f <- file.path(S2, "packages_done.txt")
    done <- if (file.exists(f)) trimws(readLines(f)) else character(0)
    done <- done[nzchar(done)]
    message("final whitelist = union over: ", paste(done, collapse = ", "))
    return(unique(unlist(WL[done])))
  }
  ids <- trimws(strsplit(w, ",")[[1]])
  if (!all(ids %in% names(WL))) stop("unknown whitelist: ", paste(setdiff(ids, names(WL)), collapse = ", "))
  unique(unlist(WL[ids]))
}

same <- function(a, b) {
  if (length(a) != length(b)) return(sprintf("length %d -> %d", length(a), length(b)))
  if (is.numeric(a) && is.numeric(b)) {
    if (!identical(is.na(a), is.na(b))) return("NA pattern differs")
    ok <- !is.na(a)
    inf <- ok & (is.infinite(a) | is.infinite(b))
    if (any(inf & (a != b))) return("infinite values differ")
    fin <- ok & !inf
    d <- if (any(fin)) max(abs(a[fin] - b[fin])) else 0
    if (d > 0) return(sprintf("max abs diff %.3g", d))     # Stage 2d: exact
    return(NULL)
  }
  if (!identical(class(a), class(b))) return(sprintf("type %s -> %s", class(a)[1], class(b)[1]))
  if (!identical(a, b)) return("values differ")
  NULL
}

# coxph reference for the FH(0,0) checks: robust score test and robust SE, Breslow ties
cox_ref <- function(d, wcol, strat = FALSE) {
  f <- if (strat) survival::Surv(tte, event) ~ treat + survival::strata(st) else survival::Surv(tte, event) ~ treat
  d$.w <- d[[wcol]]
  cp <- survival::coxph(f, data = d, weights = .w, robust = TRUE, ties = "breslow")
  list(rscore = unname(cp$rscore), se = unname(sqrt(cp$var[1, 1])))
}

if (mode == "make") {
  if (file.exists(GOLD)) stop(GOLD, " exists; it is the record and is not overwritten.")
  s <- snapshot()
  if (length(s$errors)) { cat("ERRORS:\n", paste(s$errors, collapse = "\n"), "\n"); quit(status = 1) }
  saveRDS(list(values = s$values, warnings = s$warnings, R = R.version.string,
               created = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
               sha = system("git rev-parse HEAD", intern = TRUE)), GOLD)
  nw <- sum(lengths(s$warnings) > 0)
  cat(sprintf("golden_s2c.rds: %d fields from %d datasets; %d calls recorded, %d of them raise warnings\n",
              length(s$values), length(s$data), length(s$warnings), nw))
  for (k in names(s$warnings)) if (length(s$warnings[[k]]))
    cat(sprintf("  %s: %s\n", k, paste(substr(s$warnings[[k]], 1, 110), collapse = " | ")))
  quit(status = 0)
}

if (mode == "compare") {
  w <- if (length(args) >= 2) args[2] else "none"
  wl <- whitelist_for(w)
  G0 <- readRDS(GOLD); gold <- G0$values
  s <- snapshot(); new <- s$values
  if (length(s$errors)) { cat("ERRORS:\n", paste(s$errors, collapse = "\n"), "\n"); quit(status = 1) }
  in_wl <- function(p) length(wl) > 0 && any(vapply(wl, grepl, logical(1), x = p))
  diffs <- character(0); moved_wl <- character(0); moved_paths <- character(0); contract_fail <- character(0); names_grown <- character(0)
  for (p in names(gold)) {
    if (!p %in% names(new)) {
      msg <- sprintf("REMOVED %s", p)
      if (in_wl(p)) moved_wl <- c(moved_wl, msg) else diffs <- c(diffs, msg)
      next
    }
    r <- same(gold[[p]], new[[p]])
    if (is.null(r)) next
    # a list that gains elements: its old names are an unchanged prefix of the new names (allowed)
    if (grepl("@names$", p) && is.character(gold[[p]]) && length(new[[p]]) > length(gold[[p]]) &&
        identical(new[[p]][seq_along(gold[[p]])], gold[[p]])) {
      names_grown <- c(names_grown, sprintf("%s: + %s", p, paste(new[[p]][-seq_along(gold[[p]])], collapse = ", ")))
      next
    }
    if (in_wl(p)) { moved_wl <- c(moved_wl, sprintf("%s: %s", p, r)); moved_paths <- c(moved_paths, p) }
    else diffs <- c(diffs, sprintf("%s: %s", p, r))
  }
  added <- setdiff(names(new), names(gold))
  D <- s$data
  # warnings: no call may raise a warning it did not raise in golden_s2c.rds
  new_warn <- character(0)
  for (k in names(s$warnings)) {
    extra <- setdiff(s$warnings[[k]], G0$warnings[[k]])
    if (length(extra)) new_warn <- c(new_warn, sprintf("%s: %s", k, paste(extra, collapse = " | ")))
  }
  # Contract: every moved FH(0,0) z.score / sig_bhat_asy against coxph(robust = TRUE, ties = "breslow")
  n_z <- 0L; n_se <- 0L; max_z <- 0; max_se <- 0
  for (p in moved_paths) {
    parts <- strsplit(p, "/")[[1]]; nm <- parts[1]; call <- parts[2]; fld <- paste(parts[-(1:2)], collapse = "/")
    d <- D[[nm]]; ref <- NULL; kind <- NULL
    if (call %in% c("df_counting_wv", "df_counting_weighted") && fld == "z.score") { ref <- cox_ref(d, if (call == "df_counting_wv") "wv" else "w"); kind <- "z" }
    if (call == "df_counting_strata_wv" && fld == "z.score_stratified") { ref <- cox_ref(d, "wv", strat = TRUE); kind <- "z" }
    if (call %in% c("cox_fh00_d0_wv", "cox_fh00_d200_wv") && fld == "fit/z.score") { ref <- cox_ref(d, "wv"); kind <- "z" }
    if (call %in% c("cox_fh00_d0_wv", "cox_fh00_d200_wv") && fld == "fit/sig_bhat_asy") { ref <- cox_ref(d, "wv"); kind <- "se" }
    if (is.null(kind)) next
    if (kind == "z") {
      rel <- abs(new[[p]]^2 / ref$rscore - 1); n_z <- n_z + 1L; max_z <- max(max_z, rel)
      if (!is.finite(rel) || rel > 1e-8) contract_fail <- c(contract_fail, sprintf("%s: z^2 %.10g differs from coxph rscore %.10g (rel %.3g)", p, new[[p]]^2, ref$rscore, rel))
    } else {
      rel <- abs(new[[p]] / ref$se - 1); n_se <- n_se + 1L; max_se <- max(max_se, rel)
      if (!is.finite(rel) || rel > 1e-6) contract_fail <- c(contract_fail, sprintf("%s: %.10g differs from coxph robust SE %.10g (rel %.3g)", p, new[[p]], ref$se, rel))
    }
  }
  cat(sprintf("golden compare, whitelist '%s': %d fields in golden, %d in new run, %d added\n", w, length(gold), length(new), length(added)))
  if (length(added)) cat("ADDED (allowed):\n", paste0("  ", added, collapse = "\n"), "\n", sep = "")
  if (length(names_grown)) cat(sprintf("LISTS THAT GAINED ELEMENTS (allowed, %d):\n", length(names_grown)), paste0("  ", names_grown, collapse = "\n"), "\n", sep = "")
  if (length(moved_wl)) cat(sprintf("WHITELISTED fields that moved (%d):\n", length(moved_wl)), paste0("  ", moved_wl, collapse = "\n"), "\n", sep = "")
  if (n_z + n_se) cat(sprintf("coxph check: %d moved FH(0,0) z-scores against coxph robust score test (max relative difference of z^2 %.3g); %d moved sig_bhat_asy against coxph robust SE (max relative difference %.3g)\n", n_z, max_z, n_se, max_se))
  if (length(new_warn)) cat(sprintf("NEW WARNINGS (%d):\n", length(new_warn)), paste0("  ", new_warn, collapse = "\n"), "\n", sep = "")
  else cat("Warnings: no call raises a warning it did not raise in golden_s2c.rds\n")
  if (length(contract_fail)) cat("CONTRACT FAILURES:\n", paste0("  ", contract_fail, collapse = "\n"), "\n", sep = "")
  if (length(diffs)) cat(sprintf("DIFFERENCES outside the whitelist (%d):\n", length(diffs)), paste0("  ", diffs, collapse = "\n"), "\n", sep = "")
  ok <- !length(diffs) && !length(contract_fail) && !length(new_warn)
  cat(if (ok) "GOLDEN: identical outside the whitelist\n" else "GOLDEN: FAIL\n")
  quit(status = if (ok) 0 else 1)
}
stop("usage: golden.R make | compare <whitelist>")
