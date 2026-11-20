#!/usr/bin/env Rscript
# =============================================================================
# golden.R -- golden snapshot for Stage 2e (TASK_km_fixes_stage2e.md, Step 1)
# Copy of dev/audit/stage2d/golden.R, changed as follows (the datasets and the earlier calls are unchanged):
#   * baseline is stage2e/golden_s2d.rds (the branch at the start of Stage 2e);
#   * two weight columns added to every dataset: `w01` (0/1, about 30% zeros) and `wz` (`wv` with about 30% zeros);
#   * calls added (Step 1):
#       on gbsg (tied times, strata), sim_weighted and sim_strata, with `w01` and with `wz`:
#         df_counting() FH(0,0) (zw_dc_*), with strata.name where the dataset has strata (zw_dcs_*),
#         KM_diff(draws.band = 200) (zw_km_*), cox_rhogamma() FH(0,0) draws = 0 and draws = 200 (zw_cox_d0_*, zw_cox_d200_*);
#       check_results(<df_counting result>, verbose = FALSE) without weights and with `wv` (check_results_*);
#       plotKM.band_subgroups() with neither draws nor draws.band supplied, without weights and with `wv` (pkbdef_*);
#   * whitelists are the Stage-2e items Rd, Rb, Rc (Ra and W1-W4 move nothing); "final" = union over stage2e/items_done.txt;
#   * new elements: only those listed in the task document (section 5) may appear, and only when their item is done;
#     a list may gain elements anywhere as long as its earlier names keep their order;
#   * "identical" is exact: numeric fields outside the whitelist must not differ at all (tolerance 0).
#
#   Rscript dev/audit/stage2e/golden.R make                 -> golden_s2d.rds (refuses to overwrite)
#   Rscript dev/audit/stage2e/golden.R compare <whitelist>  -> exit 0 when every field outside the whitelist is
#                                                             identical and no element other than the listed ones is new
#
# NA pattern, names, lengths and strings exact. No call may raise a warning it did not raise in golden_s2d.rds.
# Local tree via pkgload::load_all. All data seeded.
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[1] else "make"
if (!file.exists("DESCRIPTION")) stop("Run from the repo root.")
suppressMessages(pkgload::load_all(".", quiet = TRUE, export_all = FALSE, helpers = FALSE))
grDevices::pdf(NULL)
S2 <- file.path("dev", "audit", "stage2e")
GOLD <- file.path(S2, "golden_s2d.rds")

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
  # Stage 2e: a 0/1 weight column and unequal positive weights with about 30% zeros
  for (i in seq_along(D)) { set.seed(400 + i); keep <- runif(nrow(D[[i]])) < 0.7
    D[[i]]$w01 <- as.numeric(keep); D[[i]]$wz <- ifelse(keep, D[[i]]$wv, 0) }
  D
}
STRATA_OK <- c("veteran", "gbsg", "gbsg_er_gt0", "sim_strata")
# datasets whose `w` is non-constant (df_counting_weighted passes non-constant weights there)
W_NONCONST <- c("veteran", "gbsg", "gbsg_er_le0", "gbsg_er_gt0", "sim_weighted")
# Stage 2e: datasets on which the zero-weight calls run
ZW <- c("gbsg", "sim_weighted", "sim_strata")

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
  # Stage 2e additions
  out$check_results_none <- q("check_results_none", check_results(dc, verbose = FALSE))
  out$check_results_wv <- q("check_results_wv", check_results(out$df_counting_wv, verbose = FALSE))
  out$pkbdef_none <- q("pkbdef_none", plotKM.band_subgroups(d, "tte", "event", "treat"))
  out$pkbdef_wv <- q("pkbdef_wv", plotKM.band_subgroups(d, "tte", "event", "treat", weight.name = "wv"))
  if (nm %in% ZW) for (wt in c("w01", "wz")) {
    k <- function(x) paste0("zw_", x, "_", wt)
    dz <- q(k("dc"), df_counting(d, "tte", "event", "treat", weight.name = wt))
    out[[k("dc")]] <- dz
    if (nm %in% STRATA_OK) out[[k("dcs")]] <- q(k("dcs"), df_counting(d, "tte", "event", "treat", weight.name = wt, strata.name = "st"))
    out[[k("km")]] <- q(k("km"), KM_diff(d, "tte", "event", "treat", weight.name = wt, draws.band = 200))
    out[[k("cox_d0")]] <- q(k("cox_d0"), cox_rhogamma(dz, scheme = "fh", scheme_params = P00, draws = 0))
    out[[k("cox_d200")]] <- q(k("cox_d200"), cox_rhogamma(dz, scheme = "fh", scheme_params = P00, draws = 200))
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
# Rd: a call whose weight column contains a zero (zw_*): variances, standard errors, z-scores, p-value texts,
#     confidence limits, bands and resampled quantities.
# Rb: a call whose positive weights are not all equal: cox_results$score, logrank_results$chisq, $pval,
#     $logrank_text, and the output of check_results().
# Rc: a plotKM.band_subgroups() call that does not supply draws.band (pkbdef_*): the resampled matrices, the band
#     critical value and band limits.
WNC <- paste0("(", paste(W_NONCONST, collapse = "|"), ")")
RB_F <- "(cox_results/score|logrank_results/(chisq|pval|logrank_text))"
WL <- list(
  none = character(0),
  Rd = c("^[^/]+/zw_dcs?_(w01|wz)/(sig2_surv[01P](_stratum)?|sig2_lr(_stratified)?|z\\.score(_stratified)?|zlogrank_text|get_score/(sig2\\.score|z\\.score)|quantile_results/(lower|upper))$",
         "^[^/]+/zw_km_(w01|wz)/(sig2_surv[01]|sig2_dhat|lower|upper|surv[01]_star|dhat_star|Zdhat_star|c_alpha_band|sb_lower|sb_upper)$",
         "^[^/]+/zw_cox_d(0|200)_(w01|wz)/(fit/(sig2_score|sig_bhat_asy|z\\.score|wald|wald_debiased1|wald_debiased2|sig_bhat_star|bhat_debiased)|z\\.score|z\\.score_debiased|hr_ci_asy/(sig_bhat|lower|upper)|hr_ci_star/(beta|sig_bhat|hr|lower|upper)|zlogrank_text|cox_text_asy|cox_text_star|sig2_score_star|fit_resamples/(score_star|score_star_null|bhat_center_star|Z_bstar|sig_bhat_star))$"),
  Rb = c(paste0("^", WNC, "/df_counting_weighted/", RB_F, "$"),
         paste0("^[^/]+/df_counting(_fh01|_strata|_strata_fh01)?_wv/", RB_F, "$"),
         paste0("^[^/]+/zw_dcs?_wz/", RB_F, "$"),
         "^[^/]+/check_results_wv/logrank_chisq$"),
  Rc = c("^[^/]+/pkbdef_(none|wv)/fit_itt/(surv0_star|surv1_star|dhat_star|Zdhat_star)(@dim)?$",
         "^[^/]+/pkbdef_(none|wv)/fit_itt/(c_alpha_band|sb_lower|sb_upper)$"),
  Ra = character(0), W1 = character(0), W2 = character(0), W3 = character(0), W4 = character(0)
)
# New elements allowed, by item: regular expressions on the flattened path of the new field, and the names a
# list may gain (path of the list's @names -> names).
NEW_PATH <- list(
  Rb = c("/(df_counting[^/]*|zw_dcs?_[^/]+)/cox_results/se$", "/(df_counting[^/]*|zw_dcs?_[^/]+)/logrank_results/method$"),
  Ra = c("/(df_counting[^/]*|zw_dcs?_[^/]+)/ess(@names|/(group|n|ess|ratio))$"),
  W2 = c("/(df_counting[^/]*|zw_dcs?_[^/]+|KM_diff[^/]*|zw_km_[^/]+|rmst_band[^/]*)/weight_variance$", "/pkb[^/]*/fit_itt/weight_variance$"),
  W4 = c("/(cox_[^/]+|zw_cox_[^/]+)/weight_variance$")
)
NEW_NAME <- list(
  Rb = list(c("/cox_results@names$", "se"), c("/logrank_results@names$", "method")),
  Ra = list(c("/(df_counting[^/]*|zw_dcs?_[^/]+)@names$", "ess")),
  W2 = list(c("/(df_counting[^/]*|zw_dcs?_[^/]+|KM_diff[^/]*|zw_km_[^/]+|rmst_band[^/]*)@names$", "weight_variance"), c("/pkb[^/]*/fit_itt@names$", "weight_variance")),
  W4 = list(c("/(cox_[^/]+|zw_cox_[^/]+)@names$", "weight_variance"))
)
items_for <- function(w) {
  if (w == "final") {
    f <- file.path(S2, "items_done.txt")
    done <- if (file.exists(f)) trimws(readLines(f)) else character(0)
    done <- done[nzchar(done)]
    message("final whitelist = union over: ", paste(done, collapse = ", "))
    return(done)
  }
  ids <- trimws(strsplit(w, ",")[[1]])
  if (!all(ids %in% names(WL))) stop("unknown whitelist: ", paste(setdiff(ids, names(WL)), collapse = ", "))
  ids
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

if (mode == "make") {
  if (file.exists(GOLD)) stop(GOLD, " exists; it is the record and is not overwritten.")
  s <- snapshot()
  if (length(s$errors)) { cat("ERRORS:\n", paste(s$errors, collapse = "\n"), "\n"); quit(status = 1) }
  saveRDS(list(values = s$values, warnings = s$warnings, R = R.version.string,
               created = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
               sha = system("git rev-parse HEAD", intern = TRUE)), GOLD)
  nw <- sum(lengths(s$warnings) > 0)
  cat(sprintf("golden_s2d.rds: %d fields from %d datasets; %d calls recorded, %d of them raise warnings\n",
              length(s$values), length(s$data), length(s$warnings), nw))
  for (k in names(s$warnings)) if (length(s$warnings[[k]]))
    cat(sprintf("  %s: %s\n", k, paste(substr(s$warnings[[k]], 1, 110), collapse = " | ")))
  quit(status = 0)
}

if (mode == "compare") {
  w <- if (length(args) >= 2) args[2] else "none"
  items <- items_for(w)
  wl <- unique(unlist(WL[items])); new_path <- unique(unlist(NEW_PATH[items])); new_name <- do.call(c, NEW_NAME[items])
  G0 <- readRDS(GOLD); gold <- G0$values
  s <- snapshot(); new <- s$values
  if (length(s$errors)) { cat("ERRORS:\n", paste(s$errors, collapse = "\n"), "\n"); quit(status = 1) }
  in_set <- function(p, set) length(set) > 0 && any(vapply(set, grepl, logical(1), x = p))
  in_wl <- function(p) in_set(p, wl)
  name_ok <- function(p, nm) any(vapply(new_name, function(x) grepl(x[1], p) && nm == x[2], logical(1)))
  diffs <- character(0); moved_wl <- character(0); names_grown <- character(0); bad_new <- character(0)
  for (p in names(gold)) {
    if (!p %in% names(new)) {
      msg <- sprintf("REMOVED %s", p)
      if (in_wl(p)) moved_wl <- c(moved_wl, msg) else diffs <- c(diffs, msg)
      next
    }
    r <- same(gold[[p]], new[[p]])
    if (is.null(r)) next
    # a list that gains elements: its old names keep their order among the new names; each new name must be listed
    if (grepl("@names$", p) && is.character(gold[[p]]) && length(new[[p]]) > length(gold[[p]]) &&
        identical(new[[p]][new[[p]] %in% gold[[p]]], gold[[p]]) && !anyDuplicated(new[[p]])) {
      extra <- setdiff(new[[p]], gold[[p]])
      names_grown <- c(names_grown, sprintf("%s: + %s", p, paste(extra, collapse = ", ")))
      for (e in extra) if (!length(new_name) || !name_ok(p, e)) bad_new <- c(bad_new, sprintf("%s: new element '%s' is not a listed one", p, e))
      next
    }
    if (in_wl(p)) moved_wl <- c(moved_wl, sprintf("%s: %s", p, r))
    else diffs <- c(diffs, sprintf("%s: %s", p, r))
  }
  added <- setdiff(names(new), names(gold))
  for (p in added) if (!in_set(p, new_path)) bad_new <- c(bad_new, sprintf("%s: new field is not a listed one", p))
  # warnings: no call may raise a warning it did not raise in golden_s2d.rds
  new_warn <- character(0)
  for (k in names(s$warnings)) {
    extra <- setdiff(s$warnings[[k]], G0$warnings[[k]])
    if (length(extra)) new_warn <- c(new_warn, sprintf("%s: %s", k, paste(extra, collapse = " | ")))
  }
  cat(sprintf("golden compare, whitelist '%s' (%s): %d fields in golden, %d in new run, %d added\n", w, paste(items, collapse = ","), length(gold), length(new), length(added)))
  if (length(added)) { ta <- table(sub("^[^/]+/[^/]+/", "", added)); cat(sprintf("ADDED fields (%d), by element:\n", length(added)), paste0("  ", names(ta), " x ", ta, collapse = "\n"), "\n", sep = "") }
  if (length(names_grown)) { tg <- table(sub("^.*: ", "", names_grown)); cat(sprintf("LISTS THAT GAINED ELEMENTS (%d):\n", length(names_grown)), paste0("  ", names(tg), " in ", tg, " lists", collapse = "\n"), "\n", sep = "") }
  if (length(moved_wl)) cat(sprintf("WHITELISTED fields that moved (%d):\n", length(moved_wl)), paste0("  ", moved_wl, collapse = "\n"), "\n", sep = "")
  if (length(new_warn)) cat(sprintf("NEW WARNINGS (%d):\n", length(new_warn)), paste0("  ", new_warn, collapse = "\n"), "\n", sep = "")
  else cat("Warnings: no call raises a warning it did not raise in golden_s2d.rds\n")
  if (length(bad_new)) cat(sprintf("NEW ELEMENTS NOT LISTED (%d):\n", length(bad_new)), paste0("  ", bad_new, collapse = "\n"), "\n", sep = "")
  if (length(diffs)) cat(sprintf("DIFFERENCES outside the whitelist (%d):\n", length(diffs)), paste0("  ", diffs, collapse = "\n"), "\n", sep = "")
  ok <- !length(diffs) && !length(bad_new) && !length(new_warn)
  cat(if (ok) "GOLDEN: identical outside the whitelist\n" else "GOLDEN: FAIL\n")
  quit(status = if (ok) 0 else 1)
}
stop("usage: golden.R make | compare <whitelist>")
