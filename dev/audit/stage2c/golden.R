#!/usr/bin/env Rscript
# =============================================================================
# golden.R -- golden snapshot for Stage 2c (TASK_km_fixes_stage2c.md, Step 1)
# Copy of dev/audit/stage2b/golden.R, changed as follows (datasets and calls are unchanged):
#   * baseline is stage2c/golden_s2b.rds (the branch at the start of Stage 2c);
#   * whitelists are the Stage-2c items R1 and R2; "final" = union over stage2c/packages_done.txt;
#   * contract check: every analytic (draws = 0) sig2_surv0 / sig2_surv1 / sig2_survP that moved
#     equals survfit(Surv(time-fixed time, event) ~ 1, weights, robust = TRUE)$std.err^2 within
#     1e-8 relative (SE scale) wherever survfit's std.err is finite and positive;
#   * contract check: the RMST point estimate (rmst_time, rmst_maxtau_ci[1]) does not move.
#
#   Rscript dev/audit/stage2c/golden.R make                 -> golden_s2b.rds (refuses to overwrite)
#   Rscript dev/audit/stage2c/golden.R compare <whitelist>  -> exit 0 when every field outside the
#                                                             whitelist is identical and the contract checks pass
#
# Identical means: numeric tolerance 1e-12 absolute; NA pattern, names, lengths and strings exact.
# Fields present only in the new run are allowed (elements may be added) and are listed; a list's
# @names may grow when its old names are an unchanged prefix of the new ones.
# No call may raise a warning it did not raise in golden_s2b.rds.
# Local tree via pkgload::load_all. All data seeded.
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[1] else "make"
if (!file.exists("DESCRIPTION")) stop("Run from the repo root.")
suppressMessages(pkgload::load_all(".", quiet = TRUE, export_all = FALSE, helpers = FALSE))
grDevices::pdf(NULL)
S2 <- file.path("dev", "audit", "stage2c")
GOLD <- file.path(S2, "golden_s2b.rds")

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
    }
  }
  warns <- mget(sort(ls(WARN)), envir = WARN)
  list(values = acc, warnings = warns, errors = errors, data = D)
}

# ---- whitelists (regular expressions on flattened paths) -----------------------------
# R1: analytic (draws = 0) calls that pass non-constant weights; variance-derived fields only.
# R2: resampled curves and everything computed from them, in calls that pass non-constant weights:
#     KM_diff_wv_band (star matrices, band), rmst_band_wv (resampling variance, limits, band),
#     KM_diff_weighted_draws on the datasets whose `w` is non-constant (draws = 200 variance path).
WNC <- paste0("(", paste(W_NONCONST, collapse = "|"), ")")
DC_VAR <- "(sig2_surv0|sig2_surv1|sig2_survP|quantile_results/(lower|upper))"
KD_VAR <- "(sig2_surv0|sig2_surv1|sig2_dhat|lower|upper|Zdhat_star|c_alpha_band|sb_lower|sb_upper)"
KD_STAR <- "(surv0_star|surv1_star|dhat_star|Zdhat_star|c_alpha_band|sb_lower|sb_upper)"
RMST_RS <- "(sig2_rmst_time|rmst_time_lower|rmst_time_upper|rmst_maxtau_ci|rmst_text|c_alpha_band|rmst_time_sb_lower|rmst_time_sb_upper)"
WL <- list(
  none = character(0),
  R1 = c(paste0("^", WNC, "/df_counting_weighted/", DC_VAR, "$"),
         paste0("^[^/]+/df_counting_wv/", DC_VAR, "$"),
         paste0("^[^/]+/KM_diff_wv(_band)?/", KD_VAR, "$")),
  R2 = c(paste0("^[^/]+/KM_diff_wv_band/", KD_STAR, "$"),
         paste0("^[^/]+/rmst_band_wv/", RMST_RS, "$"),
         paste0("^", WNC, "/KM_diff_weighted_draws/(sig2_surv0|sig2_surv1|sig2_dhat|lower|upper)$"))
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
    if (d > 1e-12) return(sprintf("max abs diff %.3g", d))
    return(NULL)
  }
  if (!identical(class(a), class(b))) return(sprintf("type %s -> %s", class(a)[1], class(b)[1]))
  if (!identical(a, b)) return("values differ")
  NULL
}

# survfit robust SE on the time-fixed data (the package applies aeqSurv to the times)
robust_se <- function(U, D, W, at) {
  y <- survival::aeqSurv(survival::Surv(U, D))
  f <- survival::survfit(y ~ 1, weights = W, robust = TRUE)
  summary(f, times = at, extend = TRUE)$std.err
}
ANALYTIC_CALLS <- c("df_counting_weighted", "df_counting_wv", "KM_diff_wv", "KM_diff_wv_band")

if (mode == "make") {
  if (file.exists(GOLD)) stop(GOLD, " exists; it is the record and is not overwritten.")
  s <- snapshot()
  if (length(s$errors)) { cat("ERRORS:\n", paste(s$errors, collapse = "\n"), "\n"); quit(status = 1) }
  saveRDS(list(values = s$values, warnings = s$warnings, R = R.version.string,
               created = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
               sha = system("git rev-parse HEAD", intern = TRUE)), GOLD)
  nw <- sum(lengths(s$warnings) > 0)
  cat(sprintf("golden_s2b.rds: %d fields from %d datasets; %d calls recorded, %d of them raise warnings\n",
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
  # warnings: no call may raise a warning it did not raise in golden_s2b.rds
  new_warn <- character(0)
  for (k in names(s$warnings)) {
    extra <- setdiff(s$warnings[[k]], G0$warnings[[k]])
    if (length(extra)) new_warn <- c(new_warn, sprintf("%s: %s", k, paste(extra, collapse = " | ")))
  }
  # R1 contract: every moved analytic sig2_surv* equals survfit robust (SE scale, 1e-8 relative)
  m_checked <- 0L; m_maxrel <- 0
  for (p in grep("/sig2_surv(0|1|P)$", moved_paths, value = TRUE)) {
    parts <- strsplit(p, "/")[[1]]; nm <- parts[1]; call <- parts[2]; fld <- parts[3]
    if (!call %in% ANALYTIC_CALLS) next
    d <- D[[nm]]
    wcol <- if (call == "df_counting_weighted") "w" else "wv"
    arm <- sub("sig2_surv", "", fld)
    sel <- if (arm == "P") rep(TRUE, nrow(d)) else d$treat == as.numeric(arm)
    at <- new[[paste0(nm, "/", call, "/at_points")]]
    rb <- robust_se(d$tte[sel], d$event[sel], d[[wcol]][sel], at)
    got <- sqrt(pmax(new[[p]], 0))
    k <- is.finite(rb) & rb > 0
    rel <- if (any(k)) max(abs(got[k] / rb[k] - 1)) else 0
    m_checked <- m_checked + 1L; m_maxrel <- max(m_maxrel, rel)
    if (!any(k) || rel > 1e-8) contract_fail <- c(contract_fail, sprintf("R1 %s: differs from survfit robust (max rel %.3g, %d points)", p, rel, sum(k)))
  }
  # RMST point estimate must not move (rmst_maxtau_ci is whitelisted for its limits only)
  for (p in grep("/rmst_band(_wv)?/rmst_maxtau_ci$", names(gold), value = TRUE)) {
    if (!p %in% names(new) || abs(gold[[p]][1] - new[[p]][1]) > 1e-12)
      contract_fail <- c(contract_fail, sprintf("%s: RMST point estimate moved", p))
  }
  cat(sprintf("golden compare, whitelist '%s': %d fields in golden, %d in new run, %d added\n", w, length(gold), length(new), length(added)))
  if (length(added)) cat("ADDED (allowed):\n", paste0("  ", added, collapse = "\n"), "\n", sep = "")
  if (length(names_grown)) cat(sprintf("LISTS THAT GAINED ELEMENTS (allowed, %d):\n", length(names_grown)), paste0("  ", names_grown, collapse = "\n"), "\n", sep = "")
  if (length(moved_wl)) cat(sprintf("WHITELISTED fields that moved (%d):\n", length(moved_wl)), paste0("  ", moved_wl, collapse = "\n"), "\n", sep = "")
  if (m_checked) cat(sprintf("R1 check: %d moved analytic sig2_surv* fields checked against survfit robust; max relative SE difference %.3g\n", m_checked, m_maxrel))
  if (length(new_warn)) cat(sprintf("NEW WARNINGS (%d):\n", length(new_warn)), paste0("  ", new_warn, collapse = "\n"), "\n", sep = "")
  else cat("Warnings: no call raises a warning it did not raise in golden_s2b.rds\n")
  if (length(contract_fail)) cat("CONTRACT FAILURES:\n", paste0("  ", contract_fail, collapse = "\n"), "\n", sep = "")
  if (length(diffs)) cat(sprintf("DIFFERENCES outside the whitelist (%d):\n", length(diffs)), paste0("  ", diffs, collapse = "\n"), "\n", sep = "")
  ok <- !length(diffs) && !length(contract_fail) && !length(new_warn)
  cat(if (ok) "GOLDEN: identical outside the whitelist\n" else "GOLDEN: FAIL\n")
  quit(status = if (ok) 0 else 1)
}
stop("usage: golden.R make | compare <whitelist>")
