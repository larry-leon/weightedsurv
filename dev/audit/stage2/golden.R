#!/usr/bin/env Rscript
# =============================================================================
# golden.R -- golden snapshot of weightedsurv results (TASK_km_fixes_stage2.md, Step 2a)
#
#   Rscript dev/audit/stage2/golden.R make                 -> golden_before.rds (refuses to overwrite)
#   Rscript dev/audit/stage2/golden.R compare <whitelist>  -> exit 0 when every field outside the
#                                                            whitelist is identical to golden_before.rds
#   <whitelist> is one package ID (B7, B6, V, B4, B2, B1, B3, B5, B10), "none", or "final"
#   (the union over the packages listed in dev/audit/stage2/packages_done.txt).
#
# Identical means: numeric tolerance 1e-12 absolute; NA pattern, names, lengths and strings exact.
# Fields present only in the new run are allowed (elements may be added) and are listed.
# Whitelisted fields that move are listed with their extra contract check:
#   B1: surv0 / surv1 of the modify_tau fits must equal survfit at the grid (1e-10);
#   B3: each moved quantile value must equal survfit's (1e-10, NA pattern exact);
#   B5: the old/new ratio of each moved wlr_dhat_estimates() field is printed.
# Local tree via pkgload::load_all. All data seeded.
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[1] else "make"
if (!file.exists("DESCRIPTION")) stop("Run from the repo root.")
suppressMessages(pkgload::load_all(".", quiet = TRUE, export_all = FALSE, helpers = FALSE))
grDevices::pdf(NULL)
S2 <- file.path("dev", "audit", "stage2")
GOLD <- file.path(S2, "golden_before.rds")

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
make_data <- function() {
  v <- survival::veteran
  v <- data.frame(tte = v$time, event = v$status, treat = v$trt - 1, st = as.character(v$celltype))
  g <- survival::gbsg
  g <- data.frame(tte = g$rfstime / 30.4375, event = g$status, treat = g$hormon, st = paste0("G", g$grade), er = g$er)
  list(
    veteran      = add_w(v, 101),
    gbsg         = add_w(g, 102),
    gbsg_er_le0  = add_w(g[g$er <= 0, ], 103),
    gbsg_er_gt0  = add_w(g[g$er > 0, ], 104),
    sim_cont     = sim(200, 201),
    sim_ties     = sim(200, 202, ties = TRUE),
    sim_weighted = sim(300, 203, weighted = TRUE),
    sim_strata   = sim(300, 204, nstrata = 3)
  )
}
# strata column per dataset (a stratum with one arm empty would stop 0.1.0, so the
# gbsg subsets use grade; er <= 0 is left out because grade 1 has no treated patient there)
STRATA_OK <- c("veteran", "gbsg", "gbsg_er_gt0", "sim_strata")

# ---- calls ------------------------------------------------------------------------
quiet <- function(expr) suppressWarnings(suppressMessages(expr))
run_calls <- function(d) {
  out <- list()
  tz <- stats::median(d$tte); tstar <- unname(stats::quantile(d$tte, 0.5))
  dc <- quiet(df_counting(d, "tte", "event", "treat"))
  out$df_counting_default <- dc
  out$df_counting_weighted <- quiet(df_counting(d, "tte", "event", "treat", weight.name = "w"))
  out$df_counting_MB <- quiet(df_counting(d, "tte", "event", "treat", scheme = "MB", scheme_params = list(mb_tstar = tstar)))
  out$df_counting_q25 <- quiet(df_counting(d, "tte", "event", "treat", qprob = 0.25))
  out$KM_diff_default <- quiet(KM_diff(d, "tte", "event", "treat"))
  kb <- quiet(KM_diff(d, "tte", "event", "treat", draws.band = 200, show_resamples = FALSE))
  out$KM_diff_band <- kb
  out$KM_diff_band_modtau <- quiet(KM_diff(d, "tte", "event", "treat", draws.band = 200, show_resamples = FALSE, modify_tau = TRUE))
  out$KM_diff_weighted_draws <- quiet(KM_diff(d, "tte", "event", "treat", weight.name = "w", draws = 200))
  out$cox_fh00_d0 <- quiet(cox_rhogamma(dc, scheme = "fh", scheme_params = list(rho = 0, gamma = 0), draws = 0))
  out$cox_fh01_d0 <- quiet(cox_rhogamma(dc, scheme = "fh", scheme_params = list(rho = 0, gamma = 1), draws = 0))
  out$cox_fh00_d200 <- quiet(cox_rhogamma(dc, scheme = "fh", scheme_params = list(rho = 0, gamma = 0), draws = 200))
  out$cox_fh01_d200 <- quiet(cox_rhogamma(dc, scheme = "fh", scheme_params = list(rho = 0, gamma = 1), draws = 200))
  out$wlr_dhat <- quiet(wlr_dhat_estimates(dc, tzero = tz))
  out$rmst_band <- quiet(cumulative_rmst_bands(d, kb, "tte", "event", "treat", plot = FALSE))
  kd <- out$KM_diff_default
  out$km_quantile_arm0 <- quiet(km_quantile(kd$at_points, kd$surv0, sqrt(kd$sig2_surv0)))
  out$km_quantile_arm1 <- quiet(km_quantile(kd$at_points, kd$surv1, sqrt(kd$sig2_surv1)))
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
    r <- tryCatch(run_calls(d), error = function(e) { errors <<- c(errors, paste(nm, conditionMessage(e))); NULL })
    if (!is.null(r)) acc <- flatten(r, nm, acc)
    if (nm %in% STRATA_OK) {
      s <- tryCatch(quiet(df_counting(d, "tte", "event", "treat", strata.name = "st")),
                    error = function(e) { errors <<- c(errors, paste(nm, "strata:", conditionMessage(e))); NULL })
      if (!is.null(s)) acc <- flatten(list(df_counting_strata = s), nm, acc)
    }
  }
  list(values = acc, errors = errors, data = D)
}

# ---- whitelists (regular expressions on flattened paths) -----------------------------
WL <- list(
  none = character(0), B7 = character(0), B6 = character(0), V = character(0),
  B4 = character(0), B2 = character(0), B10 = character(0),
  B1 = c("^[^/]+/KM_diff_band_modtau/"),
  B3 = c("^[^/]+/df_counting_[^/]+/quantile_results/", "^[^/]+/df_counting_[^/]+/quantile_check/",
         "^[^/]+/km_quantile_arm[01]/"),
  B5 = c("^[^/]+/wlr_dhat/(sig2_dhat|cov_wlr_dhat|cor_wlr_dhat)$")
)
whitelist_for <- function(w) {
  if (w == "final") {
    f <- file.path(S2, "packages_done.txt")
    done <- if (file.exists(f)) trimws(readLines(f)) else character(0)
    done <- done[nzchar(done)]
    message("final whitelist = union over: ", paste(done, collapse = ", "))
    return(unique(unlist(WL[done])))
  }
  if (!w %in% names(WL)) stop("unknown whitelist: ", w)
  WL[[w]]
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

# extra contract checks for whitelisted fields that moved
sf_at <- function(d, g, times) {
  dd <- d[d$treat == g, ]
  summary(survival::survfit(survival::Surv(dd$tte, dd$event) ~ 1), times = times, extend = TRUE)$surv
}
sf_q <- function(d, g, qprob, cl = 0.95, weighted = FALSE) {
  dd <- d[d$treat == g, ]
  fit <- if (weighted) survival::survfit(survival::Surv(dd$tte, dd$event) ~ 1, weights = dd$w, conf.int = cl, conf.type = "log")
         else survival::survfit(survival::Surv(dd$tte, dd$event) ~ 1, conf.int = cl, conf.type = "log")
  q <- stats::quantile(fit, probs = 1 - qprob)
  c(time = unname(q$quantile), lower = unname(q$lower), upper = unname(q$upper))
}
cmpv <- function(a, b) all(is.na(a) == is.na(b)) && all(abs(a - b)[!is.na(a)] <= 1e-10)

if (mode == "make") {
  if (file.exists(GOLD)) stop(GOLD, " exists; it is the record and is not overwritten.")
  s <- snapshot()
  if (length(s$errors)) { cat("ERRORS:\n", paste(s$errors, collapse = "\n"), "\n"); quit(status = 1) }
  saveRDS(list(values = s$values, R = R.version.string, created = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
               sha = system("git rev-parse HEAD", intern = TRUE)), GOLD)
  cat(sprintf("golden_before.rds: %d fields from %d datasets\n", length(s$values), length(s$data)))
  quit(status = 0)
}

if (mode == "compare") {
  w <- if (length(args) >= 2) args[2] else "none"
  wl <- whitelist_for(w)
  gold <- readRDS(GOLD)$values
  s <- snapshot(); new <- s$values
  if (length(s$errors)) { cat("ERRORS:\n", paste(s$errors, collapse = "\n"), "\n"); quit(status = 1) }
  in_wl <- function(p) length(wl) > 0 && any(vapply(wl, grepl, logical(1), x = p))
  diffs <- character(0); moved_wl <- character(0); contract_fail <- character(0)
  for (p in names(gold)) {
    if (!p %in% names(new)) {
      msg <- sprintf("REMOVED %s", p)
      if (in_wl(p)) moved_wl <- c(moved_wl, msg) else diffs <- c(diffs, msg)
      next
    }
    r <- same(gold[[p]], new[[p]])
    if (is.null(r)) next
    if (in_wl(p)) moved_wl <- c(moved_wl, sprintf("%s: %s", p, r)) else diffs <- c(diffs, sprintf("%s: %s", p, r))
  }
  added <- setdiff(names(new), names(gold))
  D <- s$data
  # contract checks
  if (any(grepl("KM_diff_band_modtau", wl, fixed = TRUE))) {
    for (nm in names(D)) for (g in 0:1) {
      ap <- new[[paste0(nm, "/KM_diff_band_modtau/at_points")]]
      S <- new[[paste0(nm, "/KM_diff_band_modtau/surv", g)]]
      e <- max(abs(S - sf_at(D[[nm]], g, ap)))
      if (!(e <= 1e-10)) contract_fail <- c(contract_fail, sprintf("B1 %s arm %d: max|S - survfit| = %.3g", nm, g, e))
    }
  }
  if (any(grepl("quantile_results", wl, fixed = TRUE))) {
    for (nm in names(D)) {
      d <- D[[nm]]
      for (call in c("df_counting_default", "df_counting_weighted", "df_counting_MB", "df_counting_q25")) {
        qp <- if (call == "df_counting_q25") 0.25 else 0.5
        pre <- paste0(nm, "/", call, "/quantile_results/")
        qres <- data.frame(quantile = new[[paste0(pre, "quantile")]], lower = new[[paste0(pre, "lower")]], upper = new[[paste0(pre, "upper")]])
        for (g in 0:1) {
          row <- if (g == 0) 2 else 1
          ref <- sf_q(d, g, qp, weighted = call == "df_counting_weighted")
          got <- unlist(qres[row, ])
          moved <- !isTRUE(all.equal(unname(got), unname(c(
            gold[[paste0(pre, "quantile")]][row], gold[[paste0(pre, "lower")]][row], gold[[paste0(pre, "upper")]][row])), tolerance = 0))
          if (call == "df_counting_weighted") { ref <- ref[1]; got <- got[1] }  # weighted limits use the package SE (M held)
          if (moved && !cmpv(unname(got), unname(ref)))
            contract_fail <- c(contract_fail, sprintf("B3 %s %s arm %d: moved and differs from survfit (%s vs %s)", nm, call, g,
                                                      paste(signif(got, 6), collapse = "/"), paste(signif(ref, 6), collapse = "/")))
        }
      }
      for (g in 0:1) {
        pre <- paste0(nm, "/km_quantile_arm", g, "/")
        got <- c(new[[paste0(pre, "qhat")]], new[[paste0(pre, "lower")]], new[[paste0(pre, "upper")]])
        old <- c(gold[[paste0(pre, "qhat")]], gold[[paste0(pre, "lower")]], gold[[paste0(pre, "upper")]])
        ref <- sf_q(d, g, 0.5)
        if (!cmpv(got, old) && !cmpv(got, unname(ref)))
          contract_fail <- c(contract_fail, sprintf("B3 %s km_quantile arm %d: moved and differs from survfit", nm, g))
      }
    }
  }
  ratio_lines <- character(0)
  if (any(grepl("wlr_dhat", wl, fixed = TRUE))) {
    for (nm in names(D)) for (f in c("sig2_dhat", "cov_wlr_dhat", "cor_wlr_dhat")) {
      p <- paste0(nm, "/wlr_dhat/", f)
      ratio_lines <- c(ratio_lines, sprintf("  %-13s %-13s old %.6g new %.6g old/new %.6f", nm, f, gold[[p]], new[[p]], gold[[p]] / new[[p]]))
    }
  }
  cat(sprintf("golden compare, whitelist '%s': %d fields in golden, %d in new run, %d added\n", w, length(gold), length(new), length(added)))
  if (length(added)) cat("ADDED (allowed):\n", paste0("  ", added, collapse = "\n"), "\n", sep = "")
  if (length(moved_wl)) cat(sprintf("WHITELISTED fields that moved (%d):\n", length(moved_wl)), paste0("  ", moved_wl, collapse = "\n"), "\n", sep = "")
  if (length(ratio_lines)) cat("B5 old/new ratios:\n", paste(ratio_lines, collapse = "\n"), "\n", sep = "")
  if (length(contract_fail)) cat("CONTRACT FAILURES:\n", paste0("  ", contract_fail, collapse = "\n"), "\n", sep = "")
  if (length(diffs)) cat(sprintf("DIFFERENCES outside the whitelist (%d):\n", length(diffs)), paste0("  ", diffs, collapse = "\n"), "\n", sep = "")
  ok <- !length(diffs) && !length(contract_fail)
  cat(if (ok) "GOLDEN: identical outside the whitelist\n" else "GOLDEN: FAIL\n")
  quit(status = if (ok) 0 else 1)
}
stop("usage: golden.R make | compare <whitelist>")
