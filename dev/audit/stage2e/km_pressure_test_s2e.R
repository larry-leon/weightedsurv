#!/usr/bin/env Rscript
# =============================================================================
# km_pressure_test_s2e.R -- copy of dev/audit/stage2d/km_pressure_test_s2d.R for Stage 2e, Step 5.
# Only difference: the default output directory is dev/audit/stage2e. Tests are unchanged.
#
# Stage-2d header:
# km_pressure_test_s2d.R -- copy of dev/audit/stage2c/km_pressure_test_s2c.R for Stage 2d, Step 5.
# Only difference: the default output directory is dev/audit/stage2d. Tests are unchanged.
#
# Stage-2 header:
# km_pressure_test_s2.R -- copy of dev/audit/km_pressure_test.R (Stage 1) for Stage 2, Step 5.
# Differences from the Stage-1 script: an output directory (--out=<dir>, default dev/audit/stage2)
# and an optional family filter (--only=T4,T5) for development runs. Tests are unchanged.
#
# Original header:
# km_pressure_test.R -- Kaplan-Meier pressure test (TASK_km_audit.md, section 3, T1-T10)
#
# Run from the repo root:
#     Rscript dev/audit/km_pressure_test.R            # writes the CSVs below
#     Rscript dev/audit/km_pressure_test.R --no-write # prints the summary only
#
# Writes  dev/audit/km_pressure_results.csv  (family, cell, n_checks, n_fail, worst, note)
#         dev/audit/km_scaling.csv           (T9 timings; differ run to run)
#
# Loads the LOCAL tree via pkgload::load_all (exports only).  survival is NOT
# attached at the start (T7 runs once without it, then attaches it).
# Failures of the package are data: the script exits 0 when tests fail and
# exits non-zero only if the script itself breaks (a SCRIPT_ERROR row).
# Changes nothing under R/.  All randomness is seeded.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
NO_WRITE <- "--no-write" %in% args
if (!file.exists("DESCRIPTION")) stop("Run from the repo root.")
suppressMessages(pkgload::load_all(".", quiet = TRUE, export_all = FALSE, helpers = FALSE))
out_arg <- sub("^--out=", "", grep("^--out=", args, value = TRUE))
out_dir <- if (length(out_arg)) out_arg[1] else file.path("dev", "audit", "stage2e")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
only_arg <- sub("^--only=", "", grep("^--only=", args, value = TRUE))
ONLY <- if (length(only_arg)) strsplit(only_arg[1], ",")[[1]] else character(0)
grDevices::pdf(NULL)
t_start <- proc.time()[["elapsed"]]
TOL <- 1e-10
QUICK <- nzchar(Sys.getenv("KMPT_QUICK"))   # development smoke run only: reduced T1 grid and T9 sizes

# ---- infrastructure -----------------------------------------------------------
RES <- list()
add_row <- function(family, cell, n_checks, n_fail, worst, note = "") {
  RES[[length(RES) + 1]] <<- data.frame(family = family, cell = cell, n_checks = n_checks,
                                        n_fail = n_fail, worst = as.character(worst),
                                        note = note, stringsAsFactors = FALSE)
}
SCRIPT_ERRORS <- character(0)
family <- function(name, code) {
  if (length(ONLY) && !name %in% ONLY) return(invisible(NULL))
  t0 <- proc.time()[["elapsed"]]
  r <- tryCatch({ code(); NULL }, error = function(e) conditionMessage(e))
  if (!is.null(r)) {
    SCRIPT_ERRORS <<- c(SCRIPT_ERRORS, paste(name, r))
    add_row(name, "SCRIPT_ERROR", NA, NA, NA, r)
  }
  message(sprintf("[%s] done in %.0fs", name, proc.time()[["elapsed"]] - t0))
}
run <- function(expr) {
  ws <- character(0); ms <- character(0)
  val <- tryCatch(withCallingHandlers(expr,
           warning = function(w) { ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning") },
           message = function(m) { ms <<- c(ms, conditionMessage(m)); invokeRestart("muffleMessage") }),
         error = function(e) structure(conditionMessage(e), class = "run_error"))
  err <- if (inherits(val, "run_error")) unclass(val) else NULL
  list(value = if (is.null(err)) val else NULL, error = err, warnings = ws, messages = ms)
}
fmt <- function(x, d = 3) if (length(x) == 0 || all(is.na(x))) "NA" else formatC(x, digits = d, format = "g")
mx <- function(x) { x <- x[is.finite(x)]; if (length(x)) max(x) else 0 }
one_line <- function(s) gsub("[\r\n]+", " ", s)

# ---- data and references ------------------------------------------------------
# two-arm data: control rate 0.08, treatment rate 0.08*exp(-0.3); exponential censoring
gen <- function(n_arm, cens = 0.3, ties = "none", wtype = "none", last = "asis", seed = 1) {
  set.seed(seed)
  N <- 2 * n_arm
  z <- rep(0:1, each = n_arm)
  rate <- 0.08 * exp(-0.3 * z)
  Tt <- stats::rexp(N, rate)
  C <- if (cens == 0) rep(Inf, N) else stats::rexp(N, rate * cens / (1 - cens))
  time <- pmin(Tt, C); status <- as.numeric(Tt <= C)
  time <- switch(ties, none = time, moderate = ceiling(time * 10) / 10, heavy = pmin(ceiling(time / 2.5), 8))
  u <- stats::runif(N, 0.3, 3)
  w <- switch(wtype,
    none = rep(1, N),
    int = sample(1:3, N, replace = TRUE),
    unif = u,
    ipw = 1 / stats::rbeta(N, 2, 8),
    norm = u / sum(u),
    zero10 = { k <- round(0.1 * N); if (k > 0) u[sample.int(N, k)] <- 0; u })
  for (g in 0:1) {
    i <- which(z == g); j <- i[which.max(time[i])]
    if (last == "event") status[j] <- 1
    if (last == "cens") status[j] <- 0
  }
  data.frame(time = time, status = status, treat = z, w = w)
}
fix_times <- function(d) unname(survival::aeqSurv(survival::Surv(d$time, d$status))[, "time"])
sf <- function(time, status, w = NULL, robust = FALSE) {
  y <- survival::Surv(time, status)
  if (is.null(w)) survival::survfit(y ~ 1) else survival::survfit(y ~ 1, weights = w, robust = robust)
}
# step-function lookup of a survfit object (equivalent to summary(fit, times, extend = TRUE); checked in REF)
step_at <- function(fit, at) {
  idx <- findInterval(at, fit$time)
  se <- if (isTRUE(fit$logse)) fit$std.err * fit$surv else fit$std.err   # robust fits store se on the S scale
  list(S = c(1, fit$surv)[idx + 1], se = c(0, se)[idx + 1])
}
# documented weighted SE (vignettes/articles/weightedkm_methodology.Rmd):
#   SE(t) = S(t) * sqrt( sum_{s<=t} [sum_i w_i^2 dN_i(s)] / Y_w(s)^2 )
doc_se <- function(time, status, w, at) {
  ut <- sort(unique(time[status == 1]))
  if (!length(ut)) return(list(S = rep(1, length(at)), se = rep(0, length(at))))
  Yw <- colSums(outer(time, ut, ">=") * w)
  E <- outer(time, ut, "==") * (status == 1)
  dNw <- colSums(E * w); dN2 <- colSums(E * w^2)
  S <- cumprod(ifelse(Yw > 0, 1 - dNw / Yw, 1))
  v <- cumsum(ifelse(Yw > 0, dN2 / Yw^2, 0))
  idx <- findInterval(at, ut)
  Sa <- c(1, S)[idx + 1]
  list(S = Sa, se = Sa * sqrt(c(0, v)[idx + 1]))
}
arm_ref <- function(d, tf, g, weighted, at, robust = FALSE) {
  i <- d$treat == g
  if (!any(i)) return(NULL)
  fit <- sf(tf[i], d$status[i], if (weighted) d$w[i] else NULL, robust = robust)
  step_at(fit, at)
}
err_class <- function(msg) {
  cryptic <- c("missing value where TRUE/FALSE needed", "subscript out of bounds", "object '.*' not found",
               "non-conformable", "differing number of rows", "replacement has", "dim\\(X\\) must",
               "argument is of length zero", "non-numeric argument", "undefined columns", "invalid 'type'",
               "arguments imply", "NAs are not allowed", "need finite", "incorrect number of dimensions",
               "\\$ operator is invalid", "must have a positive length", "could not find function",
               "'x' must be numeric", "invalid argument", "length of 'dimnames'", "is not TRUE",
               "'to' must be a finite number", "values must be length 1", "foreign function call",
               "first argument must be a vector", "system is exactly singular", "must be an array of at least two dimensions")
  if (any(vapply(cryptic, grepl, logical(1), x = msg))) "cryptic error" else "informative error"
}

# ==============================================================================
# REF -- sanity check of the reference machinery itself
# ==============================================================================
family("REF", function() {
  worst <- 0; n <- 0
  for (s in 1:6) {
    d <- gen(60, 0.3, c("none", "moderate", "heavy")[(s - 1) %% 3 + 1], "unif", seed = 900 + s)
    tf <- fix_times(d); at <- sort(c(d$time, max(d$time) + 5))
    for (rob in c(FALSE, TRUE)) {
      fit <- sf(tf[d$treat == 0], d$status[d$treat == 0], d$w[d$treat == 0], robust = rob)
      a <- step_at(fit, at); b <- summary(fit, times = unique(at), extend = TRUE)
      bb <- list(S = b$surv[match(at, unique(at))], se = b$std.err[match(at, unique(at))])
      ok <- a$S > 0
      worst <- max(worst, abs(a$S - bb$S), abs(a$se[ok] - bb$se[ok])); n <- n + 1
    }
  }
  add_row("REF", "step_at() vs summary(survfit, times, extend=TRUE)", n, as.integer(worst > TOL), fmt(worst),
          "reference lookup used throughout; robust and non-robust weighted fits")
})

# ==============================================================================
# T1 (+ T3 shape checks, + T6 check.KM warnings) -- agreement with survfit over a design grid
# ==============================================================================
NARM <- c(1, 2, 5, 20, 100, 500); CENS <- c(0, 0.3, 0.7, 0.95); TIES <- c("none", "moderate", "heavy")
WTS <- c("none", "int", "unif", "ipw", "norm", "zero10"); LAST <- c("asis", "event", "cens")
if (QUICK) NARM <- c(1, 5, 20)
cap <- 120000
REC <- list(check = character(cap), wtype = character(cap), ties = character(cap), where = character(cap),
            fail = logical(cap), value = numeric(cap), vmin = numeric(cap), vmax = numeric(cap), msg = character(cap))
nrec <- 0
CTX <- NULL
rec <- function(check, fail, value = NA_real_, vmin = NA_real_, vmax = NA_real_, msg = "") {
  nrec <<- nrec + 1
  REC$check[nrec] <<- check; REC$wtype[nrec] <<- CTX$wtype; REC$ties[nrec] <<- CTX$ties
  REC$where[nrec] <<- sprintf("n=%g,cens=%g,last=%s,seed=%d", CTX$n, CTX$cens, CTX$last, CTX$seed)
  REC$fail[nrec] <<- fail; REC$value[nrec] <<- value; REC$vmin[nrec] <<- vmin; REC$vmax[nrec] <<- vmax
  REC$msg[nrec] <<- msg
}
# T3 shape checks on one curve. 'at' grid, S curve, arm data (fixed times tt, status st, weights ww)
shape_checks <- function(src, at, S, tt, st, ww) {
  if (!length(tt)) return(invisible())
  rec(paste0("T3 in [0,1] (", src, ")"), any(S < 0 | S > 1, na.rm = TRUE) || anyNA(S),
      value = max(c(0, -S, S - 1), na.rm = TRUE))
  rec(paste0("T3 non-increasing (", src, ")"), any(diff(S) > 0, na.rm = TRUE), value = mx(c(0, diff(S))))
  fe <- suppressWarnings(min(tt[st == 1 & ww > 0]))
  pre <- at < fe
  rec(paste0("T3 S==1 before first event (", src, ")"), any(S[pre] != 1), value = mx(abs(S[pre] - 1)))
  L <- max(tt); post <- at >= L
  if (sum(post) > 1) rec(paste0("T3 carried forward after last obs (", src, ")"),
                         any(S[post] != S[post][1]), value = mx(abs(S[post] - S[post][1])))
  atL <- tt == L
  if (all(st[atL] == 1) && all(ww[atL] > 0) && any(post)) {
    v <- S[post]
    rec(paste0("T3 exactly 0 when arm ends on event (", src, ")"), any(v != 0), value = mx(abs(v)),
        msg = if (any(v < 0)) "negative" else "")
  }
}

family("T1", function() {
  cell_id <- 0
  for (n_arm in NARM) for (cens in CENS) for (ties in TIES) for (wtype in WTS) for (last in LAST) {
    cell_id <- cell_id + 1
    for (s in if (QUICK) 1 else 1:3) {
      seed <- 100000 + cell_id * 10 + s
      CTX <<- list(n = n_arm, cens = cens, ties = ties, wtype = wtype, last = last, seed = seed)
      d <- gen(n_arm, cens, ties, wtype, last, seed)
      weighted <- wtype != "none"; wn <- if (weighted) "w" else NULL
      tf <- fix_times(d)
      kd <- run(KM_diff(d, "time", "status", "treat", weight.name = wn))
      dc <- run(df_counting(d, "time", "status", "treat", weight.name = wn, check.KM = FALSE, get.cox = FALSE))
      if (!is.null(kd$error)) rec("error: KM_diff", TRUE, msg = kd$error) else rec("error: KM_diff", FALSE)
      if (!is.null(dc$error)) rec("error: df_counting", TRUE, msg = dc$error) else rec("error: df_counting", FALSE)
      for (g in 0:1) {
        i <- d$treat == g
        refS <- run(sf(tf[i], d$status[i], if (weighted) d$w[i] else NULL))
        if (!is.null(refS$error)) { rec("error: survfit", TRUE, msg = refS$error); next }
        if (!is.null(kd$value)) {
          k <- kd$value; S <- if (g == 0) k$surv0 else k$surv1; v2 <- if (g == 0) k$sig2_surv0 else k$sig2_surv1
          r <- step_at(refS$value, k$at_points)
          e <- mx(abs(S - r$S)); rec("survival vs survfit (KM_diff)", e > TOL || anyNA(S), value = e)
          ok <- r$S > 0 & is.finite(r$se)
          if (!weighted) {
            e2 <- mx(abs(sqrt(v2[ok]) - r$se[ok])); rec("unweighted SE vs survfit (KM_diff)", e2 > TOL || anyNA(v2), value = e2)
          } else {
            rr <- run(sf(tf[i], d$status[i], d$w[i], robust = TRUE))
            if (is.null(rr$error)) {
              rb <- step_at(rr$value, k$at_points); ok2 <- ok & rb$se > 0 & v2 > 0
              q <- sqrt(v2[ok2]) / rb$se[ok2]
              if (length(q)) rec("weighted SE ratio pkg/survfit robust (report)", NA, stats::median(q), min(q), max(q))
            } else rec("error: survfit robust", TRUE, msg = rr$error)
            ds <- doc_se(tf[i], d$status[i], d$w[i], k$at_points); ok3 <- ok & ds$se > 0 & v2 > 0
            q <- sqrt(v2[ok3]) / ds$se[ok3]
            if (length(q)) rec("weighted SE ratio pkg/documented formula (report)", NA, stats::median(q), min(q), max(q))
          }
          shape_checks("KM_diff", k$at_points, S, tf[i], d$status[i], d$w[i])
        }
        if (!is.null(dc$value)) {
          dd <- dc$value; S <- if (g == 0) dd$surv0 else dd$surv1
          r <- step_at(refS$value, dd$at_points)
          e <- mx(abs(S - r$S)); rec("survival vs survfit (df_counting)", e > TOL || anyNA(S), value = e)
          shape_checks("df_counting", dd$at_points, S, tf[i], d$status[i], d$w[i])
        }
      }
      if (!is.null(kd$value)) {
        k <- kd$value
        rec("dhat == surv1 - surv0 (KM_diff)", !isTRUE(all(abs(k$dhat - (k$surv1 - k$surv0)) <= TOL)), value = mx(abs(k$dhat - (k$surv1 - k$surv0))))
        rec("sig2_dhat == sig2_surv0 + sig2_surv1", !isTRUE(all(abs(k$sig2_dhat - (k$sig2_surv0 + k$sig2_surv1)) <= TOL)),
            value = mx(abs(k$sig2_dhat - (k$sig2_surv0 + k$sig2_surv1))))
        # T3: carried forward beyond the data
        ke <- run(KM_diff(d, "time", "status", "treat", weight.name = wn, at_points = c(sort(d$time), max(d$time) + c(1, 10))))
        if (is.null(ke$error)) for (g in 0:1) {
          i <- d$treat == g; S <- if (g == 0) ke$value$surv0 else ke$value$surv1
          L <- max(tf[i]); post <- ke$value$at_points >= L
          rec("T3 carried forward beyond data (KM_diff)", any(S[post] != S[post][1]), value = mx(abs(S[post] - S[post][1])))
        }
      }
      if (!is.null(dc$value)) {
        dd <- dc$value
        for (nm in c("survP", "survG")) {
          S <- dd[[nm]]
          rec(paste0("T3 ", nm, " in [0,1]"), any(S < 0 | S > 1) || anyNA(S), value = max(c(0, -S, S - 1), na.rm = TRUE))
          rec(paste0("T3 ", nm, " non-increasing"), any(diff(S) > 0), value = mx(c(0, diff(S))))
        }
        rp <- dd$survP; ev <- tf == max(tf)
        if (all(d$status[ev] == 1) && all(d$w[ev] > 0))
          rec("T3 survP exactly 0 when data end on event", rp[length(rp)] != 0, value = abs(rp[length(rp)]),
              msg = if (rp[length(rp)] < 0) "negative" else "")
      }
      # T6: check.KM = TRUE warnings on clean T1 data (seed 1 of each cell, n per arm 5..100)
      if (s == 1 && n_arm %in% c(5, 20, 100)) {
        ck <- run(df_counting(d, "time", "status", "treat", weight.name = wn, check.KM = TRUE, get.cox = FALSE))
        bad <- grep("Discrepancy in KM|outside|non-increasing", ck$warnings, value = TRUE)
        rec("T6 check.KM=TRUE warnings on clean data", length(bad) > 0 || !is.null(ck$error), value = length(bad),
            msg = paste(c(ck$error, unique(bad)), collapse = " / "))
      }
    }
  }
  R1 <- as.data.frame(lapply(REC, `[`, seq_len(nrec)), stringsAsFactors = FALSE)
  assign("T1REC", R1, envir = .GlobalEnv)
  # aggregate: one row per check x weights x ties
  key <- paste(R1$check, R1$wtype, R1$ties, sep = "|")
  for (kk in unique(key)) {
    x <- R1[key == kk, ]
    fam <- if (grepl("^T3", x$check[1])) "T3" else if (grepl("^T6", x$check[1])) "T6" else "T1"
    cell <- sprintf("%s | w=%s | ties=%s", x$check[1], x$wtype[1], x$ties[1])
    if (all(is.na(x$fail))) {                     # report-only ratios
      j <- which.max(abs(log(pmax(c(x$vmin, x$vmax), 1e-300))))
      worst <- c(x$vmin, x$vmax)[j]
      add_row(fam, cell, nrow(x), NA, fmt(worst),
              sprintf("median of per-arm medians %s; range [%s, %s]; worst at %s", fmt(stats::median(x$value)),
                      fmt(min(x$vmin)), fmt(max(x$vmax)), x$where[(j - 1) %% nrow(x) + 1]))
    } else {
      nf <- sum(x$fail, na.rm = TRUE)
      j <- if (nf && all(is.na(x$value[x$fail]))) which(x$fail)[1] else if (nf) which(x$fail)[which.max(x$value[x$fail])] else which.max(x$value)
      worst <- if (grepl("^error", x$check[1])) nf else x$value[j]
      msgs <- unique(x$msg[x$fail & nzchar(x$msg)])
      add_row(fam, cell, nrow(x), nf, fmt(worst),
              paste0(if (nf) paste0(if (grepl("^error", x$check[1])) "first at " else "worst at ", x$where[j],
                                    if (grepl("^error", x$check[1])) paste0("; failing cells: ", paste(utils::head(unique(sub(",seed=.*", "", x$where[x$fail])), 6), collapse = " / ")) else "") else "",
                     if (length(msgs)) paste0("; msgs: ", one_line(paste(utils::head(msgs, 3), collapse = " || "))) else ""))
    }
  }
})

# ==============================================================================
# T2 -- invariants
# ==============================================================================
family("T2", function() {
  DS <- list(); k <- 0
  for (ties in TIES) for (s in 1:10) { k <- k + 1; DS[[k]] <- gen(100, 0.3, ties, "unif", "asis", seed = 200000 + k) }
  km <- function(d, wn = "w", ...) KM_diff(d, "time", "status", "treat", weight.name = wn, ...)
  dcq <- function(d, wn = "w") df_counting(d, "time", "status", "treat", weight.name = wn, check.KM = FALSE, get.cox = FALSE)
  # 1. permutation
  e <- 0; nf <- 0
  for (j in seq_along(DS)) {
    d <- DS[[j]]; set.seed(j); p <- d[sample(nrow(d)), ]
    a <- km(d); b <- km(p); ad <- dcq(d); bd <- dcq(p)
    ee <- max(abs(a$surv0 - b$surv0), abs(a$surv1 - b$surv1), abs(a$sig2_surv0 - b$sig2_surv0), abs(a$sig2_surv1 - b$sig2_surv1),
              abs(ad$z.score - bd$z.score), abs(as.matrix(ad$quantile_results[, 2:4]) - as.matrix(bd$quantile_results[, 2:4])), na.rm = TRUE)
    e <- max(e, ee); nf <- nf + (ee > TOL)
  }
  add_row("T2", "1 row permutation: S, SE, quantiles, z (weighted)", length(DS), nf, fmt(e))
  # 2. weight rescaling
  for (cc in c("0.01", "10", "1/sum(w)")) {
    eS <- 0; nfS <- 0; rat <- c()
    for (d in DS) {
      c0 <- if (cc == "1/sum(w)") 1 / sum(d$w) else as.numeric(cc)
      d2 <- d; d2$w <- d$w * c0
      a <- km(d); b <- km(d2)
      ee <- max(abs(a$surv0 - b$surv0), abs(a$surv1 - b$surv1)); eS <- max(eS, ee); nfS <- nfS + (ee > TOL)
      ok <- a$sig2_surv0 > 0 & b$sig2_surv0 > 0
      rat <- c(rat, sqrt(b$sig2_surv0[ok] / a$sig2_surv0[ok]))
    }
    add_row("T2", sprintf("2 weights w -> c*w, c=%s: survival", cc), length(DS), nfS, fmt(eS))
    add_row("T2", sprintf("2 weights w -> c*w, c=%s: SE ratio (K3)", cc), length(DS), sum(abs(rat - 1) > 1e-8 ) > 0,
            fmt(rat[which.max(abs(log(rat)))]),
            sprintf("SE(c*w)/SE(w) over arm-0 curve points: quartiles %s", paste(fmt(stats::quantile(rat, c(.25, .5, .75))), collapse = ", ")))
  }
  # 3. time rescaling
  for (cc in c(1 / 12, 30.4375)) {
    eS <- 0; nfS <- 0; eq <- 0; nfq <- 0; eK <- 0; nfK <- 0
    for (d in DS) {
      d2 <- d; d2$time <- d$time * cc
      a <- km(d); b <- km(d2)
      ee <- max(abs(a$surv0 - b$surv0), abs(a$surv1 - b$surv1)); eS <- max(eS, ee); nfS <- nfS + (ee > TOL)
      qa <- dcq(d)$quantile_results$quantile; qb <- dcq(d2)$quantile_results$quantile
      eq0 <- max(abs(qb - cc * qa) / (cc * qa), na.rm = TRUE); eq <- max(eq, eq0); nfq <- nfq + (eq0 > 1e-10)
      kt <- km(d2, draws.band = 50, modify_tau = TRUE, show_resamples = FALSE); tf <- fix_times(d2)
      ek <- max(abs(kt$surv0 - step_at(sf(tf[d2$treat == 0], d2$status[d2$treat == 0], d2$w[d2$treat == 0]), kt$at_points)$S),
                abs(kt$surv1 - step_at(sf(tf[d2$treat == 1], d2$status[d2$treat == 1], d2$w[d2$treat == 1]), kt$at_points)$S))
      eK <- max(eK, ek); nfK <- nfK + (ek > 1e-6)
    }
    add_row("T2", sprintf("3 time t -> c*t, c=%s: survival", fmt(cc, 5)), length(DS), nfS, fmt(eS))
    add_row("T2", sprintf("3 time t -> c*t, c=%s: medians scale by c (rel. err)", fmt(cc, 5)), length(DS), nfq, fmt(eq))
    add_row("T2", sprintf("3 time t -> c*t, c=%s: KM_diff(modify_tau=TRUE) vs survfit (K1)", fmt(cc, 5)), length(DS), nfK, fmt(eK),
            "max |S - survfit| on the modify_tau grid (tau.seq = 0.25 in the rescaled units)")
  }
  # 4. duplicated rows == integer weight 2
  eS <- 0; nfS <- 0; eSE <- 0
  for (d in DS) {
    dd <- rbind(d, d); dd$w <- 1; d2 <- d; d2$w <- 2; ap <- sort(unique(d$time))
    a <- km(dd, at_points = ap); b <- km(d2, at_points = ap)
    ee <- max(abs(a$surv0 - b$surv0), abs(a$surv1 - b$surv1)); eS <- max(eS, ee); nfS <- nfS + (ee > TOL)
    eSE <- max(eSE, abs(sqrt(a$sig2_surv0) - sqrt(b$sig2_surv0)))
  }
  add_row("T2", "4 duplicated rows == weight 2: survival", length(DS), nfS, fmt(eS),
          sprintf("SE difference (report only): max %s (all-equal weights take the Greenwood branch)", fmt(eSE)))
  # 5. swap arms
  e <- 0; nf <- 0
  for (d in DS) {
    d2 <- d; d2$treat <- 1 - d$treat
    a <- km(d); b <- km(d2); ad <- dcq(d); bd <- dcq(d2)
    ee <- max(abs(a$dhat + b$dhat), abs(ad$z.score + bd$z.score),
              abs(ad$quantile_results$quantile - rev(bd$quantile_results$quantile)),
              abs(ad$quantile_results$lower - rev(bd$quantile_results$lower)), na.rm = TRUE)
    e <- max(e, ee); nf <- nf + (ee > TOL)
  }
  add_row("T2", "5 swap arm labels: dhat, z, quantile rows", length(DS), nf, fmt(e))
  # 6. add non-event times; 7. remove event times
  e6 <- 0; nf6 <- 0; e7 <- 0; nf7 <- 0
  for (j in seq_along(DS)) {
    d <- DS[[j]]; ap <- sort(unique(d$time)); base <- km(d, at_points = ap)
    mids <- (utils::head(ap, -1) + utils::tail(ap, -1)) / 2
    mids <- mids[!mids %in% d$time]
    ad <- km(d, at_points = sort(c(ap, mids))); m <- match(ap, ad$at_points)
    ee <- max(abs(ad$surv0[m] - base$surv0), abs(ad$surv1[m] - base$surv1), abs(ad$sig2_surv0[m] - base$sig2_surv0))
    e6 <- max(e6, ee); nf6 <- nf6 + (ee > TOL)
    set.seed(j); evt <- unique(d$time[d$status == 1]); drop <- evt[stats::runif(length(evt)) < 0.5]
    keep <- ap[!ap %in% drop]; rm7 <- km(d, at_points = keep); m <- match(keep, ap)
    ee <- max(abs(rm7$surv0 - base$surv0[m]), abs(rm7$surv1 - base$surv1[m])); e7 <- max(e7, ee); nf7 <- nf7 + (ee > TOL)
  }
  add_row("T2", "6 add non-event times to at_points: S, SE unchanged", length(DS), nf6, fmt(e6))
  add_row("T2", "7 remove ~50% of event times from at_points: S unchanged at kept points (K1)", length(DS), nf7, fmt(e7))
  # 8. draws = 0 vs draws = 2000 (weighted)
  rat <- c()
  for (d in DS[c(1, 11, 21)]) {
    a <- km(d); b <- km(d, draws = 2000)
    ok <- a$sig2_surv0 > 0 & a$surv0 > 0.1
    rat <- c(rat, sqrt(b$sig2_surv0[ok] / a$sig2_surv0[ok]))
  }
  add_row("T2", "8 SE ratio draws=2000 / draws=0 (weighted, report)", 3, NA, fmt(rat[which.max(abs(log(rat)))]),
          sprintf("quartiles %s (3 datasets, arm 0, S>0.1)", paste(fmt(stats::quantile(rat, c(.25, .5, .75))), collapse = ", ")))
})

# ==============================================================================
# T4 -- quantiles
# ==============================================================================
family("T4", function() {
  LV <- c(0.9, 0.75, 0.5, 0.25, 0.1)
  DS <- list(); k <- 0
  for (ties in TIES) for (cens in c(0, 0.3, 0.7)) for (s in 1:23) {
    k <- k + 1; set.seed(300000 + k); n_arm <- sample(c(10, 30, 60, 100), 1)
    DS[[k]] <- gen(n_arm, cens, ties, "none", "asis", seed = 300000 + k)
  }
  DS <- DS[1:200]
  res <- data.frame()
  for (j in seq_along(DS)) {
    d <- DS[[j]]; tf <- fix_times(d); i0 <- d$treat == 0
    kd <- KM_diff(d, "time", "status", "treat")
    for (cl in c(0.8, 0.9, 0.95)) {
      fit <- survival::survfit(survival::Surv(tf[i0], d$status[i0]) ~ 1, conf.int = cl, conf.type = "log")
      qq <- stats::quantile(fit, probs = 1 - LV)
      for (li in seq_along(LV)) {
        pk <- km_quantile(kd$at_points, kd$surv0, sqrt(kd$sig2_surv0), qprob = LV[li], conf_level = cl)
        ref <- c(unname(qq$quantile[li]), unname(qq$lower[li]), unname(qq$upper[li]))
        got <- c(pk$qhat, pk$lower, pk$upper)
        pl <- c(any(abs(fit$surv - LV[li]) < 1e-8), any(abs(fit$lower - LV[li]) < 1e-8, na.rm = TRUE), any(abs(fit$upper - LV[li]) < 1e-8, na.rm = TRUE))
        zc <- stats::qnorm(1 - (1 - cl) / 2); okc <- kd$surv0 > 0
        lo_c <- exp(log(kd$surv0[okc]) - zc * sqrt(kd$sig2_surv0[okc]) / kd$surv0[okc]); up_c <- exp(log(kd$surv0[okc]) + zc * sqrt(kd$sig2_surv0[okc]) / kd$surv0[okc])
        nm <- c(FALSE, any(diff(lo_c) > 1e-12), any(diff(up_c) > 1e-12))
        res <- rbind(res, data.frame(j = j, cl = cl, lv = LV[li], what = c("quantile", "lower", "upper"), got = got, ref = ref, plateau = pl, nonmono = nm))
      }
    }
  }
  res$bad <- with(res, (is.na(got) != is.na(ref)) | (!is.na(got) & !is.na(ref) & abs(got - ref) > 1e-10))
  res$err <- with(res, ifelse(!is.na(got) & !is.na(ref), abs(got - ref), ifelse(is.na(got) == is.na(ref), 0, Inf)))
  for (w in c("quantile", "lower", "upper")) for (lv in LV) {
    x <- res[res$what == w & res$lv == lv & (w == "quantile" & res$cl == 0.95 | w != "quantile"), ]
    nna <- sum(is.na(x$got) != is.na(x$ref))
    add_row("T4", sprintf("km_quantile %s, level S=%s vs quantile(survfit, 1-level)%s", w, lv, if (w != "quantile") ", conf 0.80/0.90/0.95" else ""),
            nrow(x), sum(x$bad), fmt(mx(x$err[is.finite(x$err)])),
            sprintf("NA-pattern mismatches %d; finite mismatches %d; mismatches where the survfit curve equals the level on an interval (K2 case): %d; mismatches where the package's log-scale limit curve is non-monotone: %d",
                    nna, sum(x$bad) - nna, sum(x$bad & x$plateau), sum(x$bad & x$nonmono)))
  }
  # df_counting()$quantile_results for each level (control arm) vs survfit
  nb <- 0; nn <- 0; ew <- 0; conv <- ""
  for (j in seq(1, 200, by = 4)) {
    d <- DS[[j]]; tf <- fix_times(d); i0 <- d$treat == 0
    fit <- survival::survfit(survival::Surv(tf[i0], d$status[i0]) ~ 1)
    for (lv in LV) {
      dq <- run(df_counting(d, "time", "status", "treat", qprob = lv, check.KM = FALSE, get.cox = FALSE))
      if (!is.null(dq$error)) { nb <- nb + 1; nn <- nn + 1; next }
      got <- dq$value$quantile_results$quantile[2]; ref <- unname(stats::quantile(fit, 1 - lv)$quantile)
      nn <- nn + 1
      bad <- (is.na(got) != is.na(ref)) || (!is.na(got) && abs(got - ref) > 1e-10)
      nb <- nb + bad; if (!is.na(got) && !is.na(ref)) ew <- max(ew, abs(got - ref))
    }
  }
  add_row("T4", "df_counting()$quantile_results (control) vs quantile(survfit, 1-qprob), 5 levels", nn, nb, fmt(ew))
  # convention (I5)
  d <- DS[[5]]; kd <- KM_diff(d, "time", "status", "treat"); tf <- fix_times(d); i0 <- d$treat == 0
  fit <- survival::survfit(survival::Surv(tf[i0], d$status[i0]) ~ 1)
  pk <- km_quantile(kd$at_points, kd$surv0, qprob = 0.25)$qhat
  conv <- sprintf("qprob=0.25 -> %s; survfit probs=0.75 -> %s; probs=0.25 -> %s", fmt(pk, 5),
                  fmt(unname(stats::quantile(fit, 0.75)$quantile), 5), fmt(unname(stats::quantile(fit, 0.25)$quantile), 5))
  add_row("T4", "qprob convention (I5)", 1, NA, "time at which S drops to qprob", conv)
  # plateau cases: no censoring, n = 2, 4, ..., 60, times 1..n
  for (lv in LV) {
    nb <- 0
    for (n in seq(2, 60, by = 2)) {
      km <- KM_estimates(n:1, 1:n)
      got <- km_quantile(1:n, km$S_KM, qprob = lv)$qhat
      ref <- unname(stats::quantile(survival::survfit(survival::Surv(1:n, rep(1, n)) ~ 1), 1 - lv)$quantile)
      nb <- nb + ((is.na(got) != is.na(ref)) || (!is.na(got) && abs(got - ref) > 1e-10))
    }
    add_row("T4", sprintf("exact-plateau cases n=2..60, level %s (K2)", lv), 30, nb, "", "package left end vs survfit midpoint when S equals the level on an interval")
  }
})

# ==============================================================================
# T5 -- KM_diff bands and resampling
# ==============================================================================
family("T5", function() {
  DS <- lapply(1:5, function(s) gen(80, 0.3, c("none", "moderate", "heavy", "none", "none")[s], c("none", "none", "none", "unif", "unif")[s], seed = 400000 + s))
  km <- function(d, ...) KM_diff(d, "time", "status", "treat", weight.name = if (length(unique(d$w)) > 1) "w" else NULL, show_resamples = FALSE, ...)
  nf <- 0; vals <- c()
  for (d in DS) {
    ca <- sapply(c(0.01, 0.05, 0.10, 0.20), function(a) unname(km(d, draws.band = 500, alpha = a)$c_alpha_band))
    vals <- rbind(vals, ca); nf <- nf + !all(diff(ca) < 0)
  }
  add_row("T5", "c_alpha_band strictly decreasing in alpha {0.01,0.05,0.10,0.20} (I1)", 5, nf, fmt(vals[1, 1], 5),
          sprintf("dataset 1 c_alpha by alpha: %s", paste(fmt(vals[1, ], 5), collapse = ", ")))
  nf <- 0; nfin <- 0
  for (d in DS) {
    k <- km(d, draws.band = 500); ok <- is.finite(k$lower) & is.finite(k$sb_lower)
    nf <- nf + !all(k$sb_lower[ok] <= k$lower[ok] + 1e-12 & k$sb_upper[ok] >= k$upper[ok] - 1e-12)
    nfin <- nfin + sum(!is.finite(k$sb_lower))
  }
  add_row("T5", "simultaneous band contains pointwise CI (alpha=0.05)", 5, nf, "", sprintf("non-finite band values: %d", nfin))
  d <- DS[[4]]
  a <- km(d, draws.band = 200, seedstart = 77); b <- km(d, draws.band = 200, seedstart = 77)
  add_row("T5", "fixed seedstart: two calls identical", 1, as.integer(!identical(a, b)), "")
  set.seed(2024); s0 <- .Random.seed; invisible(km(d, draws.band = 50)); same <- identical(s0, .Random.seed)
  add_row("T5", "caller's .Random.seed unchanged after KM_diff(draws.band=50) (P2)", 1, as.integer(!same), "",
          if (!same) ".Random.seed overwritten by set.seed(seedstart) inside KM_diff" else "")
  for (db in c(1, 5, 19, 20, 21)) {
    r <- run(KM_diff(d, "time", "status", "treat", weight.name = "w", draws.band = db))
    add_row("T5", sprintf("draws.band=%d, default show_resamples=TRUE (R5)", db), 1, as.integer(!is.null(r$error)),
            if (is.null(r$error)) "ok" else "error", one_line(if (is.null(r$error)) "" else r$error))
  }
  for (wt in c("none", "unif")) {
    d <- gen(100, 0.3, "none", wt, "asis", seed = 410000 + (wt == "unif"))
    k <- km(d, draws.band = 5000)
    ok <- k$surv0 > 0.05 & k$sig2_surv0 > 0
    rat <- apply(k$surv0_star, 1, stats::var)[ok] / k$sig2_surv0[ok]
    qs <- stats::quantile(seq_len(sum(ok)), c(.1, .25, .5, .75, .9), type = 1)
    add_row("T5", sprintf("var(surv0_star) over 5000 draws / sig2_surv0, weights=%s (report)", wt), sum(ok), NA,
            fmt(rat[which.max(abs(log(rat)))]),
            sprintf("ratio at 10/25/50/75/90%% of curve: %s", paste(fmt(rat[qs]), collapse = ", ")))
  }
})

# ==============================================================================
# T6 -- internal consistency of df_counting()
# ==============================================================================
family("T6", function() {
  e <- c(s = 0, rp = 0, P = 0, G = 0, z = 0, lr = 0); nf <- c(s = 0, rp = 0, P = 0, G = 0, z = 0, lr = 0); n <- c(s = 0, rp = 0, P = 0, G = 0, z = 0, lr = 0)
  for (j in 1:40) {
    wt <- if (j %% 2) "none" else "unif"; ties <- TIES[(j %% 3) + 1]
    d <- gen(50, 0.3, ties, wt, "asis", seed = 500000 + j); wn <- if (wt == "none") NULL else "w"
    tf <- fix_times(d); W <- if (is.null(wn)) NULL else d$w
    dc <- df_counting(d, "time", "status", "treat", weight.name = wn, check.KM = FALSE, get.cox = FALSE)
    kd <- KM_diff(d, "time", "status", "treat", weight.name = wn)
    m <- match(kd$at_points, dc$at_points); ee <- max(abs(dc$surv0[m] - kd$surv0), abs(dc$surv1[m] - kd$surv1))
    e["s"] <- max(e["s"], ee); nf["s"] <- nf["s"] + (ee > TOL); n["s"] <- n["s"] + 1
    for (g in 0:1) {
      i <- d$treat == g
      fit <- if (is.null(W)) survival::survfit(survival::Surv(tf[i], d$status[i]) ~ 1) else survival::survfit(survival::Surv(tf[i], d$status[i]) ~ 1, weights = W[i])
      nr <- summary(fit, times = dc$risk.points, extend = TRUE)$n.risk
      got <- if (g == 0) dc$riskpoints0 else dc$riskpoints1
      ee <- max(abs(got - nr)); e["rp"] <- max(e["rp"], ee); nf["rp"] <- nf["rp"] + (ee > TOL); n["rp"] <- n["rp"] + 1
    }
    fP <- if (is.null(W)) survival::survfit(survival::Surv(tf, d$status) ~ 1) else survival::survfit(survival::Surv(tf, d$status) ~ 1, weights = W)
    fG <- if (is.null(W)) survival::survfit(survival::Surv(tf, 1 - d$status) ~ 1) else survival::survfit(survival::Surv(tf, 1 - d$status) ~ 1, weights = W)
    ee <- mx(abs(dc$survP - step_at(fP, dc$at_points)$S)); e["P"] <- max(e["P"], ee); nf["P"] <- nf["P"] + (ee > TOL); n["P"] <- n["P"] + 1
    ee <- mx(abs(dc$survG - step_at(fG, dc$at_points)$S)); e["G"] <- max(e["G"], ee); nf["G"] <- nf["G"] + (ee > TOL); n["G"] <- n["G"] + 1
    if (is.null(wn)) {
      cs <- survival::survdiff(survival::Surv(d$time, d$status) ~ d$treat, rho = 0)$chisq
      ee <- abs(dc$z.score^2 - cs); e["z"] <- max(e["z"], ee); nf["z"] <- nf["z"] + (ee > 1e-8); n["z"] <- n["z"] + 1
      ee <- abs(dc$lr^2 / dc$sig2_lr - cs); e["lr"] <- max(e["lr"], ee); nf["lr"] <- nf["lr"] + (ee > 1e-8); n["lr"] <- n["lr"] + 1
    }
  }
  lab <- c(s = "surv0/surv1 equal KM_diff at shared points", rp = "riskpoints0/1 vs summary(survfit, risk.points)$n.risk (unweighted + weighted)",
           P = "survP vs pooled survfit", G = "survG vs survfit on 1 - status", z = "z.score^2 vs survdiff(rho=0)$chisq (unweighted)",
           lr = "lr^2/sig2_lr vs survdiff(rho=0)$chisq (unweighted)")
  for (nm in names(lab)) add_row("T6", lab[[nm]], n[[nm]], nf[[nm]], fmt(e[[nm]]), if (nm %in% c("z", "lr")) "tolerance 1e-8 (chisq from survdiff)" else "")
})

# ==============================================================================
# T7 -- stratified output (clean session first, then with survival attached)
# ==============================================================================
family("T7", function() {
  mk <- function(K, variant, seed) {
    d <- gen(150, 0.3, "moderate", "none", "asis", seed = seed)
    set.seed(seed); d$st <- paste0("s", sample.int(K, nrow(d), replace = TRUE))
    if (variant == "arm-empty") d <- d[!(d$st == "s1" & d$treat == 1), ]
    if (variant == "no-events") d$status[d$st == "s1"] <- 0
    d
  }
  for (attached in c(FALSE, TRUE)) {
    if (attached) suppressPackageStartupMessages(library(survival))
    for (K in c(2, 3, 6)) for (variant in c("plain", "arm-empty", "no-events")) {
      d <- mk(K, variant, 600000 + K)
      r <- run(df_counting(d, "time", "status", "treat", strata.name = "st", check.KM = FALSE))
      cell <- sprintf("K=%d strata, %s, survival attached=%s", K, variant, attached)
      if (!is.null(r$error)) { add_row("T7", cell, 1, 1, "error", one_line(r$error)); next }
      dc <- r$value; e <- 0; nchk <- 0
      tfx <- dc$time
      for (st in unique(dc$strata)) for (g in 0:1) {
        i <- dc$strata == st & dc$z == g
        nm <- if (g == 0) "surv0_stratum" else "surv1_stratum"
        ii <- dc$strata == st
        if (!any(i)) { e <- max(e, if (all(is.na(dc[[nm]][ii]))) 0 else NA_real_); next }
        fit <- survival::survfit(survival::Surv(tfx[i], dc$delta[i]) ~ 1)
        ref <- step_at(fit, tfx[ii])$S
        e <- max(e, abs(dc[[nm]][ii] - ref)); nchk <- nchk + 1
      }
      for (st in unique(dc$strata)) {
        ii <- dc$strata == st
        ref <- step_at(survival::survfit(survival::Surv(tfx[ii], dc$delta[ii]) ~ 1), tfx[ii])$S
        e <- max(e, abs(dc$survP_stratum[ii] - ref)); nchk <- nchk + 1
      }
      sh <- c(dc$surv0_stratum, dc$surv1_stratum, dc$survP_stratum)
      shape_bad <- any(sh < 0 | sh > 1, na.rm = TRUE)
      add_row("T7", cell, nchk, as.integer(is.na(e) || e > TOL || shape_bad), fmt(e),
              paste0(if (length(r$warnings)) paste0("warnings: ", one_line(paste(unique(r$warnings), collapse = " / "))) else "",
                     if (shape_bad) " ; *_stratum outside [0,1]" else ""))
    }
    if (attached) detach("package:survival", unload = FALSE)
  }
})

# ==============================================================================
# T8 -- input and error-path matrix
# ==============================================================================
family("T8", function() {
  set.seed(700001)
  base <- gen(60, 0.3, "none", "unif", "asis", seed = 700001)
  # reference check: does a result agree with survfit on the 'corrected' data?
  ok_result <- function(val, fn, ref_d, wn) {
    if (is.null(ref_d)) return(FALSE)
    fit <- switch(fn, KM_diff = val, df_counting = val, plotKM.band_subgroups = val$fit_itt)
    if (is.null(fit$surv0) && is.null(fit$surv1)) return(FALSE)
    tf <- fix_times(ref_d); e <- 0
    for (g in 0:1) {
      i <- ref_d$treat == g; S <- if (g == 0) fit$surv0 else fit$surv1
      if (!any(i)) { if (length(S) && !all(is.na(S))) return(FALSE); next }
      ref <- step_at(sf(tf[i], ref_d$status[i], if (is.null(wn)) NULL else ref_d$w[i]), fit$at_points)$S
      if (length(S) != length(ref) || anyNA(S)) return(FALSE)
      e <- max(e, abs(S - ref))
    }
    e < 1e-8
  }
  classify <- function(r, fn, ref_d, wn) {
    if (!is.null(r$error)) return(list(cls = err_class(r$error), msg = r$error))
    pm <- grep("^Error", r$messages, value = TRUE)
    if (length(pm)) return(list(cls = err_class(pm[1]), msg = paste("(safe_run message)", pm[1])))
    good <- isTRUE(tryCatch(ok_result(r$value, fn, ref_d, wn), error = function(e) FALSE))
    list(cls = if (good) "correct result" else "silent wrong result",
         msg = if (length(r$warnings)) paste("warnings:", paste(unique(r$warnings), collapse = " / ")) else "")
  }
  call_fn <- function(fn, d, tte = "time", ev = "status", tr = "treat", wn = "w", strata = NULL, at = NULL) {
    switch(fn,
      KM_diff = if (is.null(at)) KM_diff(d, tte, ev, tr, weight.name = wn) else KM_diff(d, tte, ev, tr, weight.name = wn, at_points = at),
      df_counting = df_counting(d, tte, ev, tr, weight.name = wn, strata.name = strata),   # package defaults (check.KM = TRUE, get.cox = TRUE)
      plotKM.band_subgroups = plotKM.band_subgroups(d, tte, ev, tr, weight.name = wn, draws.band = 20, draws = 0))
  }
  mod <- function(f) { d <- base; f(d) }
  ev1 <- base$time[base$status == 1]
  cases <- list(
    list("empty arm", mod(function(d) d[d$treat == 0, ]), function(d) d),
    list("arm with no events", mod(function(d) { d$status[d$treat == 1] <- 0; d }), function(d) d),
    list("all censored", mod(function(d) { d$status <- 0; d }), function(d) d),
    list("one observation", base[1, ], function(d) d),
    list("NA time", mod(function(d) { d$time[3] <- NA; d }), function(d) d[-3, ]),
    list("NA event", mod(function(d) { d$status[3] <- NA; d }), function(d) d[-3, ]),
    list("NA treatment", mod(function(d) { d$treat[3] <- NA; d }), function(d) d[-3, ]),
    list("NA weight", mod(function(d) { d$w[3] <- NA; d }), function(d) d[-3, ]),
    list("NA strata", mod(function(d) { d$st <- rep(c("a", "b"), 60); d$st[3] <- NA; d }), function(d) d[-3, ], strata = "st"),
    list("Inf time", mod(function(d) { d$time[3] <- Inf; d$status[3] <- 0; d }), function(d) d),
    list("zero time", mod(function(d) { d$time[3] <- 0; d }), function(d) d),
    list("negative time", mod(function(d) { d$time[3] <- -1; d }), NULL),
    list("zero weights", mod(function(d) { d$w[d$treat == 0][1:6] <- 0; d }), function(d) d[d$w > 0, ]),
    list("negative weights", mod(function(d) { d$w[1:5] <- -1; d }), NULL),
    list("treatment 1/2", mod(function(d) { d$treat <- d$treat + 1; d }), function(d) { d$treat <- d$treat - 1; d }),
    list("treatment factor", mod(function(d) { d$treat <- factor(d$treat); d }), function(d) { d$treat <- as.numeric(as.character(d$treat)); d }),
    list("treatment logical", mod(function(d) { d$treat <- d$treat == 1; d }), function(d) { d$treat <- as.numeric(d$treat); d }),
    list("treatment character", mod(function(d) { d$treat <- as.character(d$treat); d }), function(d) { d$treat <- as.numeric(d$treat); d }),
    list("event 1/2", mod(function(d) { d$status <- d$status + 1; d }), function(d) { d$status <- d$status - 1; d }),
    list("event logical", mod(function(d) { d$status <- d$status == 1; d }), function(d) { d$status <- as.numeric(d$status); d }),
    list("unsorted at_points", base, function(d) d, at = rev(sort(unique(base$time)))),
    list("duplicated at_points", base, function(d) d, at = sort(c(base$time, base$time[1:10]))),
    list("at_points beyond last time", base, function(d) d, at = sort(c(base$time, max(base$time) + c(1, 5)))),
    list("at_points not starting at/before first event", base, function(d) d, at = sort(unique(base$time))[sort(unique(base$time)) > sort(ev1)[3]]),
    list("missing column name", base, NULL, tte = "TIME"),
    list("tibble input", tibble::as_tibble(base), function(d) as.data.frame(d)),
    list("data.table input", if (requireNamespace("data.table", quietly = TRUE)) data.table::as.data.table(base) else NULL, function(d) as.data.frame(d))
  )
  for (cs in cases) {
    nm <- cs[[1]]; d <- cs[[2]]; reff <- cs[[3]]
    for (fn in c("KM_diff", "df_counting", "plotKM.band_subgroups")) {
      if (grepl("at_points", nm) && fn != "KM_diff") { add_row("T8", paste(nm, "|", fn), 0, NA, "n/a", "argument not available"); next }
      if (nm == "NA strata" && fn != "df_counting") { add_row("T8", paste(nm, "|", fn), 0, NA, "n/a", "no strata argument"); next }
      if (is.null(d)) { add_row("T8", paste(nm, "|", fn), 0, NA, "n/a", "package not installed"); next }
      r <- run(call_fn(fn, d, tte = if (is.null(cs$tte)) "time" else cs$tte, strata = cs$strata, at = cs$at))
      ref_d <- if (is.null(reff)) NULL else reff(as.data.frame(d))
      cl <- classify(r, fn, ref_d, "w")
      if (cl$cls == "silent wrong result" && !is.null(r$value) && !is.null(cs$at) && fn == "KM_diff") {
        # for at_points cases compare at the supplied points in the supplied order
        tf <- fix_times(base); e <- 0
        for (g in 0:1) { i <- base$treat == g; S <- if (g == 0) r$value$surv0 else r$value$surv1
          e <- max(e, abs(S - step_at(sf(tf[i], base$status[i], base$w[i]), cs$at)$S)) }
        if (e < 1e-8) cl$cls <- "correct result"
        cl$msg <- paste(cl$msg, sprintf("max|S - survfit| = %s", fmt(e)))
      }
      add_row("T8", paste(nm, "|", fn), 1, as.integer(cl$cls %in% c("cryptic error", "silent wrong result")), cl$cls, one_line(cl$msg))
    }
  }
  dc <- df_counting(base, "time", "status", "treat", check.KM = FALSE, get.cox = FALSE)
  for (tz in list(c("tzero before first time", -1), c("tzero after last time", 1e6))) {
    r <- run(wlr_dhat_estimates(dc, tzero = as.numeric(tz[2])))
    want <- if (as.numeric(tz[2]) < 0) 0 else utils::tail(dc$surv1, 1) - utils::tail(dc$surv0, 1)
    cls <- if (!is.null(r$error)) err_class(r$error) else if (length(r$value$dhat) == 1 && isTRUE(abs(r$value$dhat - want) < 1e-10)) "correct result" else "silent wrong result"
    add_row("T8", paste(tz[1], "| wlr_dhat_estimates"), 1, as.integer(cls %in% c("cryptic error", "silent wrong result")), cls,
            if (!is.null(r$error)) one_line(r$error) else sprintf("length(dhat)=%d, dhat=%s, expected %s", length(r$value$dhat), fmt(r$value$dhat), fmt(want)))
  }
})

# ==============================================================================
# T9 -- scaling (separate CSV; timings vary run to run)
# ==============================================================================
SCALING <- data.frame()
family("T9", function() {
  peak_mb <- function() { g <- gc(); sum(g[, ncol(g)]) }   # last column = 'max used (Mb)'
  for (N in if (QUICK) 500 else c(500, 2000, 4000, 8000)) {
    d <- gen(N / 2, 0.3, "none", "none", "asis", seed = 800000 + N)
    for (fn in c("df_counting", "KM_diff(draws.band=1000)", "cox_rhogamma(draws=0)")) {
      prev <- SCALING[SCALING$fn == fn, ]
      if (nrow(prev) && utils::tail(prev$peak_mb, 1) * (N / utils::tail(prev$n, 1))^2 > 6000) {
        SCALING <<- rbind(SCALING, data.frame(fn = fn, n = N, elapsed_s = NA, peak_mb = NA, note = "skipped: projected peak > 6 GB")); next
      }
      invisible(gc(reset = TRUE))
      tm <- system.time(r <- run(switch(fn,
        "df_counting" = df_counting(d, "time", "status", "treat"),
        "KM_diff(draws.band=1000)" = KM_diff(d, "time", "status", "treat", draws.band = 1000, show_resamples = FALSE),
        "cox_rhogamma(draws=0)" = cox_rhogamma(df_counting(d, "time", "status", "treat", check.KM = FALSE, get.cox = FALSE), draws = 0))))[["elapsed"]]
      pk <- peak_mb()
      SCALING <<- rbind(SCALING, data.frame(fn = fn, n = N, elapsed_s = round(tm, 2), peak_mb = round(pk), note = if (is.null(r$error)) "" else one_line(r$error)))
    }
  }
  add_row("T9", "scaling cells run (timings in km_scaling.csv)", nrow(SCALING), sum(nzchar(SCALING$note)), "",
          "cox_rhogamma timing includes its own df_counting(check.KM=FALSE, get.cox=FALSE) call")
})

# ==============================================================================
# T10 -- plot paths (pdf(NULL))
# ==============================================================================
family("T10", function() {
  drew <- 0; nul <- 0; msgs <- c()
  for (s in 1:60) {
    d <- gen(40, 0.3, "none", "unif", "event", seed = 900000 + s)
    dc <- run(df_counting(d, "time", "status", "treat", weight.name = "w", check.KM = FALSE))
    if (!is.null(dc$error)) { nul <- nul + 1; msgs <- c(msgs, dc$error); next }
    p <- run(plot_weighted_km(dc$value))
    em <- grep("^Error", p$messages, value = TRUE)
    if (length(em) || !is.null(p$error)) { nul <- nul + 1; msgs <- c(msgs, em, p$error) } else drew <- drew + 1
  }
  add_row("T10", "plot_weighted_km(), 60 weighted datasets ending on an event in each arm (K6)", 60, nul,
          sprintf("drew %d / returned NULL %d", drew, nul), one_line(paste(utils::head(unique(msgs), 2), collapse = " || ")))
  d <- gen(40, 0, "none", "none", "event", seed = 910001)
  dc <- df_counting(d, "time", "status", "treat", check.KM = FALSE)
  p <- run(KM_plot_2sample_weighted_counting(dc, conf.int = TRUE))
  z <- qnorm(0.975); lo <- exp(log(dc$surv0) - z * sqrt(dc$sig2_surv0) / dc$surv0); up <- exp(log(dc$surv0) + z * sqrt(dc$sig2_surv0) / dc$surv0)
  add_row("T10", "KM_plot_2sample_weighted_counting(conf.int=TRUE) when a curve reaches 0 (S5)", 1,
          as.integer(!is.null(p$error) || sum(is.nan(lo)) > 0 || any(up > 1, na.rm = TRUE)),
          if (is.null(p$error)) "drew" else "error",
          sprintf("min S0 = %s; NaN in log-CI limits: %d; upper limit > 1 at %d points (max %s); warnings: %s",
                  fmt(min(dc$surv0)), sum(is.nan(lo) | is.nan(up)), sum(up > 1, na.rm = TRUE), fmt(max(up, na.rm = TRUE)),
                  one_line(paste(unique(p$warnings), collapse = " / "))))
  base <- gen(60, 0.3, "none", "none", "asis", seed = 920001); base$sg <- rep(c("A", "B"), 60)
  base$treat[base$sg == "B"] <- 0
  r <- run(plotKM.band_subgroups(base, "time", "status", "treat", sg_labels = c("sg == 'A'", "sg == 'B'"), sg_colors = c("red", "blue"), draws = 0))
  add_row("T10", "plotKM.band_subgroups() with a subgroup whose treatment arm is empty (R1)", 1, as.integer(!is.null(r$error)),
          if (is.null(r$error)) "ok" else "error", one_line(if (is.null(r$error)) "" else r$error))
  base <- gen(60, 0.3, "none", "none", "asis", seed = 920002)
  r <- run(plotKM.band_subgroups(base, "time", "status", "treat", modify_tau = TRUE, draws = 0, draws.band = 50))
  e <- NA
  if (is.null(r$error)) {
    f <- r$value$fit_itt; tf <- fix_times(base)
    e <- max(abs(f$surv0 - step_at(sf(tf[base$treat == 0], base$status[base$treat == 0]), f$at_points)$S),
             abs(f$surv1 - step_at(sf(tf[base$treat == 1], base$status[base$treat == 1]), f$at_points)$S))
  }
  add_row("T10", "plotKM.band_subgroups(modify_tau=TRUE): fit_itt survival vs survfit (K1, S6)", 1, as.integer(is.na(e) || e > 1e-6), fmt(e),
          if (!is.null(r$error)) one_line(r$error) else sprintf("grid of %d points", length(r$value$fit_itt$at_points)))
  k0 <- KM_diff(base, "time", "status", "treat", draws.band = 0)
  r <- run(cumulative_rmst_bands(base, k0, "time", "status", "treat", plot = FALSE, draws_sb = 100))
  add_row("T10", "cumulative_rmst_bands() given KM_diff(draws.band=0) (S4)", 1, as.integer(!is.null(r$error)),
          if (is.null(r$error)) "ok" else err_class(r$error), one_line(if (is.null(r$error)) "" else r$error))
})

# ---- output -------------------------------------------------------------------
res <- do.call(rbind, RES)
if (!NO_WRITE) {
  utils::write.csv(res, file.path(out_dir, "km_pressure_results.csv"), row.names = FALSE)
  utils::write.csv(SCALING, file.path(out_dir, "km_scaling.csv"), row.names = FALSE)
}
cat("\nweightedsurv KM pressure test -- LOCAL TREE via pkgload::load_all\n")
cat("R:", R.version.string, "| survival", as.character(utils::packageVersion("survival")), "\n\n")
fam <- split(res, factor(res$family, levels = unique(res$family)))
cat(sprintf("%-5s %6s %8s %7s\n", "fam", "cells", "checks", "failed"))
for (f in names(fam)) {
  x <- fam[[f]]
  cat(sprintf("%-5s %6d %8d %7d\n", f, nrow(x), sum(x$n_checks, na.rm = TRUE), sum(x$n_fail, na.rm = TRUE)))
}
cat("\nT9 scaling:\n"); print(SCALING, row.names = FALSE)
cat(sprintf("\nelapsed %.0fs; %s\n", proc.time()[["elapsed"]] - t_start,
            if (NO_WRITE) "--no-write: no files written" else paste("written:", file.path(out_dir, c("km_pressure_results.csv", "km_scaling.csv")), collapse = ", ")))
if (length(SCRIPT_ERRORS)) { cat("SCRIPT ERRORS:\n", paste(SCRIPT_ERRORS, collapse = "\n"), "\n"); quit(status = 1) }
