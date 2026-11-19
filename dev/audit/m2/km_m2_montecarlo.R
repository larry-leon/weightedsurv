#!/usr/bin/env Rscript
# =============================================================================
# km_m2_montecarlo.R -- M2: weighted-SE evaluation when weights are related to the hazard
# (dev/tasks/TASK_km_m2_weighted_se.md). Evaluation only; nothing under R/ is touched.
#
#   Rscript dev/audit/m2/km_m2_montecarlo.R --sanity   G1: formula cross-checks on one cell-D dataset,
#                                                       200 replicates of cell A, timing projection
#                                                       -> dev/audit/m2/km_m2_sanity.txt (via stdout)
#   Rscript dev/audit/m2/km_m2_montecarlo.R --full     all 10 arm-cells x 2,000 replicates
#                                                       -> km_m2_results.csv, km_m2_replicates.rds
#   Rscript dev/audit/m2/km_m2_montecarlo.R --check    re-runs the first 20 replicates of cell B serially
#                                                       and compares them with km_m2_replicates.rds;
#                                                       exit 0 when identical (tolerance 1e-12)
#
# Design (task section 2): C ~ U(0, 40); t = 5, 10, 20; replicate r of cell k uses
# set.seed(100000 * k + r). Cells A, B, C, Bt: one group, n = 400, X ~ Bernoulli(0.5), w = 1 + 2X.
# Cells D, E, F: two arms, N = 800, X ~ N(0, 1), true / estimated / extreme propensity weights,
# each arm evaluated separately.
# Estimators per arm-cell: pkg (analytic, draws = 0), pkg_resampling (draws = 200,
# seedstart = 816951 + r; first 500 replicates only), robust (survfit robust = TRUE),
# ij (explicit infinitesimal jackknife), xie_liu.
# The package is loaded from the local tree (branch km-fixes-0.2.0) with pkgload::load_all.
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) sub("^--", "", args[1]) else "sanity"
stopifnot(mode %in% c("sanity", "full", "check"))
if (!file.exists("DESCRIPTION")) stop("Run from the repo root.")
suppressMessages(pkgload::load_all(".", quiet = TRUE, export_all = FALSE, helpers = FALSE))
OUT <- file.path("dev", "audit", "m2")
TT <- c(5, 10, 20)
R_FULL <- 2000; R_RESAMP <- 500
CORES <- max(1L, parallel::detectCores() - 2L)

# ---- cells ---------------------------------------------------------------------------
CELLS <- data.frame(k = 1:7, cell = c("A", "B", "C", "Bt", "D", "E", "F"), stringsAsFactors = FALSE)
expit <- function(x) 1 / (1 + exp(-x))

gen <- function(k, r) {
  set.seed(100000 * k + r)
  cell <- CELLS$cell[k]
  if (cell %in% c("A", "B", "C", "Bt")) {
    n <- 400
    X <- stats::rbinom(n, 1, 0.5)
    haz <- switch(cell, A = rep(0.10, n), B = , Bt = ifelse(X == 1, 0.15, 0.05), C = ifelse(X == 1, 0.05, 0.15))
    Tt <- stats::rexp(n, haz); C <- stats::runif(n, 0, 40)
    U <- pmin(Tt, C); D <- as.numeric(Tt <= C)
    if (cell == "Bt") U <- ceiling(U * 2) / 2
    w <- 1 + 2 * X
    return(list(list(name = cell, U = U, D = D, w = w)))
  }
  N <- 800
  X <- stats::rnorm(N)
  slope <- if (cell == "F") 2 else 1
  e <- expit(-0.5 + slope * X)
  A <- stats::rbinom(N, 1, e)
  Tt <- stats::rexp(N, 0.08 * exp(0.7 * X - 0.4 * A)); C <- stats::runif(N, 0, 40)
  U <- pmin(Tt, C); D <- as.numeric(Tt <= C)
  if (cell == "E") e <- stats::fitted(stats::glm(A ~ X, family = stats::binomial))
  w <- ifelse(A == 1, 1 / e, 1 / (1 - e))
  lapply(0:1, function(a) { i <- A == a; list(name = paste0(cell, a), U = U[i], D = D[i], w = w[i]) })
}

# ---- truth: E[w 1(T > t)] / E[w] ----------------------------------------------------------
truth_of <- function(name) {
  mix <- function(h0, h1) (0.5 * 1 * exp(-h0 * TT) + 0.5 * 3 * exp(-h1 * TT)) / 2  # w = 1 (X = 0), 3 (X = 1)
  switch(substr(name, 1, 1),
    A = exp(-0.10 * TT),
    B = mix(0.05, 0.15),            # B and Bt
    C = mix(0.15, 0.05),
    { a <- as.numeric(substr(name, 2, 2))
      vapply(TT, function(t) stats::integrate(function(x) exp(-0.08 * exp(0.7 * x - 0.4 * a) * t) * stats::dnorm(x),
                                              -Inf, Inf, rel.tol = 1e-12)$value, 0) })
}

# ---- estimators for one group ------------------------------------------------------------
at_t <- function(at, v, default) { i <- findInterval(TT, at); out <- rep(default, length(TT)); out[i > 0] <- v[i[i > 0]]; out }

group_stats <- function(U, D, w, r, resamp) {
  at <- sort(unique(U))
  re <- calculate_risk_event_counts(U, D, w, at)
  km <- KM_estimates(re$ybar, re$nbar, re$sig2w_multiplier)
  S <- at_t(at, km$S_KM, 1); se_pkg <- at_t(at, sqrt(km$sig2_KM), 0)
  se_rs <- rep(NA_real_, length(TT))
  if (resamp) {
    rr <- calculate_risk_event_counts(U, D, w, at, draws = 200, seedstart = 816951 + r)
    kr <- KM_estimates(rr$ybar, rr$nbar, rr$sig2w_multiplier)
    se_rs <- at_t(at, sqrt(kr$sig2_KM), 0)
  }
  fit <- survival::survfit(survival::Surv(U, D) ~ 1, weights = w, robust = TRUE)
  sm <- summary(fit, times = TT, extend = TRUE)
  # survfit's default timefix merges observed times closer than ~1e-8 (relative); the package does
  # not. The point-estimate check is also made against survfit(timefix = FALSE).
  fit_nf <- survival::survfit(survival::Surv(U, D) ~ 1, weights = w, timefix = FALSE)
  S_nf <- summary(fit_nf, times = TT, extend = TRUE)$surv
  n_merged <- length(unique(U)) - length(unique(survival::aeqSurv(survival::Surv(U, D))[, 1]))
  # explicit IJ and Xie-Liu on the event times
  s <- sort(unique(U[D == 1]))
  Yw <- vapply(s, function(x) sum(w[U >= x]), 0)
  dNw <- vapply(s, function(x) sum(w[U == x & D == 1]), 0)
  Y2 <- vapply(s, function(x) sum(w[U >= x]^2), 0)
  okd <- Yw > 0 & (Yw - dNw) > 1e-10 * Yw
  g <- ifelse(okd, 1 / (Yw - dNw), 0)
  h <- ifelse(okd, dNw / (Yw * (Yw - dNw)), 0)
  H <- cumsum(h)
  q <- 1 - dNw / Yw; M <- Yw^2 / Y2
  xl <- ifelse(q > 1e-10, (1 - q) / (M * q), 0); XL <- cumsum(xl)
  ij <- xie <- numeric(length(TT))
  for (j in seq_along(TT)) {
    t <- TT[j]
    first <- ifelse(D == 1 & U <= t, g[match(U, s)], 0)
    first[is.na(first)] <- 0
    k <- findInterval(pmin(t, U), s)
    second <- ifelse(k > 0, H[pmax(k, 1)], 0)
    psi <- first - second
    ij[j] <- S[j] * sqrt(sum(w^2 * psi^2))
    kt <- findInterval(t, s)
    xie[j] <- S[j] * sqrt(if (kt > 0) XL[kt] else 0)
  }
  data.frame(t = TT, S = S, survfit_S = sm$surv, survfit_nofix_S = S_nf, n_merged = n_merged, pkg = se_pkg, pkg_resampling = se_rs, robust = sm$std.err,
             ij = ij, xie_liu = xie, ess = sum(w)^2 / (length(w) * sum(w^2)))
}

one_rep <- function(k, r) {
  groups <- gen(k, r)
  do.call(rbind, lapply(groups, function(gr) cbind(armcell = gr$name, r = r,
                                                    group_stats(gr$U, gr$D, gr$w, r, resamp = r <= R_RESAMP))))
}
run_cell <- function(k, reps, cores = CORES) {
  res <- parallel::mclapply(reps, function(r) one_rep(k, r), mc.cores = cores, mc.preschedule = TRUE)
  bad <- vapply(res, inherits, logical(1), "try-error")
  if (any(bad)) stop("cell ", CELLS$cell[k], ": ", as.character(res[[which(bad)[1]]]))
  do.call(rbind, res)
}

EST <- c("pkg", "pkg_resampling", "robust", "ij", "xie_liu")
summarise <- function(rep) {
  out <- list()
  for (ac in unique(rep$armcell)) for (t in TT) {
    d <- rep[rep$armcell == ac & rep$t == t, ]
    tr <- truth_of(ac)[match(t, TT)]
    sdS <- stats::sd(d$S)
    for (e in EST) {
      se <- d[[e]]; use <- !is.na(se)
      fin <- use & is.finite(se)
      out[[length(out) + 1]] <- data.frame(armcell = ac, estimator = e, t = t, truth = tr, n_reps = sum(use),
        mean_se = mean(se[fin]), emp_sd = sdS, se_over_sd = mean(se[fin]) / sdS,
        coverage = mean(abs(d$S[fin] - tr) <= 1.96 * se[fin]), n_nonfinite = sum(use & !is.finite(se)),
        bias = mean(d$S) - tr, ess_ratio = mean(d$ess),
        max_abs_diff_survfit = max(abs(d$S - d$survfit_S)), max_abs_diff_survfit_nofix = max(abs(d$S - d$survfit_nofix_S)),
        reps_timefix_merged = length(unique(d$r[d$n_merged > 0])), stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, out)
}

# ---- cross-checks on one cell-D dataset (G1) ---------------------------------------------
cross_checks <- function() {
  set.seed(100000 * 5 + 1)
  N <- 800; X <- stats::rnorm(N); e <- expit(-0.5 + X); A <- stats::rbinom(N, 1, e)
  Tt <- stats::rexp(N, 0.08 * exp(0.7 * X - 0.4 * A)); C <- stats::runif(N, 0, 40)
  U <- pmin(Tt, C); D <- as.numeric(Tt <= C); w <- ifelse(A == 1, 1 / e, 1 / (1 - e))
  lines <- character(0)
  for (tied in c(FALSE, TRUE)) {
    UU <- if (tied) ceiling(U * 2) / 2 else U
    for (a in 0:1) {
      i <- A == a
      gs <- group_stats(UU[i], D[i], w[i], 1, resamp = FALSE)
      rd <- max(abs(gs$ij - gs$robust) / gs$robust)
      lines <- c(lines, sprintf("ij vs robust, cell D dataset r = 1, arm %d, %s: max relative difference %.3g", a,
                                if (tied) "times rounded up to 0.5 (ties)" else "no ties", rd))
      if (rd > 1e-6) stop("ij differs from robust by more than 1e-6: ", rd)
    }
  }
  ac <- NA
  if (requireNamespace("adjustedCurves", quietly = TRUE)) {
    dd <- data.frame(time = U, event = D, grp = factor(A, levels = 0:1))
    adj <- suppressMessages(suppressWarnings(adjustedCurves::adjustedsurv(data = dd, variable = "grp", ev_time = "time", event = "event",
                              method = "iptw_km", treatment_model = w, conf_int = TRUE)))$adj
    ac <- 0
    for (a in 0:1) {
      i <- A == a
      gs <- group_stats(U[i], D[i], w[i], 1, resamp = FALSE)
      sa <- adj[adj$group == as.character(a), ]; sa <- sa[order(sa$time), ]
      ref <- sa$se[findInterval(TT, sa$time)]
      rd <- max(abs(gs$xie_liu - ref) / ref)
      ac <- max(ac, rd)
      lines <- c(lines, sprintf("xie_liu vs adjustedCurves::adjustedsurv(iptw_km, conf_int = TRUE), arm %d: max relative difference %.3g (xie_liu %s; adjustedCurves %s)",
                                a, rd, paste(signif(gs$xie_liu, 6), collapse = "/"), paste(signif(ref, 6), collapse = "/")))
    }
  } else lines <- c(lines, "adjustedCurves is not installed: xie_liu cross-check not run")
  # point estimate equals survfit(weights = w) on this dataset
  for (a in 0:1) { i <- A == a; gs <- group_stats(U[i], D[i], w[i], 1, FALSE)
    lines <- c(lines, sprintf("point estimate vs survfit(weights = w), arm %d: max abs difference %.3g", a, max(abs(gs$S - gs$survfit_S)))) }
  list(lines = lines, adjustedCurves_maxrel = ac)
}

# ---- modes -------------------------------------------------------------------------------
if (mode == "sanity") {
  cat(sprintf("M2 sanity run (Gate G1); %d cores; HEAD %s\n", CORES, system("git rev-parse --short HEAD", intern = TRUE)))
  cc <- cross_checks(); cat(cc$lines, sep = "\n")
  t0 <- proc.time()[["elapsed"]]
  repA <- run_cell(1, 1:200)
  tA <- proc.time()[["elapsed"]] - t0
  sA <- summarise(repA)
  options(width = 200)
  print(sA[, c("armcell", "estimator", "t", "truth", "n_reps", "se_over_sd", "coverage", "n_nonfinite", "bias", "max_abs_diff_survfit", "max_abs_diff_survfit_nofix")], row.names = FALSE)
  fails <- character(0)
  if (any(abs(sA$bias) > 0.01)) fails <- c(fails, "absolute bias above 0.01")
  if (any(sA$max_abs_diff_survfit > 1e-10)) fails <- c(fails, "point estimate differs from survfit by more than 1e-10")
  if (any(sA$max_abs_diff_survfit_nofix > 1e-10)) fails <- c(fails, "point estimate differs from survfit(timefix = FALSE) by more than 1e-10")
  core <- sA[sA$estimator %in% c("pkg", "robust", "ij", "xie_liu"), ]
  if (any(core$se_over_sd < 0.80 | core$se_over_sd > 1.20)) fails <- c(fails, "SE/SD outside 0.80-1.20")
  # timing projection: 200 reps of A (with resampling), plus 20 reps of every other cell, with
  # and without resampling, scaled to 500 resampling + 1,500 plain replicates per cell
  proj <- tA / 200 * R_FULL   # A: all 200 had resampling, an upper bound for A
  tim <- c(A = tA / 200)
  for (k in 2:7) {
    t1 <- proc.time()[["elapsed"]]; invisible(run_cell(k, 1:20)); a1 <- (proc.time()[["elapsed"]] - t1) / 20
    t2 <- proc.time()[["elapsed"]]; invisible(run_cell(k, (R_RESAMP + 1):(R_RESAMP + 20))); a2 <- (proc.time()[["elapsed"]] - t2) / 20
    tim[CELLS$cell[k]] <- a1
    proj <- proj + a1 * R_RESAMP + a2 * (R_FULL - R_RESAMP)
  }
  cat(sprintf("\nwall-clock per replicate (s, %d cores, with resampling): %s\n", CORES, paste(sprintf("%s %.4f", names(tim), tim), collapse = ", ")))
  cat(sprintf("projected wall-clock for the full run: %.1f min\n", proj / 60))
  if (proj > 30 * 60) fails <- c(fails, "projected wall-clock above 30 min")
  cat(if (length(fails)) paste0("G1: FAIL (", paste(fails, collapse = "; "), ")\n") else "G1: PASS\n")
  quit(status = if (length(fails)) 1 else 0)
}

if (mode == "full") {
  cat(sprintf("M2 full run; %d cores; HEAD %s; started %s\n", CORES, system("git rev-parse --short HEAD", intern = TRUE), format(Sys.time())))
  cc <- cross_checks(); cat(cc$lines, sep = "\n")
  t0 <- proc.time()[["elapsed"]]
  reps <- list()
  for (k in CELLS$k) {
    tk <- proc.time()[["elapsed"]]
    reps[[k]] <- run_cell(k, seq_len(R_FULL))
    el <- proc.time()[["elapsed"]] - t0
    cat(sprintf("cell %s: %d replicates in %.0f s (total %.1f min)\n", CELLS$cell[k], R_FULL, proc.time()[["elapsed"]] - tk, el / 60))
    if (el > 30 * 60 && k < 7) { cat("STOP: wall-clock passed 30 min\n"); break }
  }
  rep <- do.call(rbind, reps)
  res <- summarise(rep)
  num <- c("truth", "mean_se", "emp_sd", "se_over_sd", "coverage", "bias", "ess_ratio")
  res[num] <- lapply(res[num], signif, 6)
  utils::write.csv(res, file.path(OUT, "km_m2_results.csv"), row.names = FALSE)
  saveRDS(list(replicates = rep, cross_checks = cc$lines), file.path(OUT, "km_m2_replicates.rds"), compress = "xz")
  options(width = 220)
  print(res[, c("armcell", "estimator", "t", "truth", "n_reps", "mean_se", "emp_sd", "se_over_sd", "coverage", "n_nonfinite", "bias", "ess_ratio")], row.names = FALSE)
  cat(sprintf("\nreplicates run: %s\n", paste(sprintf("%s %d", unique(rep$armcell), tapply(rep$r, rep$armcell, function(x) length(unique(x)))[unique(rep$armcell)]), collapse = ", ")))
  cat(sprintf("max |S - survfit(weights = w)| over all replicates and times: %.3g\n", max(abs(rep$S - rep$survfit_S))))
  cat(sprintf("max |S - survfit(weights = w, timefix = FALSE)| over all replicates and times: %.3g\n", max(abs(rep$S - rep$survfit_nofix_S))))
  mg <- rep[rep$t == TT[1], ]
  cat("replicates in which survfit's timefix merged observed times (per arm-cell):",
      paste(sprintf("%s %d", unique(mg$armcell), tapply(mg$n_merged > 0, mg$armcell, sum)[unique(mg$armcell)]), collapse = ", "), "\n")
  dd <- rep[abs(rep$S - rep$survfit_S) > 1e-10, ]
  cat(sprintf("arm-cell x replicate x time values with |S - survfit| > 1e-10: %d (all in replicates with merged times: %s); max |ij - robust| / robust there %.3g, elsewhere %.3g\n",
              nrow(dd), all(dd$n_merged > 0), if (nrow(dd)) max(abs(dd$ij - dd$robust) / dd$robust) else NA,
              max((abs(rep$ij - rep$robust) / rep$robust)[rep$n_merged == 0 & rep$robust > 0])))
  cat(sprintf("wall-clock: %.1f min\n", (proc.time()[["elapsed"]] - t0) / 60))
  quit(status = 0)
}

if (mode == "check") {
  stored <- readRDS(file.path(OUT, "km_m2_replicates.rds"))$replicates
  old <- stored[stored$armcell == "B" & stored$r <= 20, ]
  new <- do.call(rbind, lapply(1:20, function(r) one_rep(2, r)))
  old <- old[order(old$r, old$t), ]; new <- new[order(new$r, new$t), ]
  cols <- c("S", "survfit_S", "survfit_nofix_S", "n_merged", "pkg", "pkg_resampling", "robust", "ij", "xie_liu", "ess")
  d <- max(vapply(cols, function(cn) max(abs(old[[cn]] - new[[cn]]), na.rm = TRUE), 0))
  na_same <- identical(unname(is.na(as.matrix(old[cols]))), unname(is.na(as.matrix(new[cols]))))
  cat(sprintf("check: first 20 replicates of cell B re-run serially; max abs difference %.3g; NA pattern %s\n",
              d, if (na_same) "identical" else "DIFFERS"))
  quit(status = if (d <= 1e-12 && na_same && nrow(old) == nrow(new)) 0 else 1)
}
