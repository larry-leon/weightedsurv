#!/usr/bin/env Rscript
# =============================================================================
# km_m4_null.R -- M4: weighted log-rank and Cox under the null, measurement only
# (TASK_km_fixes_stage2c.md, section 7). Nothing under R/ follows from it.
#
#   Rscript dev/audit/stage2c/km_m4_null.R [reps]   -> km_m4_results.csv, km_m4_replicates.rds, table on stdout
#
# Cells U, D, E, F as in M3 (km_m3_bands.R) but with no treatment effect: T ~ Exp(0.08 exp(0.7 X)).
# Replicate r of cell k (k = 1..4 for U, D, E, F) uses set.seed(300000 * k + r).
# Per replicate: dc <- df_counting(df, ..., weight.name); cox_rhogamma(dc, scheme = "fh",
# scheme_params = list(rho = 0, gamma = 0), draws = 0); reference coxph(Surv(time, status) ~ treat,
# weights = w, robust = TRUE) (cell U: no weights). In cells D, E, F the two package calls are repeated
# with the weights multiplied by 10 and divided by their sum. Package loaded from HEAD (load_all).
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
R_REPS <- if (length(args)) as.integer(args[1]) else 2000L
if (!file.exists("DESCRIPTION")) stop("Run from the repo root.")
suppressMessages(pkgload::load_all(".", quiet = TRUE, export_all = FALSE, helpers = FALSE))
grDevices::pdf(NULL)
OUT <- file.path("dev", "audit", "stage2c")
CELLS <- c("U", "D", "E", "F")
CORES <- max(1L, parallel::detectCores() - 2L)
expit <- function(x) 1 / (1 + exp(-x))

gen <- function(k, r) {
  set.seed(300000 * k + r)
  cell <- CELLS[k]; N <- 800
  X <- stats::rnorm(N)
  e <- switch(cell, U = rep(0.5, N), D = , E = expit(-0.5 + X), F = expit(-0.5 + 2 * X))
  A <- stats::rbinom(N, 1, e)
  Tt <- stats::rexp(N, 0.08 * exp(0.7 * X)); C <- stats::runif(N, 0, 40)
  if (cell == "E") e <- stats::fitted(stats::glm(A ~ X, family = stats::binomial))
  data.frame(time = pmin(Tt, C), status = as.numeric(Tt <= C), treat = A, w = ifelse(A == 1, 1 / e, 1 / (1 - e)))
}

pkg_calls <- function(df, wn) {
  nw <- 0L
  res <- withCallingHandlers({
    dc <- df_counting(df, "time", "status", "treat", weight.name = wn)
    cr <- cox_rhogamma(dc, scheme = "fh", scheme_params = list(rho = 0, gamma = 0), draws = 0)
    c(z = dc$z.score, bhat = if (is.null(cr$fit)) NA_real_ else cr$fit$bhat,
      se_bhat = if (is.null(cr$fit)) NA_real_ else cr$fit$sig_bhat_asy)
  }, warning = function(w) { nw <<- nw + 1L; invokeRestart("muffleWarning") })
  c(res, n_warn = nw)
}

one_rep <- function(k, r) {
  df <- gen(k, r)
  out <- list()
  if (CELLS[k] == "U") {
    out$w <- pkg_calls(df, NULL)
    cx <- survival::coxph(survival::Surv(time, status) ~ treat, data = df, robust = TRUE)
  } else {
    out$w <- pkg_calls(df, "w")
    d10 <- df; d10$w <- 10 * df$w; out$w_x10 <- pkg_calls(d10, "w")
    dn <- df; dn$w <- df$w / sum(df$w); out$w_norm <- pkg_calls(dn, "w")
    cx <- survival::coxph(survival::Surv(time, status) ~ treat, data = df, weights = w, robust = TRUE)
  }
  b <- unname(stats::coef(cx)); se_r <- sqrt(cx$var[1, 1]); se_m <- sqrt(cx$naive.var[1, 1])
  rows <- lapply(names(out), function(s) data.frame(cell = CELLS[k], scaling = s, r = r, t(out[[s]])))
  cbind(do.call(rbind, rows), coxph_b = b, coxph_se_robust = se_r, coxph_se_model = se_m)
}

cat(sprintf("M4: HEAD %s; %d cores; %d replicates per cell; started %s\n", system("git rev-parse --short HEAD", intern = TRUE), CORES, R_REPS, format(Sys.time())))
t0 <- proc.time()[["elapsed"]]; reps <- list()
for (k in seq_along(CELLS)) {
  tk <- proc.time()[["elapsed"]]
  res <- parallel::mclapply(seq_len(R_REPS), function(r) one_rep(k, r), mc.cores = CORES, mc.preschedule = TRUE)
  bad <- vapply(res, inherits, logical(1), "try-error")
  if (any(bad)) stop("cell ", CELLS[k], ": ", as.character(res[[which(bad)[1]]]))
  reps[[k]] <- do.call(rbind, res)
  el <- proc.time()[["elapsed"]] - t0
  cat(sprintf("cell %s: %d replicates in %.0f s (total %.1f min)\n", CELLS[k], R_REPS, proc.time()[["elapsed"]] - tk, el / 60))
  if (el > 60 * 60 && k < length(CELLS)) { cat("STOP: wall-clock passed 60 min\n"); break }
}
rep <- do.call(rbind, reps)
saveRDS(rep, file.path(OUT, "km_m4_replicates.rds"), compress = "xz")

fin_stats <- function(x) { f <- is.finite(x); c(mean = mean(x[f]), sd = stats::sd(x[f]), rej = mean(abs(x[f]) > 1.96), n_nonfinite = sum(!f)) }
out <- list()
for (k in seq_along(reps)) for (s in unique(reps[[k]]$scaling)) {
  d <- reps[[k]][reps[[k]]$scaling == s, ]
  z <- fin_stats(d$z)
  fb <- is.finite(d$bhat) & is.finite(d$se_bhat)
  zr <- d$coxph_b / d$coxph_se_robust; zm <- d$coxph_b / d$coxph_se_model
  out[[length(out) + 1]] <- data.frame(cell = CELLS[k], scaling = s, n_reps = nrow(d),
    z_mean = z[["mean"]], z_sd = z[["sd"]], z_reject = z[["rej"]], z_nonfinite = z[["n_nonfinite"]],
    cox_bhat_sd = stats::sd(d$bhat[fb]), cox_mean_se = mean(d$se_bhat[fb]), cox_se_over_sd = mean(d$se_bhat[fb]) / stats::sd(d$bhat[fb]),
    cox_nonfinite = sum(!fb), cox_reject = mean(abs(d$bhat[fb] / d$se_bhat[fb]) > 1.96),
    coxph_robust_z_mean = mean(zr), coxph_robust_z_sd = stats::sd(zr), coxph_robust_reject = mean(abs(zr) > 1.96),
    coxph_b_sd = stats::sd(d$coxph_b), coxph_mean_se_robust = mean(d$coxph_se_robust),
    coxph_model_z_mean = mean(zm), coxph_model_z_sd = stats::sd(zm), coxph_model_reject = mean(abs(zm) > 1.96),
    coxph_mean_se_model = mean(d$coxph_se_model), mean_n_warn = mean(d$n_warn), stringsAsFactors = FALSE)
}
res <- do.call(rbind, out)
utils::write.csv(res, file.path(OUT, "km_m4_results.csv"), row.names = FALSE)
options(width = 260); print(res, row.names = FALSE, digits = 4)
cat(sprintf("\nwall-clock: %.1f min\n", (proc.time()[["elapsed"]] - t0) / 60))
