# Reference implementation of the per-subject score contributions of the subject-weighted
# log-rank / Cox score (TASK_km_fixes_stage2d.md, section 2), used as given.
# time, delta, z (0/1), a = subject weights, K = time weight per observation (aligned with `time`), beta
ws_score_resid_ref <- function(time, delta, z, a, K = rep(1, length(time)), beta = 0) {
  o <- order(time); tt <- time[o]; D <- delta[o]; Z <- z[o]; w <- a[o]; Kk <- K[o]; n <- length(tt)
  first_idx <- which(!duplicated(tt)); grp <- cumsum(!duplicated(tt))        # distinct times and tie-group id
  eb <- exp(beta * Z)
  S0 <- rev(cumsum(rev(w * eb)))[first_idx]                                   # risk-set sums at distinct times
  zb <- rev(cumsum(rev(w * eb * Z)))[first_idx] / S0                          # weighted proportion treated in the risk set
  dNw <- as.vector(rowsum(w * D, grp))                                        # weighted events at each distinct time
  Kd <- Kk[first_idx]                                                         # time weight at each distinct time
  dL <- dNw / S0; cK <- cumsum(Kd * dL); cKz <- cumsum(Kd * zb * dL)
  psi <- numeric(n)
  psi[o] <- w * (D * Kd[grp] * (Z - zb[grp]) - eb * (Z * cK[grp] - cKz[grp]))
  list(psi = psi, U = sum(psi), info = sum(Kd * zb * (1 - zb) * dNw))        # package score = -U ; info = i_bhat
}

# Confounded two-arm data with inverse-probability weights and two strata (seeded; RNG state restored).
ws_sim_ipw <- function(n = 300, seed = 4, ties = FALSE) {
  old <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) get(".Random.seed", envir = globalenv()) else NULL
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else assign(".Random.seed", old, envir = globalenv()))
  set.seed(seed)
  x <- stats::rnorm(n); e <- stats::plogis(-0.5 + x); a <- stats::rbinom(n, 1, e)
  tt <- stats::rexp(n, 0.05 * exp(0.6 * x - 0.4 * a)); cc <- stats::runif(n, 0, 60)
  data.frame(time = if (ties) ceiling(pmin(tt, cc)) else pmin(tt, cc), status = as.numeric(tt <= cc), treat = a,
             w = ifelse(a == 1, 1 / e, 1 / (1 - e)), s = ifelse(x > 0, "hi", "lo"), one = 1, stringsAsFactors = FALSE)
}
ws_quiet <- function(expr) { r <- NULL; invisible(utils::capture.output(r <- suppressWarnings(expr))); r }
ws_schemes <- function(d) list(
  fh00 = list("fh", list(rho = 0, gamma = 0)), fh01 = list("fh", list(rho = 0, gamma = 1)),
  MB = list("MB", list(mb_tstar = unname(stats::quantile(d$time, 0.5)))), schemper = list("schemper", list()))
# Evaluate expr, discarding printed output but letting warnings through (for expect_no_warning()).
ws_capture <- function(expr) { r <- NULL; invisible(utils::capture.output(r <- expr)); r }
