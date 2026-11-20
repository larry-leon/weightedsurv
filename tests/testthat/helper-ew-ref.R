# Reference implementation for estimated weights (TASK_km_fixes_stage2e.md, section 2), used as given.

# Propensity model: weights, influence of the coefficients (n x p), derivative of each weight (n x p)
ws_ps_ref <- function(formula, data, type = c("ipw", "stabilized")) {
  type <- match.arg(type)
  g  <- stats::glm(formula, family = stats::binomial(), data = data)
  xm <- stats::model.matrix(g); e <- as.numeric(stats::fitted(g)); A <- as.numeric(g$y)
  IF <- (xm * (A - e)) %*% solve(crossprod(xm * sqrt(e * (1 - e))))        # coefficient estimate minus its limit ~ colSums(IF)
  w  <- ifelse(A == 1, 1 / e, 1 / (1 - e))
  dw <- ifelse(A == 1, -(1 - e) / e, e / (1 - e)) * xm
  if (type == "stabilized") { p1 <- mean(A); sc <- ifelse(A == 1, p1, 1 - p1); w <- sc * w; dw <- sc * dw }
  list(weights = w, IF = IF, dw = dw, treat = A)
}

# Kaplan-Meier of the group `rows` (row numbers), on a sorted grid that holds the group's event times.
# time, delta, w: all n subjects; IF, dw: n x p in the same row order. Returns S, Phi (n x m) and the variance.
ws_km_ew_ref <- function(time, delta, w, IF, dw, rows, grid) {
  U <- time[rows]; D <- delta[rows]; W <- w[rows]
  ybar <- colSums(outer(U, grid, ">=") * W)
  dN   <- diff(c(0, colSums(outer(U[D == 1], grid, "<=") * W[D == 1])))
  ok   <- ybar > 0 & (ybar - dN) > 1e-10 * ybar
  psi  <- t(apply(sweep(outer(U, grid, "==") * D, 2, ifelse(ok, 1 / (ybar - dN), 0), "*") -
                  sweep(outer(U, grid, ">="), 2, ifelse(ok, dN / (ybar * (ybar - dN)), 0), "*"), 1, cumsum))   # n_g x m
  S    <- cumprod(1 - ifelse(ybar > 0, dN / ybar, 0))
  Dg   <- crossprod(psi, dw[rows, , drop = FALSE])                          # m x p
  Phi  <- IF %*% t(Dg); Phi[rows, ] <- Phi[rows, ] + W * psi                # n x m
  list(S = S, Dg = Dg, Phi = Phi, sig2 = S^2 * colSums(Phi^2))
}
# Difference of the arms, k1 and k0 from ws_km_ew_ref(), n = nrow(IF):
ws_km_ew_diff_ref <- function(k1, k0) { n <- nrow(k1$Phi); colSums((k1$Phi * rep(k1$S, each = n) - k0$Phi * rep(k0$S, each = n))^2) }

# Score contributions with estimated weights: psi from ws_score_resid_ref() (within stratum), a = weights, all n subjects
ws_score_ew_ref <- function(psi, a, IF, dw) as.numeric(psi + IF %*% colSums((psi / a) * dw))

# Confounded two-arm data with two covariates for a fitted propensity model (seeded; RNG state restored).
ws_sim_ps <- function(n = 400, seed = 11, ties = FALSE) {
  old <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) get(".Random.seed", envir = globalenv()) else NULL
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else assign(".Random.seed", old, envir = globalenv()))
  set.seed(seed)
  x1 <- stats::rnorm(n); x2 <- stats::rbinom(n, 1, 0.4); A <- stats::rbinom(n, 1, stats::plogis(-0.4 + 0.9 * x1 - 0.5 * x2))
  tt <- stats::rexp(n, 0.04 * exp(0.5 * x1 - 0.4 * A)); cc <- stats::runif(n, 0, 72)
  data.frame(time = if (ties) ceiling(pmin(tt, cc)) else pmin(tt, cc), status = as.numeric(tt <= cc), treat = A,
             x1 = x1, x2 = x2, s = ifelse(x2 == 1, "hi", "lo"), stringsAsFactors = FALSE)
}
ws_rel <- function(a, b) max(abs(a - b)) / max(abs(b))
