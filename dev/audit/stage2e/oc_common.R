# oc_common.R -- shared by oc_sim.R, oc_boot.R, oc_summarise.R (TASK_km_oc_evaluation.md)
# Data-generating model, cells, seeds, truth, and the script-side reference functions of
# sections 6 and 7 (cox_robust, ew_parts, ew_km), used as given in the task document.

# ---- section 3: data-generating model ----
Lam    <- function(t, a, sc) if (a == 0 || sc == "null") 0.04 * t else
                             if (sc == "ph") 0.04 * exp(-0.4) * t else 0.04 * pmin(t, 6) + 0.02 * pmax(t - 6, 0)
LamInv <- function(h, a, sc) if (a == 0 || sc == "null") h / 0.04 else
                             if (sc == "ph") h / (0.04 * exp(-0.4)) else ifelse(h <= 0.24, h / 0.04, 6 + (h - 0.24) / 0.02)

gen <- function(N, sc, assign) {                       # assign = "rand", "conf" or "strong"
  X <- rnorm(N)
  H <- -log(pnorm(0.6 * X + 0.8 * rnorm(N)))           # standard exponential, larger X -> earlier event
  e <- switch(assign, rand = rep(0.5, N), conf = plogis(-0.5 + X), strong = plogis(-0.5 + 2 * X))
  A <- rbinom(N, 1, e)
  T <- ifelse(A == 1, LamInv(H, 1, sc), LamInv(H, 0, sc))
  C <- runif(N, 0, 72)
  data.frame(X = X, e = e, treat = A, time = pmin(T, C), status = as.numeric(T <= C))
}

tt <- c(1, 3, 6, 9, 12, 18, 24, 30, 36, 42, 48)
band_range <- c(3, 42)

# ---- cells ----
cells <- data.frame(
  cell_id = 1:19,
  cell = c("U-400", "R-400", "S-400", "D-200", "D-400", "D-800", "E-200", "E-400", "E-800", "F-400", "F-800",
           "Dt-400", "U-400-null", "D-400-null", "E-400-null", "F-400-null", "U-400-delay", "D-400-delay", "E-400-delay"),
  sc = c(rep("ph", 12), rep("null", 4), rep("delay", 3)),
  wt = c("U", "R", "S", "D", "D", "D", "E", "E", "E", "F", "F", "Dt", "U", "D", "E", "F", "U", "D", "E"),
  N = c(400, 400, 400, 200, 400, 800, 200, 400, 800, 400, 800, 400, 400, 400, 400, 400, 400, 400, 400),
  stringsAsFactors = FALSE)
cells$assign <- c(U = "rand", R = "rand", S = "rand", D = "conf", E = "conf", F = "strong", Dt = "conf")[cells$wt]
cells$loghr_truth <- ifelse(cells$sc == "ph" & cells$cell_id != 3, -0.4, ifelse(cells$sc == "null", 0, NA))
cells$do_rmst <- cells$cell_id != 12
cells$do_cox <- cells$cell_id != 12
cells$ew <- cells$wt == "E"
cells$d20 <- cells$cell_id %in% c(5, 8)

# Replicate r of cell k: data and weights, complete before any package call. Returns the data
# frame (column w always present; w = 1 in U cells), the weight column name to pass (NULL in U
# cells), the fitted propensity model (E cells) and the checksum.
make_data <- function(k, r) {
  cl <- cells[k, ]
  set.seed(1e6 * k + r)
  df <- gen(cl$N, cl$sc, cl$assign)
  g <- NULL
  df$w <- switch(cl$wt,
    U = rep(1, cl$N),
    R = runif(cl$N, 0.5, 2),
    S = exp(0.5 * df$X - 0.125),
    D = , F = ifelse(df$treat == 1, 1 / df$e, 1 / (1 - df$e)),
    E = { g <- glm(treat ~ X, binomial, data = df); eh <- fitted(g); ifelse(df$treat == 1, 1 / eh, 1 / (1 - eh)) },
    Dt = { df$time <- ceiling(df$time); ifelse(df$treat == 1, 1 / df$e, 1 / (1 - df$e)) })
  list(df = df, wname = if (cl$wt == "U") NULL else "w", g = g,
       chk = sprintf("%.10f", sum(df$time) + sum(df$status) + sum(df$w)))
}

# ---- truth ----
S_true <- function(t, a, k) {
  sc <- cells$sc[k]
  s <- exp(-Lam(t, a, sc))
  if (k == 3) s <- pnorm(qnorm(exp(-Lam(t, a, "ph"))) - 0.3)
  s
}
delta_true <- function(t, k) S_true(t, 1, k) - S_true(t, 0, k)
# R(t) = integral of delta over [0, t], stats::integrate with rel.tol = 1e-10; the range is split at
# month 6 (the kink of the "delay" scenario) so that each piece is smooth.
R_true1 <- function(t, k) {
  f <- function(u) delta_true(u, k)
  if (t <= 0) return(0)
  a <- stats::integrate(f, 0, min(t, 6), rel.tol = 1e-10, abs.tol = 0)$value
  if (t > 6) a <- a + stats::integrate(f, 6, t, rel.tol = 1e-10, abs.tol = 0)$value
  a
}
R_true <- function(t, k) {
  if (cells$sc[k] == "null") return(rep(0, length(t)))
  vapply(t, R_true1, numeric(1), k = k)
}

# ---- section 6: proposed variance for the subject-weighted log-rank / Cox (as given) ----
cox_robust <- function(fit, beta = fit$bhat) {         # fit = cox_rhogamma(...)$fit ; psi is returned in the order of fit$time
  o <- order(fit$time); D <- fit$delta[o]; Z <- fit$z[o]; a <- fit$w_hat[o]; K <- fit$wt_rg[o]
  eb <- exp(beta * Z); S0 <- rev(cumsum(rev(a * eb))); S1 <- rev(cumsum(rev(a * eb * Z))); zb <- S1 / S0
  dL <- a * D / S0; cK <- cumsum(K * dL); cKz <- cumsum(K * zb * dL)
  psi <- numeric(length(o)); psi[o] <- a * (D * K * (Z - zb) - eb * (Z * cK - cKz))
  U <- sum(a * D * K * (Z - zb)); info <- sum(a * D * K * zb * (1 - zb))
  list(U = U, info = info, se = sqrt(sum(psi^2)) / info, z = U / sqrt(sum(psi^2)), psi = psi)
}

# ---- section 7: estimated weights carried through (as given) ----
ew_parts <- function(g, treat) {                       # g = glm(treat ~ ..., binomial) fitted on the analysis rows
  xm <- model.matrix(g); e <- fitted(g)
  list(IF = (xm * (treat - e)) %*% solve(crossprod(xm * sqrt(e * (1 - e)))),   # n x p, influence of the coefficients
       dw = ifelse(treat == 1, -(1 - e) / e, e / (1 - e)) * xm)                # n x p, derivative of the weights
}

ew_km <- function(P0, P1, S0, S1, idx0, idx1, w, IF, dw) {     # S0, S1 = surv0, surv1 on the grid
  n <- nrow(IF)
  Dg0 <- crossprod(P0, dw[idx0, , drop = FALSE]); Dg1 <- crossprod(P1, dw[idx1, , drop = FALSE])   # m x p
  Phi0 <- IF %*% t(Dg0); Phi0[idx0, ] <- Phi0[idx0, ] + w[idx0] * P0                                 # n x m
  Phi1 <- IF %*% t(Dg1); Phi1[idx1, ] <- Phi1[idx1, ] + w[idx1] * P1
  list(Dg0 = Dg0, Dg1 = Dg1, sig2_surv0 = S0^2 * colSums(Phi0^2), sig2_surv1 = S1^2 * colSums(Phi1^2),
       sig2_dhat = colSums((Phi1 * rep(S1, each = n) - Phi0 * rep(S0, each = n))^2))
}

# ---- helpers ----
# index of the last grid point <= t (NA when there is none)
last_le <- function(grid, t) { j <- findInterval(t, grid); ifelse(j >= 1, j, NA_integer_) }

parse_args <- function(defaults) {
  a <- commandArgs(trailingOnly = TRUE)
  for (x in a) {
    kv <- strsplit(sub("^--", "", x), "=", fixed = TRUE)[[1]]
    defaults[[kv[1]]] <- if (length(kv) > 1) paste(kv[-1], collapse = "=") else TRUE
  }
  defaults
}
oc_dir <- function() {
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) dirname(normalizePath(f[1])) else normalizePath("dev/audit/oc")
}

# =============================================================================
# Stage 2e, M7 (TASK_km_fixes_stage2e.md, section 6): additions to the copy of dev/audit/oc/oc_common.R.
# Everything above is unchanged. The M7 cells, the generator gen_v(), the data of a replicate and the truth
# by scenario follow.
# =============================================================================
gen_v <- function(N, sc) {            # shaped like the vignette examples: about 14% treated, treated-arm ESS about a third of the arm
  X <- rnorm(N); H <- -log(pnorm(0.6 * X + 0.8 * rnorm(N)))
  e <- plogis(-2.3 + 1.3 * X); A <- rbinom(N, 1, e)
  T <- ifelse(A == 1, LamInv(H, 1, sc), LamInv(H, 0, sc)); C <- runif(N, 0, 72)
  data.frame(X = X, e = e, treat = A, time = pmin(T, C), status = as.numeric(T <= C))
}

m7 <- data.frame(
  k    = c(7, 8, 9, 15, 19, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40),
  cell = c("E-200", "E-400", "E-800", "E-400-null", "E-400-delay", "Es-400", "Es-400-null", "Es-400-null-strata",
           "Et-400-null", "Fe-400", "Fe-800", "Fe-400-null", "V-1450-null", "V-1450-ph", "Sub-800-null"),
  sc   = c("ph", "ph", "ph", "null", "delay", "ph", "null", "null", "null", "ph", "ph", "null", "null", "ph", "null"),
  type = c(rep("ipw", 5), "stabilized", "stabilized", "stabilized", "ipw", "ipw", "ipw", "ipw", "stabilized", "stabilized", "ipw"),
  reps = c(rep(500, 5), rep(2000, 10)),
  stringsAsFactors = FALSE)
m7$repro  <- m7$k %in% c(7, 8, 9, 15, 19)          # reproduction of the evaluation's script values (requirement 1)
m7$strat  <- m7$k == 33                             # stratified log-rank only
m7$sub    <- m7$k == 40                             # subgroup analyses A, B, C
m7$band   <- !(m7$k %in% c(33, 40))
m7$rmst   <- !(m7$k %in% c(33, 34, 40))
m7$fh01   <- !(m7$k %in% c(38, 39))                 # FH(0,1) as well as FH(0,0)
m7$both   <- m7$k %in% c(31, 32)                    # also the "ipw" weights of the same model (requirement 2), and wlr_dhat_estimates()
m7$coxph  <- m7$k == 34                             # requirement 4
m7$loghr_truth <- ifelse(m7$sc == "ph", -0.4, ifelse(m7$sc == "null", 0, NA))

# Data of replicate r of cell k, complete before any package call. Cells 7, 8, 9, 15, 19 are the evaluation's
# cells with their seeds; cells 31-33 take the data of the evaluation's cells 8 and 15 with their seeds; the
# others are seeded 1e6 * k + r. `g` is the evaluation's own glm where there is one (ps_weights(model = g)).
m7_data <- function(k, r) {
  if (k %in% c(7, 8, 9, 15, 19)) { d <- make_data(k, r); return(list(df = d$df, g = d$g, chk = d$chk)) }
  if (k %in% c(31, 32, 33)) {
    d <- make_data(if (k == 31) 8 else 15, r); df <- d$df
    if (k == 33) df$s <- ifelse(df$X > 0, "hi", "lo")
    return(list(df = df, g = d$g, chk = d$chk))
  }
  set.seed(1e6 * k + r)
  df <- switch(as.character(k),
    "34" = { x <- gen(400, "null", "conf"); x$time <- ceiling(x$time); x },
    "35" = gen(400, "ph", "strong"), "36" = gen(800, "ph", "strong"), "37" = gen(400, "null", "strong"),
    "38" = gen_v(1450, "null"), "39" = gen_v(1450, "ph"), "40" = gen(800, "null", "conf"))
  list(df = df, g = NULL, chk = sprintf("%.10f", sum(df$time) + sum(df$status) + sum(df$X)))
}

# Truth by scenario (the evaluation's: S_a(t) = exp(-Lam(t, a, sc)); the same for gen_v(), whose H is standard
# exponential). R(t), the integral of the difference over [0, t]: 0 under "null"; under "ph" the closed form of
# the evaluation's integral (both curves are exponential); under "delay" the evaluation's stats::integrate().
delta_sc <- function(t, sc) exp(-Lam(t, 1, sc)) - exp(-Lam(t, 0, sc))
R_sc <- function(t, sc) {
  if (sc == "null") return(rep(0, length(t)))
  if (sc == "ph") { c1 <- 0.04 * exp(-0.4); c0 <- 0.04; return((1 - exp(-c1 * t)) / c1 - (1 - exp(-c0 * t)) / c0) }
  f <- function(u) delta_sc(u, sc)
  vapply(t, function(x) {
    if (x <= 0) return(0)
    a <- stats::integrate(f, 0, min(x, 6), rel.tol = 1e-10, abs.tol = 0)$value
    if (x > 6) a <- a + stats::integrate(f, 6, x, rel.tol = 1e-10, abs.tol = 0)$value
    a
  }, numeric(1))
}
