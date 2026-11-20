# m6_common.R -- copied from dev/audit/oc/oc_common.R for M6 (TASK_km_fixes_stage2d.md, section 6):
# the data-generating model, the cells, make_data() and cox_robust(), verbatim; nothing under oc/ is edited.
# Added below the copy: the M6 cells 21-23.

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

# ---- section 6: proposed variance for the subject-weighted log-rank / Cox (as given) ----
cox_robust <- function(fit, beta = fit$bhat) {         # fit = cox_rhogamma(...)$fit ; psi is returned in the order of fit$time
  o <- order(fit$time); D <- fit$delta[o]; Z <- fit$z[o]; a <- fit$w_hat[o]; K <- fit$wt_rg[o]
  eb <- exp(beta * Z); S0 <- rev(cumsum(rev(a * eb))); S1 <- rev(cumsum(rev(a * eb * Z))); zb <- S1 / S0
  dL <- a * D / S0; cK <- cumsum(K * dL); cKz <- cumsum(K * zb * dL)
  psi <- numeric(length(o)); psi[o] <- a * (D * K * (Z - zb) - eb * (Z * cK - cKz))
  U <- sum(a * D * K * (Z - zb)); info <- sum(a * D * K * zb * (1 - zb))
  list(U = U, info = info, se = sqrt(sum(psi^2)) / info, z = U / sqrt(sum(psi^2)), psi = psi)
}

# ---- M6 cells (added) ----
# 21 Dt-400-null: gen(400, "null", "conf"), true weights, ceiling(time); seed 1e6 * 21 + r
# 22 D-400-null-strata: the data of cell 14 with s = ifelse(X > 0, "hi", "lo")
# 23 Dt-400-null-strata: the data of cell 21 with the same strata
m6_cells <- data.frame(
  k = c(1, 2, 5, 8, 10, 14, 16, 12, 21, 22, 23),
  cell = c("U-400", "R-400", "D-400", "E-400", "F-400", "D-400-null", "F-400-null", "Dt-400", "Dt-400-null",
           "D-400-null-strata", "Dt-400-null-strata"),
  truth = c(-0.4, -0.4, -0.4, -0.4, -0.4, 0, 0, -0.4, 0, 0, 0),
  null = c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE),
  strata = c(rep(FALSE, 9), TRUE, TRUE), stringsAsFactors = FALSE)
m6_data <- function(k, r) {
  if (k <= 19) return(make_data(k, r))
  if (k == 21) {
    set.seed(1e6 * 21 + r)
    df <- gen(400, "null", "conf")
    df$time <- ceiling(df$time)
    df$w <- ifelse(df$treat == 1, 1 / df$e, 1 / (1 - df$e))
    return(list(df = df, wname = "w"))
  }
  d <- if (k == 22) make_data(14, r) else m6_data(21, r)
  d$df$s <- ifelse(d$df$X > 0, "hi", "lo")
  d
}
