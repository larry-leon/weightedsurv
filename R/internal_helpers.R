# Internal helpers (not exported).

# Evaluate `expr` after set.seed(seed) and then restore the caller's RNG state:
# .Random.seed is put back exactly as it was, or removed again if it did not exist.
# `expr` is evaluated lazily in the caller's frame, after the seed is set.
#' @noRd
ws_with_seed <- function(seed, expr) {
  genv <- globalenv()
  had <- exists(".Random.seed", envir = genv, inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = genv, inherits = FALSE) else NULL
  on.exit({
    if (had) {
      assign(".Random.seed", old, envir = genv)
    } else if (exists(".Random.seed", envir = genv, inherits = FALSE)) {
      rm(".Random.seed", envir = genv)
    }
  }, add = TRUE)
  set.seed(seed)
  expr
}

# Are the weights "not all equal"? Decided on the positive weights only (Stage 2e, Rd): zero-weight
# rows contribute nothing, so they do not decide which variance formula runs. A 0/1 weight column
# therefore counts as equal weights, and gives the results of the positive-weight rows without weights.
#' @noRd
ws_weights_unequal <- function(w) length(unique(w[w > 0])) > 1

# Shared input validation for KM_diff(), df_counting() and plotKM.band_subgroups().
# Errors name the column and the problem; rows are never dropped. Returns the event
# indicator recoded to numeric 0/1 (accepted codings: 0/1, logical, and 1/2 read as
# survival::Surv reads it, 2 = event), plus the time, treatment, weight and strata vectors.
# `context` prefixes every message (e.g. "Subgroup 'er <= 0'").
#' @noRd
ws_validate_survdata <- function(df, tte.name, event.name, treat.name, weight.name = NULL,
                                 strata.name = NULL, context = NULL, require_arms = TRUE) {
  pre <- if (is.null(context)) "" else paste0(context, ": ")
  fail <- function(...) stop(pre, ..., call. = FALSE)
  cols <- c(tte.name, event.name, treat.name, weight.name, strata.name)
  missing_cols <- setdiff(cols, names(df))
  if (length(missing_cols) > 0) {
    fail("Missing required column(s) in df: ", paste(missing_cols, collapse = ", "))
  }
  for (nm in unique(cols)) {
    k <- sum(is.na(df[[nm]]))
    if (k > 0) fail(sprintf("Column '%s' has %d missing (NA) value(s); rows are not dropped, so remove or impute them first.", nm, k))
  }
  time <- df[[tte.name]]
  if (!is.numeric(time)) fail(sprintf("Column '%s' (time to event) must be numeric.", tte.name))
  if (any(!is.finite(time))) fail(sprintf("Column '%s' (time to event) has non-finite values.", tte.name))
  if (any(time < 0)) fail(sprintf("Column '%s' (time to event) has negative values; times must be >= 0.", tte.name))

  ev <- df[[event.name]]
  if (is.logical(ev)) {
    event <- as.numeric(ev)
  } else if (is.numeric(ev) && all(ev %in% c(0, 1))) {
    event <- as.numeric(ev)
  } else if (is.numeric(ev) && all(ev %in% c(1, 2)) && max(ev) == 2) {
    event <- as.numeric(ev) - 1
  } else {
    bad <- utils::head(sort(unique(as.character(ev[!(ev %in% c(0, 1))]))), 5)
    fail(sprintf("Column '%s' (event) must be coded 0/1 (1 = event), logical, or 1/2 (2 = event, as survival::Surv reads it); found: %s.",
                 event.name, paste(bad, collapse = ", ")))
  }

  treat <- df[[treat.name]]
  if (!all(treat %in% c(0, 1))) {
    fail(sprintf("Column '%s' (treatment) must be coded 0 = control, 1 = experimental.", treat.name))
  }
  weight <- NULL
  if (!is.null(weight.name)) {
    weight <- df[[weight.name]]
    if (!is.numeric(weight)) fail(sprintf("Column '%s' (weights) must be numeric.", weight.name))
    if (any(!is.finite(weight))) fail(sprintf("Column '%s' (weights) has non-finite values.", weight.name))
    if (any(weight < 0)) fail(sprintf("Column '%s' (weights) has negative values; weights must be >= 0.", weight.name))
  }
  if (require_arms) {
    pos <- if (is.null(weight)) rep(TRUE, length(time)) else weight > 0
    for (a in c(0, 1)) {
      if (!any(treat == a & pos)) {
        fail(sprintf("No observations%s in the %s arm (%s == %d).",
                     if (is.null(weight)) "" else " with positive weight",
                     if (a == 0) "control" else "treatment", treat.name, a))
      }
    }
  }
  list(time = time, event = event, treat = treat, weight = weight,
       strata = if (is.null(strata.name)) NULL else df[[strata.name]])
}

# Kaplan-Meier for one arm at arbitrary time points. The calculations run on the grid
# G = sort(unique(c(U, at_points))), which holds every event time of the arm, and the rows at
# match(at_points, G) are returned (duplicates in at_points are returned as given). On a grid
# that already holds every event time the values are those of the direct calculation.
#' @noRd
km_arm_on_grid <- function(U, D, W, at_points, draws = 0, seedstart = 816951) {
  G <- sort(unique(c(U, at_points)))
  re <- calculate_risk_event_counts(U, D, W, at_points = G, draws = draws, seedstart = seedstart)
  km <- KM_estimates(ybar = re$ybar, nbar = re$nbar, sig2w_multiplier = re$sig2w_multiplier)
  idx <- match(at_points, G)
  list(grid = G, idx = idx, S_grid = km$S_KM, surv = km$S_KM[idx], sig2 = km$sig2_KM[idx])
}

# Quantile time of a survival step curve, following the rule of survival::quantile.survfit:
# work on F = 1 - S, starting from F = 0 at time 0, with level p = 1 - qprob; keep the first
# time of each distinct value of F; the quantile is the midpoint between the first time F
# reaches p and the first time F exceeds p (both within tol), so that a plateau at the level
# gives the midpoint of the plateau; a plateau that runs to the end gives the midpoint between
# its start and the last time point. NA when the curve never reaches the level. NA values in
# survival_probs (e.g. limits where S = 0) are skipped.
#' @noRd
ws_km_quantile_time <- function(time_points, survival_probs, qprob, tol = sqrt(.Machine$double.eps)) {
  x <- c(min(0, time_points[1]), time_points)
  f <- c(0, 1 - survival_probs)
  p <- 1 - qprob
  if (all(is.na(f)) || max(f, na.rm = TRUE) < p) return(NA)
  x_last <- x[length(x)]
  first <- !duplicated(f)
  x <- x[first]
  f <- f[first]
  k <- length(f)
  i_reach <- stats::approx(f + tol, seq_len(k), p, method = "constant", f = 1)$y
  i_pass <- stats::approx(f - tol, seq_len(k), p, method = "constant", f = 1)$y
  q <- (x[i_reach] + x[i_pass]) / 2
  if (p == 0) q <- x[1]
  if (!is.na(f[k]) && abs(p - f[k]) < tol) q <- (x[i_reach] + x_last) / 2
  if (is.na(q)) NA else q
}

# Vertices of a Kaplan-Meier curve as plot_km_curves_counting() draws it (B10).
# The curve is the right-continuous step function from time.zero (S = 1) through the group's
# event times to its last observed time: at_points[idv] are the event times followed by the
# last observed time (df_counting()'s idv0 / idv1). x and y are drawn with type = "s".
# When se is given, lower / upper are the log-transform pointwise limits at the same vertices:
# 1 before the first event, upper capped at 1, and equal to the curve (0) where S = 0; they are
# step functions too. band_x / band_y trace the staircase polygon of the band.
#' @noRd
km_step_vertices <- function(at_points, S, idv, se = NULL, conf_level = 0.95, time.zero = 0) {
  x <- c(time.zero, at_points[idv])
  y <- c(1, S[idv])
  out <- list(x = x, y = y, lower = NULL, upper = NULL, band_x = NULL, band_y = NULL)
  if (!is.null(se)) {
    z <- stats::qnorm(1 - (1 - conf_level) / 2)
    s <- S[idv]; e <- se[idv]
    lo <- up <- s
    ok <- !is.na(s) & s > 0
    lo[ok] <- exp(log(s[ok]) - z * e[ok] / s[ok])
    up[ok] <- pmin(exp(log(s[ok]) + z * e[ok] / s[ok]), 1)
    out$lower <- c(1, lo)
    out$upper <- c(1, up)
    stair <- function(xx, v) {
      k <- length(xx)
      list(x = c(xx[1], rep(xx[-1], each = 2)), y = c(rep(v[-k], each = 2), v[k]))
    }
    su <- stair(x, out$upper)
    sl <- stair(x, out$lower)
    out$band_x <- c(su$x, rev(sl$x))
    out$band_y <- c(su$y, rev(sl$y))
  }
  out
}

# Cumulative integral of a curve d over a time grid with increments dt = diff(at_points), as
# used by cumulative_rmst_bands() for the estimate and for every resampled curve (B8).
# "trapezoid": R(t_j) = sum_{k <= j} (d(t_{k-1}) + d(t_k)) / 2 * (t_k - t_{k-1});
# "rectangle": R(t_j) = sum_{k <= j} d(t_{k-1}) * (t_k - t_{k-1}) (left rectangles: the exact area
# under the right-continuous step function d). R = 0 at the first grid point.
#' @noRd
ws_cumulative_mean <- function(d, dt, rule = c("trapezoid", "rectangle")) {
  rule <- match.arg(rule)
  h <- if (rule == "trapezoid") (utils::head(d, -1) + utils::tail(d, -1)) / 2 else utils::head(d, -1)
  c(0, cumsum(h * dt))
}

# Per-subject influence of a weighted Kaplan-Meier curve on a sorted grid (Stage 2c, R1 and R2).
# For subject i with time U_i, event indicator D_i and weight w_i,
#   psi_i(t) = D_i 1(U_i <= t) / (Y_w(U_i) - dN_w(U_i))
#              - sum over event times s <= min(t, U_i) of dN_w(s) / (Y_w(s) (Y_w(s) - dN_w(s))),
# a term being 0 where Y_w - dN_w <= 1e-10 Y_w (or Y_w = 0). Then S(t)^2 sum_i w_i^2 psi_i(t)^2 is
# the robust (infinitesimal-jackknife) variance of survfit(weights = w, robust = TRUE), and
# -S(t) sum_i G_i w_i psi_i(t) is one resampled deviation for multipliers G_i. Valid when the grid
# holds every event time of the group.
# On the grid g_1 <= ... <= g_m, psi_i(g_j) = D_i 1(k_i <= j) a[k_i] - H[min(j, k_i)], with k_i the
# first grid point >= U_i for an event and the last grid point <= U_i for a censored time (0 before
# the grid, m + 1 after it). psi_i is constant from j = k_i on, at `exit` = D_i a[k_i] - H[k_i].
#' @noRd
ws_km_psi_parts <- function(U, D, at_points, ybar, nbar) {
  m <- length(at_points)
  dN <- diff(c(0, nbar))
  ok <- ybar > 0 & (ybar - dN) > 1e-10 * ybar
  a <- ifelse(ok, 1 / (ybar - dN), 0)
  H <- cumsum(ifelse(ok, dN / (ybar * (ybar - dN)), 0))
  ev <- D == 1
  k <- ifelse(ev, findInterval(U, at_points, left.open = TRUE) + 1L, findInterval(U, at_points))
  in_grid <- k >= 1 & k <= m
  kk <- pmin(pmax(k, 1L), max(m, 1L))
  exit <- ifelse(in_grid, ifelse(ev, a[kk], 0) - H[kk], 0)
  list(m = m, k = k, in_grid = in_grid, H = H, exit = exit)
}

# Sum the rows of x (vector or matrix with one row per subject) by grid index idx in 1..m;
# returns an m-row matrix (zero rows where no subject has that index).
#' @noRd
ws_rowsum_grid <- function(x, idx, m) {
  x <- as.matrix(x)
  out <- matrix(0, m, ncol(x))
  if (length(idx)) {
    s <- rowsum(x, idx, reorder = TRUE)
    out[as.integer(rownames(s)), ] <- s
  }
  out
}

# Column-wise cumulative sums, keeping the matrix shape.
#' @noRd
ws_colcumsum <- function(x) matrix(apply(x, 2, cumsum), nrow = nrow(x))

# V(g_j) = sum_i w_i^2 psi_i(g_j)^2 for every grid point (R1), in O(n + m):
# subjects with k_i <= j contribute their exit value, the others H_j^2 w_i^2.
#' @noRd
ws_km_ij_var <- function(parts, W) {
  m <- parts$m
  if (m == 0) return(numeric(0))
  w2 <- W^2
  ig <- parts$in_grid
  C <- cumsum(ws_rowsum_grid(w2[ig] * parts$exit[ig]^2, parts$k[ig], m)[, 1])
  Q <- sum(w2[parts$k >= 1]) - cumsum(ws_rowsum_grid(w2[ig], parts$k[ig], m)[, 1])
  C + parts$H^2 * Q
}

# M[j, b] = sum_i G[i, b] w_i psi_i(g_j) for every grid point and draw (R2), in O((n + m) B).
#' @noRd
ws_km_ij_draws <- function(parts, W, G) {
  m <- parts$m
  G <- as.matrix(G)
  if (m == 0) return(matrix(0, 0, ncol(G)))
  ig <- parts$in_grid
  Gw <- G * W
  E <- ws_colcumsum(ws_rowsum_grid(Gw[ig, , drop = FALSE] * parts$exit[ig], parts$k[ig], m))
  X <- ws_colcumsum(ws_rowsum_grid(Gw[ig, , drop = FALSE], parts$k[ig], m))
  tot <- colSums(Gw[parts$k >= 1, , drop = FALSE])
  still <- -sweep(X, 2, tot, "-")  # sum over subjects with k_i > j
  E - parts$H * still
}

# The robust variance or resampled deviations for one group on any at_points: computed on the
# sorted grid and returned in the order of at_points.
#' @noRd
ws_km_ij_sorted <- function(U, D, W, at_points, ybar, nbar, G = NULL) {
  o <- order(at_points)
  srt <- !is.unsorted(at_points)
  g <- if (srt) at_points else at_points[o]
  yb <- if (srt) ybar else ybar[o]
  nb <- if (srt) nbar else nbar[o]
  parts <- ws_km_psi_parts(U, D, g, yb, nb)
  if (is.null(G)) {
    v <- ws_km_ij_var(parts, W)
    if (srt) v else v[order(o)]
  } else {
    M <- ws_km_ij_draws(parts, W, G)
    if (srt) M else M[order(o), , drop = FALSE]
  }
}

# Per-subject score contributions of the subject-weighted log-rank / Cox score (Stage 2d, L1-L4).
# Subject i: time U_i, event D_i, arm Z_i (0/1), subject weight a_i; K(s) the time weight of the
# scheme; sums over s run over distinct times.
#   S0(s, b)   = sum_j a_j 1(U_j >= s) exp(b Z_j)
#   zbar(s, b) = sum_j a_j 1(U_j >= s) Z_j exp(b Z_j) / S0(s, b)
#   dL(s, b)   = sum_j a_j D_j 1(U_j = s) / S0(s, b)
#   psi_i(b)   = a_i [ D_i K(U_i) (Z_i - zbar(U_i, b)) - exp(b Z_i) sum_{s <= U_i} K(s) (Z_i - zbar(s, b)) dL(s, b) ]
# The package's score is -sum_i psi_i(b); sum_i psi_i(b)^2 estimates its variance (for FH(0,0) the
# robust variance of coxph(weights = , robust = TRUE, ties = "breslow")), and -sum_i G_i psi_i(b) is
# one resampled score. `K` holds one time weight per observation; the weight of a tied time is the
# one at its first observation (the one the score uses). A term whose S0 is 0 is 0.
# Returns psi in the order of the input, and, for the cumulative version, the distinct times with
# cK(t) = sum_{s <= t} K(s) dL(s) and cKz(t) = sum_{s <= t} K(s) zbar(s) dL(s).
#' @noRd
ws_score_psi <- function(time, delta, z, a, K = rep(1, length(time)), beta = 0) {
  n <- length(time)
  if (n == 0) return(list(psi = numeric(0), times = numeric(0), cK = numeric(0), cKz = numeric(0)))
  o <- order(time)
  tt <- time[o]; D <- delta[o]; Z <- z[o]; w <- a[o]
  first <- !duplicated(tt)
  first_idx <- which(first)
  grp <- cumsum(first)
  eb <- exp(beta * Z)
  S0 <- rev(cumsum(rev(w * eb)))[first_idx]
  S1 <- rev(cumsum(rev(w * eb * Z)))[first_idx]
  ok <- S0 > 0
  zb <- ifelse(ok, S1 / S0, 0)
  dNw <- as.vector(rowsum(w * D, grp))
  dL <- ifelse(ok, dNw / S0, 0)
  Kd <- K[o][first_idx]
  Kd[dL == 0] <- 0                      # a time without (positive-weight) events contributes nothing
  cK <- cumsum(Kd * dL)
  cKz <- cumsum(Kd * zb * dL)
  psi <- numeric(n)
  psi[o] <- w * (D * Kd[grp] * (Z - zb[grp]) - eb * (Z * cK[grp] - cKz[grp]))
  list(psi = psi, times = tt[first_idx], cK = cK, cKz = cKz)
}

# Subject-level data for the weighted log-rank variance (L1): from the attribute set by
# df_counting(), or from a df_counting() result passed directly (elements time, delta, z, w_hat).
# NULL when there are none, when they do not fit the grid, or when the weights are all equal --
# the earlier formula then runs unchanged.
#' @noRd
ws_lr_subject_data <- function(x, at_points) {
  s <- attr(x, "ws_subject", exact = TRUE)
  if (is.null(s)) {
    nms <- c("time", "delta", "z", "w_hat")
    if (!is.list(x) || !all(nms %in% names(x))) return(NULL)
    s <- list(time = x[["time"]], delta = x[["delta"]], z = x[["z"]], w_hat = x[["w_hat"]])
    s$unequal <- ws_weights_unequal(s$w_hat)
    # a df_counting() result made with weight_variance = "estimated": the pieces of the propensity model
    # (rows in the order of `time`); the influence-based formulas then run whatever the weights
    if (identical(x[["weight_variance"]], "estimated") && is.list(x[["ew_parts"]])) {
      s$ew <- x[["ew_parts"]]
      s$unequal <- TRUE
    }
  }
  if (!is.null(s$ew) && !(is.matrix(s$ew$IF) && is.matrix(s$ew$dw) && nrow(s$ew$IF) == length(s$time) &&
                          nrow(s$ew$dw) == length(s$time))) s$ew <- NULL
  ok <- is.numeric(s$time) && is.numeric(s$delta) && is.numeric(s$z) && is.numeric(s$w_hat) &&
    length(s$delta) == length(s$time) && length(s$z) == length(s$time) && length(s$w_hat) == length(s$time)
  if (!ok || !isTRUE(s$unequal)) return(NULL)
  if (is.unsorted(at_points) || anyNA(match(s$time, at_points))) return(NULL)
  s
}

# V(t) = sum_i (contribution of subject i to the score up to t)^2 on a sorted grid (L1, cumulative
# version, b = 0). The contribution is psi_i when U_i <= t and -a_i (Z_i cK(t) - cKz(t)) when U_i > t, so
#   V(t) = sum_{U_i <= t} psi_i^2 + (cK(t) - cKz(t))^2 sum_{U_i > t, Z_i = 1} a_i^2 + cKz(t)^2 sum_{U_i > t, Z_i = 0} a_i^2.
#' @noRd
ws_lr_cumvar <- function(sp, time, z, a, at_points) {
  m <- length(at_points)
  k <- match(time, at_points)
  A <- cumsum(ws_rowsum_grid(sp$psi^2, k, m)[, 1])
  a2 <- a^2
  R1 <- sum(a2[z == 1]) - cumsum(ws_rowsum_grid(a2 * (z == 1), k, m)[, 1])
  R0 <- sum(a2[z == 0]) - cumsum(ws_rowsum_grid(a2 * (z == 0), k, m)[, 1])
  j <- findInterval(at_points, sp$times)
  cK <- c(0, sp$cK)[j + 1]
  cKz <- c(0, sp$cKz)[j + 1]
  A + (cK - cKz)^2 * R1 + cKz^2 * R0
}
