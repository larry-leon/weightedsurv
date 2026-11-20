# Internal helpers for weight_variance = "estimated" (Stage 2e, W2-W4). Not exported.
#
# The weights come from a logistic propensity model (ps_weights()): IF (n x p) is the influence of
# each subject on the fitted coefficients, dw (n x p) the derivative of each weight in the
# coefficients. For a Kaplan-Meier curve of a group g, with psi_i(t) the influence value of
# ws_km_psi_parts() computed from the group's own weighted counts,
#   Dg(t)    = sum_{k in g} psi_k(t) dw_k
#   Phi_i(t) = 1(i in g) w_i psi_i(t) + IF_i' Dg(t)          for every subject i = 1..n
#   var S_g(t) = S_g(t)^2 sum_i Phi_i(t)^2.
# For the score of the log-rank test and Cox estimate, with psi_i(b) of ws_score_psi(),
#   B(b) = sum_j (psi_j(b) / a_j) dw_j,  psi^ew_i(b) = psi_i(b) + IF_i' B(b),  var = sum_i psi^ew_i(b)^2.

# Validation of the two arguments, shared by df_counting(), KM_diff(), cumulative_rmst_bands() and
# plotKM.band_subgroups(). Returns NULL for "known" (weight_model is then ignored) and
# list(IF, dw, CIF = crossprod(IF)) for "estimated", rows in the order of df.
#' @noRd
ws_ew_setup <- function(weight_variance, weight_model, df, treat.name, weight.name, draws = 0, fn = "KM_diff") {
  if (weight_variance != "estimated") return(NULL)
  fail <- function(...) stop(fn, "(): weight_variance = \"estimated\": ", ..., call. = FALSE)
  if (is.null(weight.name)) fail("weight.name must name the column that holds the weights of weight_model.")
  if (is.null(weight_model)) fail("weight_model is missing; it must be the object returned by ps_weights().")
  if (!inherits(weight_model, "ps_weights")) {
    fail("weight_model must be the object returned by ps_weights(); it has class '", class(weight_model)[1], "'.")
  }
  n <- length(weight_model$weights)
  if (nrow(df) != n) {
    fail(sprintf("df has %d rows but weight_model was fitted on %d; the propensity model must be fitted on exactly the rows analysed.", nrow(df), n))
  }
  if (!is.matrix(weight_model$IF) || !is.matrix(weight_model$dw) || nrow(weight_model$IF) != n || nrow(weight_model$dw) != n ||
      ncol(weight_model$IF) != ncol(weight_model$dw)) {
    fail("weight_model is not a complete ps_weights object (IF and dw must be matrices with one row per subject).")
  }
  if (!isTRUE(all(df[[treat.name]] == weight_model$treat))) {
    fail(sprintf("the treatment column '%s' differs from the treatment of weight_model (same rows in the same order are needed).", treat.name))
  }
  if (!isTRUE(all.equal(as.numeric(df[[weight.name]]), as.numeric(weight_model$weights), tolerance = 1e-8))) {
    fail(sprintf("the weight column '%s' differs from the weights of weight_model (same rows in the same order, and the same type, are needed).", weight.name))
  }
  if (!is.numeric(draws) || length(draws) != 1 || draws != 0) {
    fail("draws must be 0 (the resampling-based pointwise variance is not available with estimated weights).")
  }
  list(IF = weight_model$IF, dw = weight_model$dw, CIF = crossprod(weight_model$IF))
}

# Reorder the rows of IF and dw (the data were sorted by time).
#' @noRd
ws_ew_reorder <- function(ew, ord) {
  if (is.null(ew)) return(NULL)
  ew$IF <- ew$IF[ord, , drop = FALSE]
  ew$dw <- ew$dw[ord, , drop = FALSE]
  ew
}

# One group (rows `rows` of the n subjects) on a sorted grid that holds the group's event times.
# ybar, nbar: the group's weighted counts on the grid. Returns, per grid point,
#   Dg (m x p), A = sum_{i in g} w_i psi_i IF_i (m x p), sumsq = sum_i Phi_i^2, and the parts.
# sum_i Phi_i^2 = sum_{i in g} w_i^2 psi_i^2 + 2 A . Dg + Dg' {sum_i IF_i IF_i'} Dg.
#' @noRd
ws_km_ew_group <- function(U, D, W, rows, grid, ybar, nbar, ew) {
  parts <- ws_km_psi_parts(U, D, grid, ybar, nbar)
  Dg <- ws_km_ij_draws(parts, rep(1, length(U)), ew$dw[rows, , drop = FALSE])
  A <- ws_km_ij_draws(parts, W, ew$IF[rows, , drop = FALSE])
  V <- ws_km_ij_var(parts, W)
  sumsq <- V + 2 * rowSums(A * Dg) + rowSums((Dg %*% ew$CIF) * Dg)
  list(parts = parts, Dg = Dg, A = A, sumsq = pmax(sumsq, 0))
}

# Variance of the Kaplan-Meier curve S of one group at at_points (sorted, holding its event times).
#' @noRd
ws_km_ew_sig2 <- function(U, D, W, rows, at_points, ybar, nbar, S, ew) {
  v <- S^2 * ws_km_ew_group(U, D, W, rows, at_points, ybar, nbar, ew)$sumsq
  v[S == 0] <- 0
  v
}
