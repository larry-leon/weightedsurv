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
