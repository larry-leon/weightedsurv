# ---- Utility Functions ----
#' Safe execution wrapper
#'
#' Executes an R expression safely, returning NULL and printing an error message if an error occurs.
#'
#' @param expr An R expression to evaluate.
#' @return The result of expr, or NULL if an error occurs.
#' @export

safe_run <- function(expr) {
  tryCatch(expr, error = function(e) {
    message("Error: ", e$message)
    NULL
  })
}


#' Get df_counting results
#'
#' Wrapper for df_counting, safely executes and returns results for survival analysis.
#'
#' @param ... Arguments passed to \code{\link{df_counting}}.
#' @return Result from df_counting or NULL if error.
#' @seealso \code{\link{df_counting}}
#' @export
get_dfcounting <- function(...) {
  safe_run({
    df_counting(...)
  })
}

#' Check and Compare Statistical Test Results
#'
#' Calculates and compares three equivalent test statistics from weighted survival analysis:
#' the squared standardized weighted log-rank statistic, the log-rank chi-squared statistic,
#' and the squared Cox model z-score. These should be approximately equal under correct implementation.
#'
#' @param dfcount A list or data frame from \code{\link{df_counting}} containing:
#'   \describe{
#'     \item{lr}{Weighted log-rank statistic}
#'     \item{sig2_lr}{Variance of the weighted log-rank statistic}
#'     \item{z.score}{Standardized z-score from Cox model}
#'     \item{logrank_results}{List containing \code{chisq}, the log-rank chi-squared statistic}
#'   }
#' @param verbose Logical; if \code{TRUE}, prints the comparison table to the console.
#'   Default: \code{TRUE}.
#'
#' @return A data frame with one row and three columns, returned invisibly:
#' \describe{
#'   \item{zlr_sq}{Squared standardized weighted log-rank: \eqn{(lr / \sqrt{sig2\_lr})^2}}
#'   \item{logrank_chisq}{Chi-squared statistic from log-rank test}
#'   \item{zCox_sq}{Squared z-score from Cox model: \eqn{z.score^2}}
#' }
#'
#' @details
#' This function serves as a diagnostic tool to verify computational consistency.
#' The three statistics should be numerically equivalent (within rounding error):
#' \deqn{(lr / \sqrt{sig2\_lr})^2 \approx logrank\_chisq \approx z.score^2}
#'
#' Discrepancies between these values may indicate:
#' \itemize{
#'   \item Numerical instability in variance estimation
#'   \item Incorrect weighting scheme application
#'   \item Data processing errors
#' }
#'
#' @note This function is primarily used for package development and validation.
#'   End users typically don't need to call it directly.
#'
#' @examples
#' \dontrun{
#' # After running df_counting
#' library(survival)
#' data(veteran)
#' veteran$treat <- as.numeric(veteran$trt) - 1
#'
#' result <- df_counting(
#'   df = veteran,
#'   tte.name = "time",
#'   event.name = "status",
#'   treat.name = "treat"
#' )
#'
#' # Check consistency of test statistics
#' check_results(result)
#'
#' # Store results without printing
#' stats_comparison <- check_results(result, verbose = FALSE)
#' print(stats_comparison)
#' }
#'
#' # Simple example with constructed data
#' dfcount_example <- list(
#'   lr = 2.5,
#'   sig2_lr = 1.0,
#'   z.score = 2.5,
#'   logrank_results = list(chisq = 6.25)
#' )
#' check_results(dfcount_example)
#'
#' @seealso
#' \code{\link{df_counting}} for generating the input object
#'
#' @family diagnostic_functions
#' @export
check_results <- function(dfcount, verbose = TRUE) {
  # Check required columns/elements
  required_cols <- c("lr", "sig2_lr", "z.score", "logrank_results")
  missing_cols <- setdiff(required_cols, names(dfcount))

  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required elements: %s",
                 paste(missing_cols, collapse = ", ")))
  }

  # Verify logrank_results structure
  if (is.null(dfcount$logrank_results$chisq)) {
    stop("logrank_results must contain a 'chisq' element")
  }

  # Calculate statistics
  zlr_sq <- with(dfcount, (lr^2) / sig2_lr)
  zCox_sq <- with(dfcount, z.score^2)
  logrank_chisq <- dfcount$logrank_results$chisq

  # Create summary data frame
  result <- data.frame(
    zlr_sq = zlr_sq,
    logrank_chisq = logrank_chisq,
    zCox_sq = zCox_sq
  )

  if (verbose) {
    cat("\nTest Statistic Comparison:\n")
    cat("(These should be approximately equal)\n\n")
    print(result, row.names = FALSE)

    # Calculate and display differences
    max_val <- max(result)
    min_val <- min(result)
    rel_diff <- (max_val - min_val) / mean(c(max_val, min_val))

    cat("\nRelative difference: ", format(rel_diff * 100, digits = 4), "%\n", sep = "")

    if (rel_diff > 0.01) {
      warning("Statistics differ by more than 1%. Check for computational issues.",
              call. = FALSE)
    }
  }

  invisible(result)
}

#' Plot Kaplan-Meier curves
#'
#' Plots Kaplan-Meier survival curves for groups in the data.
#'
#' @param df Data frame containing survival data.
#' @param tte.name Name of time-to-event column.
#' @param event.name Name of event indicator column.
#' @param treat.name Name of treatment/group column.
#' @param weights Optional; name of weights column.
#' @param ... Additional arguments passed to plot().
#' @importFrom survival Surv survfit
#' @return Kaplan-Meier fit object (invisible).
#' @export

plot_km <- function(df, tte.name, event.name, treat.name, weights=NULL, ...) {
  safe_run({
    surv_obj <- Surv(df[[tte.name]], df[[event.name]])
    formula <- as.formula(paste("surv_obj ~", treat.name))
    if (!is.null(weights)) {
      km_fit <- survfit(formula, data=df, weights=df[[weights]])
    } else {
      km_fit <- survfit(formula, data=df)
    }
    plot(km_fit, mark.time=TRUE, ...)
    invisible(km_fit)
  })
}

#' Plot weighted Kaplan-Meier curves
#'
#' Plots weighted Kaplan-Meier curves using a custom function.
#'
#' @param dfcount Result object from df_counting.
#' @param ... Additional arguments passed to KM_plot_2sample_weighted_counting.
#' @return None. Plots the curves.
#' @export

plot_weighted_km <- function(dfcount, ...) {
  safe_run({
    KM_plot_2sample_weighted_counting(dfcount, ...
    )
  })
}


#' Extract time, event, and weight data for a group
#'
#' Extracts time, event, and weight vectors for a specified group.
#'
#' @param time Numeric vector of times.
#' @param delta Numeric vector of event indicators (1=event, 0=censored).
#' @param wgt Numeric vector of weights.
#' @param z Numeric vector of group indicators.
#' @param group Value of group to extract (default 1).
#' @return List with U (times), D (events), W (weights).
#' @export

extract_group_data <- function(time, delta, wgt, z, group = 1) {
  list(
    U = time[z == group],
    D = delta[z == group],
    W = wgt[z == group]
  )
}

#' Calculate risk set and event counts at time points
#'
#' Calculates risk set and event counts for a group at specified time points, with variance estimation.
#'
#' @param U Numeric vector of times for group.
#' @param D Numeric vector of event indicators for group.
#' @param W Numeric vector of weights for group.
#' @param at_points Numeric vector of time points.
#' @param draws Number of draws for variance estimation (default 0).
#' @param seedstart Random seed for draws (default 816951).
#' @return List with ybar (risk set counts), nbar (event counts), sig2w_multiplier (variance term).
#' @export

calculate_risk_event_counts <- function(U, D, W, at_points, draws = 0, seedstart = 816951) {
  ybar <- colSums(outer(U, at_points, FUN = ">=") * W)
  nbar <- colSums(outer(U[D == 1], at_points, FUN = "<=") * W[D == 1])
  dN <- diff(c(0, nbar))
  # For un-weighted (all weights equal) return standard variance term
  if(length(unique(W)) == 1){
    # Greenwood
    sig2w_multiplier <- ifelse(ybar > 0 & ybar > dN, dN / (ybar * (ybar-dN)), 0.0)
    # Alternative with dN / (ybar^2)
  }
  if(length(unique(W)) > 1){
    n <- length(U)
    event_mat <- outer(U, at_points, FUN = "<=")
    risk_mat  <- outer(U, at_points, FUN = ">=")
    risk_w <- colSums(risk_mat *  W)
    result <- switch(
      as.integer(draws > 0) + 1,
      {
        # draws == 0
        counting <- colSums(event_mat * (D * W))
        dN_w <- diff(c(0, counting))
        dJ <- ifelse(risk_w == 1, 0, (dN_w - 1) / (risk_w - 1))
        dL <- ifelse(risk_w == 0, 0, dN_w / risk_w)
        h2 <- ifelse(risk_w > 0, (1 / (risk_w)), 0)
        sig2w_multiplier <- (h2 * (dN_w - dL))^2
        sig2w_multiplier
      },
      {
        # draws > 0
        set.seed(seedstart)
        G.draws <- matrix(rnorm(draws * n), ncol = draws)
        counting_star_all <- t(event_mat * W) %*% (D * G.draws)
        dN_star_all <- apply(counting_star_all, 2, function(x) diff(c(0, x)))
        drisk_star <- sweep(dN_star_all, 1, risk_w, "/")
        drisk_star[is.infinite(drisk_star) | is.nan(drisk_star)] <- 0
        sig2w_multiplier <- apply(drisk_star, 1, var)
        sig2w_multiplier
      }
    )
  }
  list(ybar = ybar, nbar = nbar, sig2w_multiplier = sig2w_multiplier)
}

#' Get censoring and event times and their indices
#'
#' Extracts censoring and event times and their indices for a group at specified time points.
#'
#' @param time Numeric vector of times.
#' @param delta Numeric vector of event indicators.
#' @param z Numeric vector of group indicators.
#' @param group Value of group to extract.
#' @param censoring_allmarks Logical; if FALSE, remove events from censored.
#' @param at_points Numeric vector of time points.
#' @return List with cens (censored times), ev (event times), idx_cens, idx_ev, idx_ev_full.
#' @export

get_censoring_and_events <- function(time, delta, z, group, censoring_allmarks, at_points) {
  cens <- time[z == group & delta == 0]
  ev <- sort(unique(time[z == group & delta == 1]))
  if (!censoring_allmarks) cens <- setdiff(cens, ev)
  idx_cens <- match(cens, at_points)
  idx_ev <- match(ev, at_points)
  ev <- c(ev, max(time[z == group]))
  idx_ev_full <- match(ev, at_points)
  list(
    cens = cens,
    ev = ev,
    idx_cens = idx_cens,
    idx_ev = idx_ev,
    idx_ev_full = idx_ev_full
  )
}

#' Get risk set counts at specified risk points
#'
#' Returns risk set counts at specified risk points.
#'
#' @param ybar Numeric vector of risk set counts.
#' @param risk_points Numeric vector of risk points.
#' @param at_points Numeric vector of time points.
#' @return Numeric vector of risk set counts at risk points.
#' @export

get_riskpoints <- function(ybar, risk_points, at_points) {
  ybar[match(risk_points, at_points)]
}
