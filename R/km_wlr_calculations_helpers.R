
#' Format p-value for display
#'
#' Formats a p-value for display, showing "<eps" for small values.
#'
#' @param pval Numeric p-value.
#' @param eps Threshold for small p-values.
#' @param digits Number of digits to display.
#' @return Formatted p-value as character.
#' @export

format_pval <- function(pval, eps = 0.001, digits = 3) {
  if (is.na(pval)) return(NA)
  if (pval < eps) return(paste0("<", eps))
  format(round(pval, digits), nsmall = digits)
}

#' Validate required columns in a data frame
#'
#' Checks that all required columns are present in a data frame.
#'
#' @param df Data frame to check.
#' @param required_cols Character vector of required column names.
#' @return NULL if all columns present, otherwise error.
#' @export

validate_input <- function(df, required_cols) {
  missing <- setdiff(required_cols, names(df))
  if (length(missing) > 0) stop(paste("Missing required columns:", paste(missing, collapse = ", ")))
  invisible(NULL)
}


#' Weighted counting process
#'
#' Computes the weighted count of events up to a specified time.
#'
#' @param x Time point.
#' @param y Vector of event/censoring times.
#' @param w Weights (default 1).
#' @return Weighted count of events up to time x.
#' @export

count_weighted <- function(x, y, w = rep(1, length(y))) {
  sum(w * (y <= x))
}

#' Weighted risk set
#'
#' Computes the weighted number at risk at a specified time.
#'
#' @param x Time point.
#' @param y Vector of event/censoring times.
#' @param w Weights (default 1).
#' @return Weighted number at risk at time x.
#' @export

risk_weighted <- function(x, y, w = rep(1, length(y))) {
  sum(w * (y >= x))
}

#' Kaplan-Meier quantile calculation
#'
#' Calculates the quantile time for a Kaplan-Meier curve: the time at which survival drops to
#' \code{qprob} (\code{qprob = 0.5} is the median; \code{qprob = 0.25} is the time at which
#' 25\% are still event-free, i.e. \code{quantile(survfit(...), probs = 0.75)}).
#'
#' With \code{type = "midpoint"} the rule of \code{survival::quantile.survfit} is used: the first
#' time at which S(t) <= qprob (comparisons within \code{sqrt(.Machine$double.eps)}); when S(t)
#' equals \code{qprob} on an interval, the midpoint of that interval; when that plateau runs to the
#' last time point, the midpoint between its start and the last time point. NA when the curve
#' never reaches \code{qprob}. With \code{type = "min"}, the first time at which S(t) <= qprob.
#' The curve is taken to start at S = 1 at time 0 and to end at the last element of
#' \code{time_points}; values may be NA (they are skipped).
#'
#' @param time_points Vector of time points (non-decreasing).
#' @param survival_probs Vector of survival probabilities at \code{time_points}.
#' @param qprob Survival level (default 0.5): the quantile is the time at which S drops to \code{qprob}.
#' @param type Calculation type (midpoint or min).
#' @return Estimated quantile time.
#' @importFrom stats quantile
#' @export

kmq_calculations <- function(time_points, survival_probs, qprob = 0.5, type = "midpoint") {
  if (type != "midpoint") {
    tq <- suppressWarnings(min(time_points[which(survival_probs <= qprob)]))
    return(if (is.finite(tq)) tq else NA)
  }
  ws_km_quantile_time(time_points, survival_probs, qprob)
}

#' Kaplan-Meier quantile and confidence interval
#'
#' Calculates the quantile and confidence interval for a Kaplan-Meier curve. The quantile is the
#' time at which survival drops to \code{qprob} (see \code{\link{kmq_calculations}}); the limits
#' are the same quantity read from the log-transform pointwise limits of the curve (upper limit
#' capped at 1, NA where S = 0), as \code{survival::quantile.survfit} does with
#' \code{conf.type = "log"}.
#'
#' @param time_points Vector of time points.
#' @param survival_probs Vector of survival probabilities.
#' @param se_probs Standard errors of survival probabilities.
#' @param qprob Survival level (default 0.5): the quantile is the time at which S drops to \code{qprob}.
#' @param type Calculation type (midpoint or min).
#' @param conf_level Confidence level (default 0.95).
#' @return List with quantile and confidence interval.
#' @importFrom stats quantile
#' @export

km_quantile <- function(time_points, survival_probs, se_probs = NULL, qprob = 0.5, type = c("midpoint","min"), conf_level = 0.95) {
  type <- match.arg(type)
  z <- qnorm(1 - (1 - conf_level) / 2)
  qhat <- kmq_calculations(time_points = time_points, survival_probs = survival_probs, qprob = qprob, type = type)
  qhat[is.infinite(qhat)] <- NA
  qhat_lower <- qhat_upper <- NA
  if (!is.null(se_probs)) {
    # log transform for CI; where S = 0 the limits are NA (as in survfit)
    lower_probs <- upper_probs <- rep(NA_real_, length(survival_probs))
    ok <- !is.na(survival_probs) & survival_probs > 0
    lower_probs[ok] <- exp(log(survival_probs[ok]) - z * se_probs[ok] / survival_probs[ok])
    upper_probs[ok] <- pmin(exp(log(survival_probs[ok]) + z * se_probs[ok] / survival_probs[ok]), 1)
    qhat_lower <- kmq_calculations(time_points = time_points, survival_probs = lower_probs, qprob = qprob, type = type)
    qhat_lower[is.infinite(qhat_lower)] <- NA
    qhat_upper <- kmq_calculations(time_points = time_points, survival_probs = upper_probs, qprob = qprob, type = type)
    qhat_upper[is.infinite(qhat_upper)] <- NA
  }
  return(list(qhat = qhat, lower = qhat_lower, upper = qhat_upper))
}

#' Table of KM quantiles for two groups
#'
#' Returns a data frame of quantiles and confidence intervals for two groups.
#'
#' @param time_points Vector of time points.
#' @param surv0 Survival probabilities for group 0.
#' @param se0 Standard errors for group 0.
#' @param surv1 Survival probabilities for group 1.
#' @param se1 Standard errors for group 1.
#' @param arms Group labels.
#' @param qprob Quantile probability.
#' @param type Calculation type.
#' @param conf_level Confidence level.
#' @return Data frame of quantiles and CIs for each group.
#' @export

km_quantile_table <- function(time_points, surv0, se0, surv1, se1, arms = c("treat", "control"), qprob = 0.5, type = c("midpoint","min"), conf_level = 0.95) {
  type <- match.arg(type)
  kmq0 <- km_quantile(time_points = time_points, survival_probs = surv0, se_probs = se0, qprob = qprob, type = type, conf_level = conf_level)
  df0 <- data.frame(group = arms[2], quantile = kmq0$qhat, lower = kmq0$lower, upper = kmq0$upper)
  kmq1 <- km_quantile(time_points = time_points, survival_probs = surv1, se_probs = se1, qprob = qprob, type = type, conf_level = conf_level)
  df1 <- data.frame(group = arms[1], quantile = kmq1$qhat, lower = kmq1$lower, upper = kmq1$upper)
  quantiles_df <- rbind(df1, df0)
  return(quantiles_df)
}


#' Weighted Log-Rank and Difference Estimate at a Specified Time
#'
#' Computes the weighted log-rank statistic, its variance, the difference in survival at a specified time (`tzero`),
#' the variance of the difference, their covariance, and correlation, using flexible time-dependent weights.
#' The weighting scheme is selected via the \code{scheme} argument and is calculated using \code{wt.rg.S}.
#'
#' @param dfcounting List output from \code{df_counting} containing risk sets, event counts, and survival estimates.
#' @param scheme_params Named list with numeric weighting parameters \code{rho}
#'   and \code{gamma} (used for "fh" and "custom_code" schemes).
#' @param tzero Time point at which to evaluate the difference in survival (default: NULL, no
#'   difference computed). The curves and the pooled survival are read at the last grid point
#'   \code{<= tzero}; before the first grid point the difference is 0 with variance 0, and at or after
#'   the last grid point the values at the last grid point are used.
#' @param scheme Character string specifying weighting scheme. One of:
#'   "fh" (Fleming-Harrington), "schemper", "XO", "MB", "custom_time", or "custom_code".
#'
#' @return A list with elements:
#'   \item{lr}{Weighted log-rank test statistic.}
#'   \item{sig2_lr}{Variance of the log-rank statistic.}
#'   \item{dhat}{Difference in survival at \code{tzero}.}
#'   \item{cov_wlr_dhat}{Covariance between log-rank and difference at \code{tzero}.}
#'   \item{sig2_dhat}{Variance of the difference at \code{tzero}.}
#'   \item{cor_wlr_dhat}{Correlation between log-rank and difference at \code{tzero}.}
#'
#' @details
#' The weighting scheme is selected via the \code{scheme} argument and calculated using \code{wt.rg.S}.
#' Supports standard Fleming-Harrington, Schemper, XO (Xu & O'Quigley), MB (Maggir-Burman), custom time-based, and custom code weights.
#'
#' With no subject weights or all subject weights equal the three variances are the formulas used
#' before, and these results are unchanged. When \code{dfcounting} is a \code{\link{df_counting}} result
#' whose subject weights are not all equal, the log-rank statistic, the survival difference and their
#' covariance rest on the same per-subject contributions: \code{sig2.score} is
#' \eqn{\sum_i \psi_i(0)^2}{sum_i psi_i(0)^2} (below; the same as \code{df_counting()$sig2_lr}),
#' \code{sig2_dhat} is \eqn{\sum_i d_i^2}{sum_i d_i^2} and \code{cov_wlr_dhat} is
#' \eqn{\sum_i \{-\psi_i(0)\} d_i}{sum_i (-psi_i(0)) d_i}, where \eqn{d_i} is subject \eqn{i}'s
#' contribution to the survival difference at the grid point read for \code{tzero}:
#' \eqn{d_i = -\hat S_1 a_i \psi^{KM}_i}{d_i = -S1 a_i psiKM_i} in the treatment arm and
#' \eqn{d_i = \hat S_0 a_i \psi^{KM}_i}{d_i = S0 a_i psiKM_i} in the control arm, with
#' \eqn{\psi^{KM}_i}{psiKM_i} the Kaplan-Meier influence contribution of that arm (see
#' \code{\link{KM_diff}}, section "Variance of the Kaplan-Meier estimate"). \code{sig2_dhat} then
#' equals \code{KM_diff(draws = 0)$sig2_dhat} at that grid point. \code{score} and \code{dhat} are unchanged.
#'
#' With a \code{df_counting()} result made with \code{weight_variance = "estimated"}, \code{sig2.score},
#' \code{z.score}, \code{sig2_dhat}, \code{cov_wlr_dhat} and \code{cor_wlr_dhat} carry the estimation of
#' the weights (section "Estimated weights"); \code{sig2_dhat} then equals
#' \code{KM_diff(weight_variance = "estimated")$sig2_dhat} at that grid point.
#'
#' @inheritSection df_counting Variance of the log-rank and Cox score with subject weights
#' @inheritSection df_counting Estimated weights
#' @export

wlr_dhat_estimates <- function(dfcounting,
                               scheme = "fh", scheme_params = list(rho = 0, gamma = 0), tzero = NULL
) {
  # Check for required columns
  required_cols <- c("at_points", "nbar0", "nbar1", "ybar0", "ybar1", "survP")
  if(scheme == "schemper") required_cols <- c(required_cols,"survG")
  missing_cols <- setdiff(required_cols, names(dfcounting))
  if (length(missing_cols) > 0) {
    stop("dfcounting is missing required elements: ", paste(missing_cols, collapse = ", "),
         if ("survG" %in% missing_cols) " (scheme = 'schemper' needs the censoring curve survG at at_points)" else "")
  }

  # exact names only ([[ ]]): no partial matching (e.g. survG must not pick up survG_all)
  at_points <- dfcounting[["at_points"]]
  nbar0 <- dfcounting[["nbar0"]]
  nbar1 <- dfcounting[["nbar1"]]
  ybar0 <- dfcounting[["ybar0"]]
  ybar1 <- dfcounting[["ybar1"]]
  S.pool <- dfcounting[["survP"]]
  G.pool <- dfcounting[["survG"]]

  # weights (the censoring curve is only needed, and only added, when present)
  dfwlr <- data.frame(at_points, S.pool, Ybar = ybar0 + ybar1)
  if (!is.null(G.pool)) dfwlr$G.pool <- G.pool

  wt_rg <- get_validated_weights(dfwlr, scheme, scheme_params)

  S1 <- dfcounting[["surv1"]]
  S0 <- dfcounting[["surv0"]]

  dN0 <- diff(c(0, nbar0))
  dN1 <- diff(c(0, nbar1))
  dN <- dN0 + dN1
  ybar <- ybar0 + ybar1

  if(!is.null(tzero)){
    # values at tzero: those at the last grid point <= tzero (step functions);
    # tzero before the first grid point: dhat = 0 (both curves are 1) and variance 0
    le_tzero <- which(at_points <= tzero)
    if (length(le_tzero) > 0) {
      loc_tzero <- max(le_tzero)
      dhat_tzero <- S1[loc_tzero] - S0[loc_tzero]
      Sp_tzero <- S.pool[loc_tzero]
    } else {
      dhat_tzero <- 0
      Sp_tzero <- 1
    }
  }

  K <- ifelse(ybar > 0, wt_rg * (ybar0 * ybar1) / ybar, 0.0)
  drisk0 <- sum(ifelse(ybar0 > 0, (K / ybar0) * dN0, 0.0))
  drisk1 <- sum(ifelse(ybar1 > 0, (K / ybar1) * dN1, 0.0))
  lr <- drisk0 - drisk1
  h0 <- ifelse(ybar0 == 0, 0, (K^2 / ybar0))
  h1 <- ifelse(ybar1 == 0, 0, (K^2 / ybar1))
  dJ <- ifelse(ybar == 1, 0, (dN - 1) / (ybar - 1))
  dL <- ifelse(ybar == 0, 0, dN / ybar)
  sig2_lr <- sum((h0 + h1) * (1 - dJ) * dL)

  if(!is.null(tzero)){
    w_tzero <- wt_rg * ifelse(at_points <= tzero, 1, 0)
    w_integral_t0 <- sum(w_tzero * (1 - dJ) * dL)
    cov_wlr_dhat <- Sp_tzero * w_integral_t0
    h2 <- ifelse(ybar0 * ybar1 > 0, (ybar / (ybar0 * ybar1)), 0)
    h2 <- h2 * ifelse(at_points <= tzero, 1, 0)
    sig2_dhat <- (Sp_tzero^2) * sum(h2 * (1 - dJ) * dL)
    # Correlation between log-rank statistic and dhat at time tzero
    cor_wlr_dhat <- if (sig2_lr > 0 && sig2_dhat > 0) cov_wlr_dhat / (sqrt(sig2_lr) * sqrt(sig2_dhat)) else NA_real_
  }

  # Subject weights not all equal (a df_counting() result): the three variances rest on the same
  # per-subject contributions. Score: -psi_i(0) (see ws_score_psi()). Survival difference at the grid
  # point read for tzero: d_i = -S1 a_i psiKM_i (treated), +S0 a_i psiKM_i (control), psiKM being the
  # Kaplan-Meier influence of ws_km_psi_parts(). The score and dhat are unchanged.
  subj <- ws_lr_subject_data(dfcounting, at_points)
  if (!is.null(subj)) {
    psi <- ws_score_psi(subj$time, subj$delta, subj$z, subj$w_hat, K = wt_rg[match(subj$time, at_points)], beta = 0)$psi
    ew <- subj$ew
    # Estimated weights: psi_i + IF_i' B for the score, d_i + IF_i' (S0 Dg0 - S1 Dg1) for the difference
    if (!is.null(ew)) psi <- ws_score_ew(psi, subj$w_hat, ew$IF, ew$dw)
    sig2_lr <- sum(psi^2)
    if (!is.null(tzero)) {
      d <- numeric(length(psi))
      if (length(le_tzero) > 0) {
        dDg <- 0
        for (arm in c(0, 1)) {
          ia <- which(subj$z == arm)
          parts <- ws_km_psi_parts(subj$time[ia], subj$delta[ia], at_points,
                                   if (arm == 1) ybar1 else ybar0, if (arm == 1) nbar1 else nbar0)
          psi_km <- ifelse(parts$k <= loc_tzero, parts$exit, -parts$H[loc_tzero])
          sgnS <- if (arm == 1) -S1[loc_tzero] else S0[loc_tzero]
          d[ia] <- sgnS * subj$w_hat[ia] * psi_km
          if (!is.null(ew)) dDg <- dDg + sgnS * colSums(psi_km * ew$dw[ia, , drop = FALSE])
        }
        if (!is.null(ew)) d <- d + as.numeric(ew$IF %*% dDg)
      }
      sig2_dhat <- sum(d^2)
      cov_wlr_dhat <- sum(-psi * d)
      cor_wlr_dhat <- if (sig2_lr > 0 && sig2_dhat > 0) cov_wlr_dhat / (sqrt(sig2_lr) * sqrt(sig2_dhat)) else NA_real_
    }
  }

  if(!is.null(tzero)){
    ans <- list(
      score = lr, sig2.score = sig2_lr, z.score = lr / sqrt(sig2_lr), dhat = dhat_tzero,
      cov_wlr_dhat = cov_wlr_dhat, sig2_dhat = sig2_dhat, cor_wlr_dhat = cor_wlr_dhat
    )
  } else{
    ans <- list(
      score = lr, sig2.score  = sig2_lr, z.score = lr / sqrt(sig2_lr)
    )
  }
  return(ans)
}




#' Kaplan-Meier Survival Estimates and Variance
#'
#' Computes Kaplan-Meier survival estimates and their variances given risk and event counts.
#' Survival is returned in \[0, 1\]; when everything still at risk fails at a time point, survival
#' is exactly 0 from there on and so is its variance.
#'
#' @param ybar Vector of risk set sizes at each time point.
#' @param nbar Vector of cumulative event counts at each time point.
#' @details The result is the Kaplan-Meier estimate only when the time points behind \code{ybar}
#'   and \code{nbar} include every event time of the group (see \code{\link{KM_diff}}, which
#'   arranges this for any \code{at_points}).
#' @param sig2w_multiplier Optional vector for variance calculation. If NULL, calculated internally.
#' @return List with survival estimates and variances.
#' @export

KM_estimates <- function(ybar, nbar, sig2w_multiplier = NULL){
  dN <- diff(c(0, nbar))
  dN_risk <- ifelse(ybar > 0, dN / ybar, 0.0)
  # Everything still at risk fails (dN equals ybar up to relative tolerance 1e-10): the
  # factor is exactly 0, so S is exactly 0 from there on (no floating-point residue).
  all_fail <- ybar > 0 & abs(dN - ybar) <= 1e-10 * ybar
  dN_risk[all_fail] <- 1
  S_KM <- pmin(pmax(cumprod(1 - dN_risk), 0), 1)
  if(is.null(sig2w_multiplier)){
    sig2w_multiplier  <- ifelse(ybar > 0 & ybar > dN, dN / (ybar * (ybar - dN)), 0.0)
  }
  var_KM <- (S_KM^2) * cumsum(sig2w_multiplier)
  var_KM[S_KM == 0] <- 0
  list(S_KM = S_KM, sig2_KM = var_KM)
}

#' Event and Risk Matrices for Survival Analysis
#'
#' Constructs matrices indicating event and risk status for each subject at specified time points.
#'
#' @param U Vector of observed times (e.g., time-to-event).
#' @param at.points Vector of time points at which to evaluate events and risk.
#' @return A list with event and risk matrices.
#' @export

get_event_risk_matrices <- function(U, at.points) {
  event_mat <- outer(U, at.points, FUN = "<=")
  risk_mat  <- outer(U, at.points, FUN = ">=")
  list(event_mat = event_mat, risk_mat = risk_mat)
}



#' Resampling Survival Curves for Confidence Bands
#'
#' Performs resampling to generate survival curves for a group, used for constructing confidence bands.
#'
#' @details Each column is one resampled deviation \eqn{\hat S^*(t) - \hat S(t)} at \code{at.points},
#' for one column of multipliers \code{G_draws} (one row per subject). With no weights or all weights
#' equal, the multipliers perturb the events only:
#' \eqn{-\hat S(t) \sum_{s \le t} \sum_i G_i w_i dN_i(s) / Y_w(s)}; these results are unchanged.
#' With weights that are not all equal, the deviation is built from the per-subject influence
#' contributions \eqn{\psi_i(t)} of the robust variance (see \code{\link{KM_diff}}, section
#' "Variance of the Kaplan-Meier estimate"):
#' \eqn{-\hat S(t) \sum_i G_i w_i \psi_i(t)}, so its variance over draws is the robust variance
#' \eqn{\hat S(t)^2 \sum_i w_i^2 \psi_i(t)^2}, which equals that of
#' \code{survival::survfit(weights = w, robust = TRUE)}. The weights are treated as fixed and known.
#' Both forms need \code{at.points} to hold every event time of the group
#' (\code{\link{KM_diff}} arranges this).
#'
#' @param U Vector of observed times.
#' @param W Vector of weights.
#' @param D Vector of event indicators (0/1).
#' @param at.points Vector of time points for evaluation.
#' @param draws.band Number of resampling draws.
#' @param surv Vector of survival estimates.
#' @param G_draws Matrix of random draws for resampling.
#' @return Matrix of resampled survival curves.
#' @export

resampling_survival <- function(U, W, D, at.points, draws.band, surv, G_draws) {
  if (ws_weights_unequal(W)) {
    # Weights not all equal: per-subject influence contributions, -S(t) * sum_i G_ib w_i psi_i(t)
    # (see ws_km_psi_parts()); their variance over draws is the robust variance of KM_estimates().
    ybar <- colSums(outer(U, at.points, FUN = ">=") * W)
    nbar <- colSums(outer(U[D == 1], at.points, FUN = "<=") * W[D == 1])
    M <- ws_km_ij_sorted(U, D, W, at.points, ybar, nbar, G = G_draws)
    return((-1) * surv * M)
  }
  # No weights or all weights equal: events-only perturbation, as before
  mats <- get_event_risk_matrices(U, at.points)
  risk_w <- colSums(mats$risk_mat * W)
  counting_star_all <- t(mats$event_mat * W) %*% (D * G_draws)
  dN_star_all <- apply(counting_star_all, 2, function(x) diff(c(0, x)))
  drisk_star <- sweep(dN_star_all, 1, risk_w, "/")
  drisk_star[is.infinite(drisk_star) | is.nan(drisk_star)] <- 0
  surv_star <- (-1) * surv * apply(drisk_star, 2, cumsum)
  return(surv_star)
}

#' Kaplan-Meier Difference Between Groups
#'
#' Calculates the difference in Kaplan-Meier curves between two groups, with confidence
#' intervals and optional resampling-based simultaneous confidence bands.
#'
#' @param df Data frame containing survival data.
#' @param tte.name Character; name of time-to-event variable in \code{df}.
#' @param event.name Character; name of event indicator variable in \code{df} (1=event, 0=censored).
#'   Logical (TRUE = event) and 1/2 coding (2 = event, as \code{survival::Surv} reads it) are also accepted.
#' @param treat.name Character; name of treatment group variable in \code{df} (0=control, 1=treatment).
#' @param weight.name Character or NULL; name of weights variable in \code{df}.
#' @param at_points Numeric vector; time points at which the curves are returned: finite and
#'   non-decreasing (duplicates are allowed and returned as given). The Kaplan-Meier value at each
#'   point is right-continuous, whether or not the points include the event times.
#'   Default: the sorted observed times.
#' @param alpha Numeric; significance level for the pointwise intervals and the simultaneous band
#'   (the band's critical value is the \code{1 - alpha} quantile of the resampled suprema). Default: 0.05.
#' @param seedstart Integer; random seed for reproducibility. Default: 8316951.
#' @param draws Integer; number of draws for pointwise variance estimation. Default: 0.
#' @param risk.points Numeric vector; time points for risk table display.
#' @param draws.band Integer; number of draws for simultaneous confidence bands. Default: 0.
#' @param tau.seq Numeric; step size for tau sequence when \code{modify_tau=TRUE}. Default: 0.25.
#' @param qtau Numeric; quantile for tau range restriction. Default: 0.025.
#' @param show_resamples Logical; whether to plot the first 20 standardised resampled difference curves
#'   (only when \code{draws.band > 0}). Default: FALSE.
#' @param modify_tau Logical; whether to restrict time range for simultaneous bands. Default: FALSE.
#' @param weight_variance How the weights are treated in standard errors, intervals, bands and tests.
#'   \code{"known"} (default): as fixed, known numbers. \code{"estimated"}: as estimated by the logistic
#'   propensity model given in \code{weight_model}, whose estimation is carried into the inference. Estimates
#'   never change. See the section "Estimated weights".
#' @param weight_model With \code{weight_variance = "estimated"}, the object returned by
#'   \code{\link{ps_weights}} for the rows of \code{df}, in the same order; its weights must be the column
#'   \code{weight.name}. Ignored with \code{"known"}.
#'
#' @return A list containing:
#' \describe{
#'   \item{at_points}{Time points used in calculations}
#'   \item{surv0, surv1}{Survival estimates for control and treatment groups}
#'   \item{sig2_surv0, sig2_surv1}{Variance estimates for survival curves}
#'   \item{dhat}{Survival difference (S1 - S0) at each time point}
#'   \item{sig2_dhat}{Variance of survival difference}
#'   \item{lower, upper}{Pointwise confidence limits (1 - alpha/2)}
#'   \item{sb_lower, sb_upper}{Simultaneous band limits (if draws.band > 0)}
#'   \item{c_alpha_band}{Critical value for simultaneous band (if draws.band > 0)}
#'   \item{dhat_star}{Matrix of resampled differences (if draws.band > 0)}
#'   \item{Zdhat_star}{Standardized resampled differences (if draws.band > 0)}
#'   \item{weight_variance}{\code{"known"} or \code{"estimated"}: how the weights were treated}
#' }
#'
#' @details
#' This function computes the difference in Kaplan-Meier survival curves, delta(t) = S_1(t) - S_0(t),
#'
#' along with variance estimates and confidence intervals.
#'
#' When \code{draws.band > 0}, simultaneous confidence bands are constructed using
#' the supremum distribution of the standardized difference process. These bands
#' maintain the specified coverage probability across all time points simultaneously.
#'
#' The pointwise variance of each arm's curve is Greenwood's formula for unweighted data (or all
#' weights equal) and the robust variance for weights that are not all equal (see the section below),
#' or resampling-based when \code{draws > 0} and the weights are not all equal.
#'
#' For weights that are not all equal, the resampled curves (\code{surv0_star}, \code{surv1_star},
#' and from them \code{dhat_star}, \code{Zdhat_star}, the simultaneous band and
#' \code{\link{cumulative_rmst_bands}}) are built from the same per-subject contributions as the
#' robust variance (\code{\link{resampling_survival}}), so the variance across draws matches the
#' analytic variance. With no weights or all weights equal they are unchanged.
#'
#' @section Variance of the Kaplan-Meier estimate:
#' With no weights, or with all weights equal, the variance of a Kaplan-Meier estimate is
#' Greenwood's formula, and these results are unchanged. When the weights are not all equal it is
#' the robust (infinitesimal-jackknife) variance
#' \deqn{\widehat{Var}[\hat S(t)] = \hat S(t)^2 \sum_i w_i^2 \psi_i(t)^2,}{Var[S(t)] = S(t)^2 * sum_i w_i^2 psi_i(t)^2,}
#' \deqn{\psi_i(t) = \frac{D_i 1(U_i \le t)}{Y_w(U_i) - dN_w(U_i)} - \sum_{s \le \min(t, U_i)} \frac{dN_w(s)}{Y_w(s)\{Y_w(s) - dN_w(s)\}},}{psi_i(t) = D_i 1(U_i <= t) / (Y_w(U_i) - dN_w(U_i)) - sum_{s <= min(t, U_i)} dN_w(s) / [Y_w(s) (Y_w(s) - dN_w(s))],}
#' where \eqn{U_i}, \eqn{D_i} and \eqn{w_i} are subject \eqn{i}'s time, event indicator and weight,
#' \eqn{Y_w(s)} is the weighted number at risk, \eqn{dN_w(s)} the weighted number of events at
#' \eqn{s}, the sum runs over event times, and a term is 0 where
#' \eqn{Y_w(s) - dN_w(s) \le 10^{-10} Y_w(s)}. It equals the standard error of
#' \code{survival::survfit(weights = w, robust = TRUE)} (on the time-fixed data), reduces to
#' Greenwood's formula at equal weights, and does not depend on the scale of the weights.
#' With \code{draws > 0} and unequal weights the variance is instead the variance over resampling
#' draws of \eqn{\sum_i G_i w_i \psi_i(t)} (standard-normal multipliers \eqn{G_i}), which estimates
#' the same quantity; resampled curves for weighted data are built from the same per-subject
#' contributions (see \code{\link{resampling_survival}}).
#' By default the weights are treated as fixed and known. For weights estimated by a logistic
#' propensity model, \code{weight_variance = "estimated"} carries their estimation into the standard
#' error (section "Estimated weights" of \code{\link{df_counting}}). When the weights are extreme
#' (low effective sample size) the standard error can be too small under either option.
#'
#' @inheritSection df_counting Estimated weights
#'
#' @section Confidence Intervals vs Bands:
#' \itemize{
#'   \item Pointwise CIs (\code{lower}, \code{upper}): Cover the true difference at each time point with probability 1-alpha
#'   \item Simultaneous bands (\code{sb_lower}, \code{sb_upper}): Cover the entire difference curve with probability 1-alpha
#' }
#'
#' @note Treatment must be coded as 0=control, 1=experimental, and both arms must have observations.
#'   Missing values, non-finite or negative times and negative weights are errors (rows are never
#'   dropped silently); zero weights are allowed and those rows contribute nothing: they do not
#'   affect which variance formula runs either (whether the weights are all equal is decided on the
#'   positive weights).
#'
#'   With \code{weight_variance = "estimated"} the variance of the difference is not the sum of the
#'   two arm variances, and \code{draws} must be 0.
#'
#' @examples
#' library(survival)
#' str(veteran)
#' veteran$treat <- as.numeric(veteran$trt) - 1
#'
#' # Basic KM difference
#' result <- KM_diff(
#'   df = veteran,
#'   tte.name = "time",
#'   event.name = "status",
#'   treat.name = "treat"
#' )
#'
#' # Plot the difference
#' plot(result$at_points, result$dhat, type = "s",
#'      xlab = "Time", ylab = "Survival Difference")
#'
#' # With simultaneous confidence bands
#' result_band <- KM_diff(
#'   df = veteran,
#'   tte.name = "time",
#'   event.name = "status",
#'   treat.name = "treat",
#'   draws.band = 1000,
#'   modify_tau = TRUE,
#'   show_resamples = TRUE
#' )
#'
#' @seealso
#' \code{\link{df_counting}} for full survival analysis
#' \code{\link{plotKM.band_subgroups}} for visualization
#' \code{\link{cumulative_rmst_bands}} for RMST analysis
#'
#' @family survival_analysis
#' @family plotting_functions
#' @importFrom survival aeqSurv Surv
#' @importFrom stats quantile
#' @importFrom graphics matplot
#' @export

KM_diff <- function(df, tte.name = "tte", event.name = "event", treat.name = "treat", weight.name=NULL, at_points = sort(df[[tte.name]]), alpha = 0.05, seedstart = 8316951, draws = 0,
                    risk.points = NULL, draws.band = 0, tau.seq = 0.25, qtau = 0.025, show_resamples = FALSE, modify_tau = FALSE,
                    weight_variance = c("known", "estimated"), weight_model = NULL) {

  weight_variance <- match.arg(weight_variance)
  # Columns present; no NA; times finite and >= 0; event 0/1, logical or 1/2; weights >= 0;
  # both arms non-empty (errors name the column and the problem)
  vin <- ws_validate_survdata(df, tte.name, event.name, treat.name, weight.name)
  # "estimated": the pieces of the propensity model (NULL for "known"), rows in the order of df
  ew <- ws_ew_setup(weight_variance, weight_model, df, treat.name, weight.name, draws, "KM_diff")

  tfixed <- aeqSurv(Surv(vin$time, vin$event))
  time<- tfixed[,"time"]
  delta <- tfixed[,"status"]
  z <- df[[treat.name]]
  wgt <- if (!is.null(weight.name)) vin$weight else rep(1, length(time))

  if (!all(z %in% c(0, 1))) stop("Treatment must be numerical indicator: 0=control, 1=experimental")

  if (is.unsorted(time)) {
    ord <- order(time)
    time <- time[ord]
    delta <- delta[ord]
    z <- z[ord]
    wgt <- wgt[ord]
    ew <- ws_ew_reorder(ew, ord)
  }

  # For simultaneous bands restrict time range if modify_tau (otherwise align with RMST 'max(tau)')
  if(draws.band > 0 && modify_tau){
    taus <- quantile(time[delta ==1], c(qtau, 1-qtau))
    at_points<-seq(taus[1],taus[2], by =tau.seq)
    riskp <- risk.points[which(risk.points <= taus[2])]
    at_points <- sort(unique(c(at_points, riskp)))
  }
  if (!is.numeric(at_points) || anyNA(at_points) || any(!is.finite(at_points))) {
    stop("at_points must be numeric and finite.", call. = FALSE)
  }
  if (is.unsorted(at_points)) {
    stop("at_points must be non-decreasing (sorted); duplicates are allowed.", call. = FALSE)
  }

  # Each arm is estimated on its own grid G = (its observed times) U at_points, so that the
  # values returned at at_points are the Kaplan-Meier curve whatever at_points holds.
  # Treatment arm
  group_data <- extract_group_data(time, delta, wgt, z, group = 1)
  km1 <- km_arm_on_grid(group_data$U, group_data$D, group_data$W, at_points, draws = draws, seedstart = seedstart)
  group_data1 <- group_data
  surv1 <- km1$surv
  sig2_surv1 <- km1$sig2

  # Control arm
  group_data <- extract_group_data(time, delta, wgt, z, group = 0)
  km0 <- km_arm_on_grid(group_data$U, group_data$D, group_data$W, at_points, draws = draws, seedstart = seedstart)
  group_data0 <- group_data
  surv0 <- km0$surv
  sig2_surv0 <- km0$sig2

  dhat <- surv1 - surv0
  sig2_dhat <- sig2_surv0 + sig2_surv1

  # Estimated weights: the influence of every subject on each arm's curve includes the term through
  # the propensity model, Phi_i = 1(i in arm) w_i psi_i + IF_i' Dg (see ws_km_ew_group()). The two arms
  # share the model's coefficients, so the variance of the difference is not the sum of the arm variances.
  if (!is.null(ew)) {
    i1 <- which(z == 1); i0 <- which(z == 0)
    arm_ew <- function(gd, rows, km) {
      re <- calculate_risk_event_counts(gd$U, gd$D, gd$W, at_points = km$grid)
      ws_km_ew_group(gd$U, gd$D, gd$W, rows, km$grid, re$ybar, re$nbar, ew)
    }
    e1 <- arm_ew(group_data1, i1, km1)
    e0 <- arm_ew(group_data0, i0, km0)
    sig2_surv1 <- surv1^2 * e1$sumsq[km1$idx]; sig2_surv1[surv1 == 0] <- 0
    sig2_surv0 <- surv0^2 * e0$sumsq[km0$idx]; sig2_surv0[surv0 == 0] <- 0
    Dg1 <- e1$Dg[km1$idx, , drop = FALSE]; Dg0 <- e0$Dg[km0$idx, , drop = FALSE]
    cross <- rowSums(e1$A[km1$idx, , drop = FALSE] * Dg0) + rowSums(e0$A[km0$idx, , drop = FALSE] * Dg1) +
      rowSums((Dg1 %*% ew$CIF) * Dg0)
    sig2_dhat <- pmax(sig2_surv1 + sig2_surv0 - 2 * surv1 * surv0 * cross, 0)
  }

  surv0_star <- surv1_star <- dhat_star <- NULL
  c_alpha_band <- sb_lower <- sb_upper <- NULL
  Zdhat_star <- NULL
  if(draws.band > 0){
    # draws > 0: drawn under set.seed(seedstart); the caller's RNG state is restored afterwards
    n0 <- length(group_data0$U)
    n1 <- length(group_data1$U)
    G_draws <- ws_with_seed(seedstart, list(
      G0 = matrix(rnorm(draws.band * n0), ncol = draws.band),
      G1 = matrix(rnorm(draws.band * n1), ncol = draws.band)
    ))
    G0.draws <- G_draws$G0
    G1.draws <- G_draws$G1

    if (is.null(ew)) {
    # Control group (resampled on the arm's grid, read off at at_points)
    surv0_star <- resampling_survival(group_data0$U, group_data0$W, group_data0$D, km0$grid, draws.band, km0$S_grid, G0.draws)[km0$idx, , drop = FALSE]
    # Treatment group
    surv1_star <- resampling_survival(group_data1$U, group_data1$W, group_data1$D, km1$grid, draws.band, km1$S_grid, G1.draws)[km1$idx, , drop = FALSE]
    } else {
      # Estimated weights: -S_a(t) sum_i G_i Phi^a_i(t), with G the multipliers of all subjects (row i is
      # that subject's row of G0.draws or G1.draws): the influence-based deviation of the arm minus
      # (S_a Dg_a) (IF' G).
      IFG <- crossprod(ew$IF[i0, , drop = FALSE], G0.draws) + crossprod(ew$IF[i1, , drop = FALSE], G1.draws)
      arm_star <- function(gd, km, e, Ga) {
        M <- ws_km_ij_draws(e$parts, gd$W, Ga)
        ((-1) * km$S_grid * M - (km$S_grid * e$Dg) %*% IFG)[km$idx, , drop = FALSE]
      }
      surv0_star <- arm_star(group_data0, km0, e0, G0.draws)
      surv1_star <- arm_star(group_data1, km1, e1, G1.draws)
    }

    dhat_star <- (surv1_star - surv0_star)

    Zdhat_star <- dhat_star / sqrt(sig2_dhat)

    # simultaneous band
    sups <- apply(abs(Zdhat_star), 2, max, na.rm = TRUE)
    c_alpha_band <- quantile(sups, 1 - alpha)
    # Show first 20 (or all of them when draws.band < 20)
    if(show_resamples){
      matplot(at_points, Zdhat_star[, seq_len(min(20, draws.band)), drop = FALSE], type="s", lty=2, lwd = 1, xlab="time",
              ylab = sprintf("Normalized survival differences (1st %d)", min(20, draws.band)),
              main = sprintf("c_alpha (simul. band): %.2f", c_alpha_band)
      )
    }
    # simulataneous band
    sb_lower <- dhat - c_alpha_band * sqrt(sig2_dhat)
    sb_upper <- dhat + c_alpha_band * sqrt(sig2_dhat)
  }
  # Standard point-wise CIs
  c_alpha <- qnorm(1 - alpha / 2)
  lower <- dhat - c_alpha * sqrt(sig2_dhat)
  upper <- dhat + c_alpha * sqrt(sig2_dhat)
  list(
    at_points = at_points, surv0 = surv0, sig2_surv0 = sig2_surv0,
    surv1 = surv1, sig2_surv1 = sig2_surv1, dhat = dhat, sig2_dhat = sig2_dhat,
    lower = lower, upper = upper,
    dhat_star = dhat_star,
    Zdhat_star = Zdhat_star,
    surv0_star = surv0_star, surv1_star = surv1_star,
    c_alpha_band = c_alpha_band, sb_lower = sb_lower, sb_upper = sb_upper,
    weight_variance = weight_variance
  )
}


#' Weighted log-rank cumulative statistics
#'
#' Calculates cumulative weighted log-rank statistics for survival data.
#'
#' @param df_weights Data frame with the time grid, risk and event counts and pooled survival
#'   (\code{at_points}, \code{nbar0}, \code{nbar1}, \code{ybar0}, \code{ybar1}, \code{survP}, and
#'   \code{survG} for Schemper weights), or a \code{\link{df_counting}} result, which holds the same
#'   elements together with the subject-level data (\code{time}, \code{delta}, \code{z}, \code{w_hat}).
#' @param scheme Weighting scheme.
#' @param scheme_params List of scheme parameters.
#' @param return_cumulative Logical; whether to return cumulative statistics.
#' @return List with \code{z.score}, \code{time} (the grid), \code{score} and \code{sig2.score} (the
#'   variance of the score); with \code{return_cumulative = TRUE} each is a vector along the grid.
#' @details With no subject weights, with all subject weights equal, and for a bare table of counts
#'   (a data frame without the subject-level data) the variance of the score is the formula used
#'   before, computed from the counts; these results are unchanged. When \code{df_weights} is a
#'   \code{\link{df_counting}} result whose subject weights are not all equal (and inside
#'   \code{df_counting()} itself, unstratified and per stratum), \code{sig2.score} is the sum of squared
#'   per-subject score contributions \eqn{\sum_i \psi_i(0)^2}{sum_i psi_i(0)^2} described below, and
#'   \code{z.score} follows; \code{score} is unchanged. With \code{return_cumulative = TRUE} the variance
#'   at grid time \eqn{t} is the sum of squared contributions to the score up to \eqn{t}: a subject with
#'   \eqn{U_i \le t}{U_i <= t} contributes \eqn{\psi_i(0)}{psi_i(0)}, a subject still at risk contributes
#'   \eqn{-a_i \{Z_i c_K(t) - c_{Kz}(t)\}}{-a_i (Z_i cK(t) - cKz(t))}, with
#'   \eqn{c_K(t) = \sum_{s \le t} K(s) d\hat\Lambda(s)}{cK(t) = sum_{s <= t} K(s) dL(s)} and
#'   \eqn{c_{Kz}(t) = \sum_{s \le t} K(s) \bar z(s) d\hat\Lambda(s)}{cKz(t) = sum_{s <= t} K(s) zbar(s) dL(s)};
#'   its last value is the non-cumulative variance. A table of counts taken from a weighted analysis
#'   does not carry the subject-level data, so it gets the earlier formula: pass the
#'   \code{df_counting()} result itself to get the variance described below. With a
#'   \code{df_counting()} result made with \code{weight_variance = "estimated"}, \code{sig2.score} and
#'   \code{z.score}, cumulative or not, carry the estimation of the weights (section "Estimated
#'   weights"); on a bare table of counts nothing changes.
#' @inheritSection df_counting Variance of the log-rank and Cox score with subject weights
#' @inheritSection df_counting Estimated weights
#' @export

wlr_cumulative <- function(df_weights, scheme, scheme_params = list(rho = 0, gamma = 0), return_cumulative = FALSE) {
# Check for required columns
  required_cols <- c("at_points", "nbar0", "nbar1", "ybar0", "ybar1", "survP")
  if(scheme == "schemper") required_cols <- c(required_cols,"survG")
  missing_cols <- setdiff(required_cols, names(df_weights))
  if (length(missing_cols) > 0) {
    stop("df_weights is missing required columns: ", paste(missing_cols, collapse = ", "),
         if ("survG" %in% missing_cols) " (scheme = 'schemper' needs the censoring curve survG at at_points)" else "")
  }

  # exact names only ([[ ]]): no partial matching
  at_points <- df_weights[["at_points"]]
  nbar0 <- df_weights[["nbar0"]]
  nbar1 <- df_weights[["nbar1"]]
  ybar0 <- df_weights[["ybar0"]]
  ybar1 <- df_weights[["ybar1"]]
  S.pool <- df_weights[["survP"]]
  G.pool <- df_weights[["survG"]]

  # weights (the censoring curve is only needed, and only added, when present)
  dfwlr <- data.frame(at_points, S.pool, Ybar = ybar0 + ybar1)
  if (!is.null(G.pool)) dfwlr$G.pool <- G.pool
  wt_rg <- get_validated_weights(dfwlr, scheme, scheme_params)

  dN0 <- diff(c(0, nbar0))
  dN1 <- diff(c(0, nbar1))
  dN <- dN0 + dN1
  ybar <- ybar0 + ybar1
  K <- ifelse(ybar > 0, wt_rg * (ybar0 * ybar1) / ybar, 0.0)
  h0 <- ifelse(ybar0 == 0, 0, (K^2 / ybar0))
  h1 <- ifelse(ybar1 == 0, 0, (K^2 / ybar1))
  dJ <- ifelse(ybar == 1, 0, (dN - 1) / (ybar - 1))
  dL <- ifelse(ybar == 0, 0, dN / ybar)

  if(return_cumulative){
    drisk0 <- cumsum(ifelse(ybar0 > 0, (K / ybar0) * dN0, 0.0))
    drisk1 <- cumsum(ifelse(ybar1 > 0, (K / ybar1) * dN1, 0.0))
    sig2_lr <- cumsum((h0 + h1) * (1 - dJ) * dL)
  } else {
    drisk0 <- sum(ifelse(ybar0 > 0, (K / ybar0) * dN0, 0.0))
    drisk1 <- sum(ifelse(ybar1 > 0, (K / ybar1) * dN1, 0.0))
    sig2_lr <- sum((h0 + h1) * (1 - dJ) * dL)
  }
  lr <- drisk0 - drisk1
  # Subject weights that are not all equal, subject-level data present: the variance of the score
  # is the sum of squared per-subject score contributions (see ws_score_psi()); the score is unchanged.
  subj <- ws_lr_subject_data(df_weights, at_points)
  if (!is.null(subj)) {
    sp <- ws_score_psi(subj$time, subj$delta, subj$z, subj$w_hat, K = wt_rg[match(subj$time, at_points)], beta = 0)
    sig2_lr <- if (return_cumulative) ws_lr_cumvar(sp, subj$time, subj$z, subj$w_hat, at_points) else sum(sp$psi^2)
    # Estimated weights (a df_counting() result made with weight_variance = "estimated"): the contributions
    # carry the term through the propensity model, psi_i + IF_i' B (see ws_score_ew(), ws_lr_cumvar_ew()).
    if (!is.null(subj$ew)) {
      sig2_lr <- if (return_cumulative) ws_lr_cumvar_ew(sp, subj$time, subj$z, subj$w_hat, at_points, subj$ew$IF, subj$ew$dw)
                 else sum(ws_score_ew(sp$psi, subj$w_hat, subj$ew$IF, subj$ew$dw)^2)
    }
  }
  z <- lr / sqrt(sig2_lr)
  list(z.score = z, time = at_points, score = lr, sig2.score = sig2_lr)
}

