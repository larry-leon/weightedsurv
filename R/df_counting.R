#' Weighted and Stratified Survival Analysis
#'
#' Performs comprehensive weighted and/or stratified survival analysis, including
#' Cox proportional hazards model, logrank/Fleming-Harrington tests, and calculation
#' of risk/event sets, Kaplan-Meier curves, quantiles, and variance estimates.
#'
#' @param df Data frame containing survival data.
#' @param tte.name Character; name of the time-to-event column in \code{df}.
#' @param event.name Character; name of the event indicator column in \code{df} (1=event, 0=censored).
#'   Logical (TRUE = event) and 1/2 coding (2 = event, as \code{survival::Surv} reads it) are also accepted.
#' @param treat.name Character; name of the treatment/group column in \code{df} (must be coded as 0=control, 1=experimental).
#' @param weight.name Character or NULL; name of the weights column in \code{df}. If NULL, equal weights are used.
#' @param strata.name Character or NULL; name of the strata column in \code{df} for stratified analysis.
#' @param arms Character vector of length 2; group labels. Default: \code{c("treat","control")}.
#' @param time.zero Numeric; time value to use as zero. Default: 0.
#' @param tpoints.add Numeric vector; additional time points to include in calculations. Default: \code{c(0)}.
#' @param by.risk Numeric; interval for risk set time points. Default: 6.
#' @param time.zero.label Numeric; label for time zero in output. Default: 0.0.
#' @param risk.add Numeric vector or NULL; additional specific risk points to include.
#' @param get.cox Logical; whether to fit Cox proportional hazards model. Default: TRUE.
#' @param cox.digits Integer; number of decimal places for Cox output formatting. Default: 2.
#' @param lr.digits Integer; number of decimal places for logrank output formatting. Default: 2.
#' @param cox.eps Numeric; threshold for Cox p-value formatting (values below shown as "<eps"). Default: 0.001.
#' @param lr.eps Numeric; threshold for logrank p-value formatting. Default: 0.001.
#' @param verbose Logical; whether to print warnings and diagnostic messages. Default: FALSE.
#' @param qprob Numeric in (0,1); survival level for the KM quantile table: the quantile reported is the
#'   time at which survival drops to \code{qprob}, i.e. \code{quantile(survfit(...), probs = 1 - qprob)}.
#'   Default: 0.5 (median).
#' @param scheme Character; weighting scheme for logrank/Fleming-Harrington test.
#'   Options: "fh", "schemper", "XO", "MB", "custom_time", "fh_exp1", "fh_exp2". Default: "fh".
#' @param scheme_params List; parameters for the selected weighting scheme. Default: \code{list(rho = 0, gamma = 0)}.
#'   \itemize{
#'     \item For "fh": \code{rho} and \code{gamma} (Fleming-Harrington parameters)
#'     \item For "MB": \code{mb_tstar} (cutoff time)
#'     \item For "custom_time": \code{t.tau}, \code{w0.tau}, \code{w1.tau}
#'   }
#' @param conf_level Numeric in (0,1); confidence level for quantile intervals. Default: 0.95.
#' @param check.KM Logical; whether to check KM curve validity against \code{survival::survfit}. Default: TRUE.
#' @param check.seKM Logical; whether to check KM standard error estimates. Default: FALSE.
#' @param draws Integer; number of draws for resampling-based variance estimation. Default: 0 (no resampling).
#' @param seedstart Integer; random seed for reproducible resampling. Default: 8316951.
#' @param stop.onerror Logical; whether to stop execution on errors (TRUE) or issue warnings (FALSE). Default: FALSE.
#' @param censoring_allmarks Logical; if FALSE, removes events from censored time points. Default: TRUE.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{cox_results}{List with Cox model results including hazard ratio, confidence interval, p-value, and formatted text}
#'   \item{logrank_results}{List with log-rank test chi-square statistic, p-value, and formatted text}
#'   \item{z.score}{Standardized weighted log-rank test statistic}
#'   \item{at_points}{Vector of all time points used in calculations}
#'   \item{surv0, surv1}{Kaplan-Meier survival estimates for control and treatment groups}
#'   \item{sig2_surv0, sig2_surv1}{Variance estimates for survival curves}
#'   \item{survP}{Pooled survival estimates}
#'   \item{survG}{Censoring distribution estimates}
#'   \item{quantile_results}{Data frame with median survival and confidence intervals by group}
#'   \item{lr, sig2_lr}{Weighted log-rank statistic and its variance}
#'   \item{riskpoints0, riskpoints1}{Risk set counts at specified time points}
#'   \item{z.score_stratified}{Stratified z-score (if stratified analysis)}
#' }
#'
#' @details
#' This function implements a comprehensive survival analysis framework supporting:
#' \itemize{
#'   \item Weighted observations via \code{weight.name}
#'   \item Stratified analysis via \code{strata.name}
#'   \item Multiple weighting schemes for log-rank tests
#'   \item Resampling-based variance estimation
#'   \item Automatic validation against \code{survival} package results
#' }
#'
#' The function performs time-fixing using \code{survival::aeqSurv} to handle tied event times.
#'
#' Input checks: missing values (in any named column, including \code{strata.name}), non-finite or
#' negative times, negative weights and an empty arm are errors; rows are never dropped silently.
#' Zero weights are allowed: those rows contribute nothing, and the Cox model and log-rank test are
#' fitted on the positive-weight rows. A stratum in which one arm is empty contributes nothing to the
#' stratified statistics (that arm's stratum curve is \code{NA}). On degenerate data (an arm without events, no events at all, one subject
#' per arm) a result is still returned; quantities that cannot be computed are \code{NA}, with one
#' warning that says why.
#' For stratified analyses, stratum-specific estimates are computed and combined appropriately.
#'
#' @section Variance of the Kaplan-Meier estimate:
#' Unless resampling is requested (\code{draws > 0} with weights that are not all equal), the
#' variance of each arm's Kaplan-Meier estimate is weighted Greenwood:
#' \deqn{\widehat{Var}[\hat S(t)] = \hat S(t)^2 \sum_{s \le t} \frac{\sum_i w_i^2 dN_i(s)}{Y_w(s)\{Y_w(s) - dN_w(s)\}},}{Var[S(t)] = S(t)^2 * sum_{s <= t} [sum_i w_i^2 dN_i(s)] / [Y_w(s) (Y_w(s) - dN_w(s))],}
#' with \eqn{Y_w(s)} the weighted number at risk and \eqn{dN_w(s) = \sum_i w_i dN_i(s)} the
#' weighted number of events at \eqn{s}; a term is 0 where \eqn{Y_w(s) = 0} or
#' \eqn{Y_w(s) - dN_w(s) \le 10^{-10} Y_w(s)}. With all weights equal (or no weights) this is
#' Greenwood's formula. The weights are treated as fixed and known: there is no adjustment for
#' estimated weights (for example estimated propensity scores). The standard error does not depend
#' on the scale of the weights, so weights are not frequency weights: a weight of 2 is not two
#' copies of the row.
#'
#' @section Weighting Schemes:
#' \describe{
#'   \item{fh}{Fleming-Harrington: w(t) = S(t)^rho * (1-S(t))^gamma}
#'   \item{MB}{Magirr-Burman: w(t) = 1/max(S(t), S(t*))}
#'   \item{schemper}{Schemper: w(t) = S(t)/G(t) where G is the censoring distribution}
#'   \item{XO}{Xu-O'Quigley: w(t) = S(t)/Y(t) where Y is risk set size}
#' }
#'
#' @examples
#' # Basic survival analysis
#' library(survival)
#' str(veteran)
#' veteran$treat <- as.numeric(veteran$trt) - 1
#'
#' result <- df_counting(
#'   df = veteran,
#'   tte.name = "time",
#'   event.name = "status",
#'   treat.name = "treat",
#'   arms = c("Treatment", "Control")
#' )
#'
#' # Print results
#' print(result$cox_results$cox_text)
#' print(result$zlogrank_text)
#'
#' # Fleming-Harrington (0,1) weights (emphasizing late differences)
#' result_fh <- df_counting(
#'   df = veteran,
#'   tte.name = "time",
#'   event.name = "status",
#'   treat.name = "treat",
#'   scheme = "fh",
#'   scheme_params = list(rho = 0, gamma = 1)
#' )
#'
#' # Stratified analysis
#' result_strat <- df_counting(
#'   df = veteran,
#'   tte.name = "time",
#'   event.name = "status",
#'   treat.name = "treat",
#'   strata.name = "celltype"
#' )
#'
#' @seealso
#' \code{\link[survival]{coxph}}, \code{\link[survival]{survdiff}}, \code{\link[survival]{survfit}}
#' \code{\link{cox_rhogamma}} for weighted Cox models
#' \code{\link{KM_diff}} for Kaplan-Meier differences
#'
#' @references
#' Fleming, T. R. and Harrington, D. P. (1991). Counting Processes and Survival Analysis. Wiley.
#'
#' Magirr, D. and Burman, C. F. (2019). Modestly weighted logrank tests.
#' Statistics in Medicine, 38(20), 3782-3790.
#'
#' @family survival_analysis
#' @family weighted_tests
#' @importFrom stats as.formula pchisq
#' @importFrom survival aeqSurv Surv coxph survfit survdiff strata
#' @export

df_counting <- function(df, tte.name = "tte", event.name = "event", treat.name = "treat", weight.name=NULL, strata.name = NULL, arms=c("treat","control"),
                        time.zero=0, tpoints.add=c(0),
                        by.risk=6, time.zero.label = 0.0, risk.add=NULL, get.cox=TRUE, cox.digits=2, lr.digits=2,
                        cox.eps = 0.001, lr.eps = 0.001, verbose = FALSE,
                        qprob=0.5,
                        scheme = "fh", scheme_params = list(rho = 0, gamma = 0),
                        conf_level = 0.95, check.KM = TRUE, check.seKM = FALSE, draws = 0, seedstart = 8316951,
                        stop.onerror=FALSE,censoring_allmarks=TRUE) {

  # Columns present (strata too); no NA; times finite and >= 0; event 0/1, logical or 1/2;
  # weights >= 0; both arms non-empty. Errors name the column and the problem.
  vin <- ws_validate_survdata(df, tte.name, event.name, treat.name, weight.name, strata.name)
  # 1/2 event coding: use the 0/1 recoding everywhere, so that subsets cannot be misread
  if (any(vin$event != df[[event.name]])) df[[event.name]] <- vin$event
  # Zero weights are accepted: those rows contribute nothing. Cox and log-rank are fitted on the
  # positive-weight rows; per-observation output keeps every row.
  pos_w <- if (is.null(weight.name)) rep(TRUE, nrow(df)) else vin$weight > 0
  df_fit <- if (all(pos_w)) df else df[which(pos_w), ]
  events_by_arm <- c(sum(vin$event[vin$treat == 0 & pos_w]), sum(vin$event[vin$treat == 1 & pos_w]))
  # Quantities that cannot be computed on degenerate data are set to NA; the reasons are
  # collected here and reported in one warning at the end.
  na_reasons <- character(0)

  # Validate scheme and parameters
  supported_schemes <- c("fh", "schemper", "XO", "MB", "custom_time", "fh_exp1", "fh_exp2")
  if (!(scheme %in% supported_schemes)) {
    stop("scheme must be one of: ", paste(supported_schemes, collapse = ", "))
  }

  ans <- list()

  cox_results <- NULL
  cox_na <- list(cox_fit = NA, hr = NA_real_, hr_ci = c(NA_real_, NA_real_), pval = NA_real_, score = NA_real_,
                 cox_text = "HR = NA")
  if (get.cox && any(events_by_arm == 0)) {
    na_reasons <- c(na_reasons, sprintf("Cox model: no events in the %s; the hazard ratio is not identified",
                                        paste(c("control arm", "treatment arm")[events_by_arm == 0], collapse = " and ")))
    cox_results <- cox_na
  } else if (get.cox) {
    treat.name.strata <- treat.name
    if (!is.null(strata.name)) {
      treat.name.strata <- paste(treat.name,"+",paste("strata(",eval(strata.name),")"))
    }
    cox_formula <- as.formula(paste0("survival::Surv(", tte.name, ",", event.name, ") ~ ", treat.name.strata))
    cox_fit <- tryCatch({
      if (!is.null(weight.name)) {
        survival::coxph(cox_formula, data = df_fit, weights = df_fit[[weight.name]], robust = TRUE)
      } else {
        survival::coxph(cox_formula, data = df_fit, robust=TRUE)
      }
    }, error = function(e) {
      na_reasons <<- c(na_reasons, paste("Cox model failed:", conditionMessage(e)))
      NULL
    })
  }
  if (get.cox && is.null(cox_results) && is.null(cox_fit)) cox_results <- cox_na
  if (get.cox && is.null(cox_results)) {
    cox_summary <- summary(cox_fit)
    cox_score <- cox_summary$sctest[1]
    hr <- exp(cox_fit$coef)
    hr_ci <- exp(confint(cox_fit))
    pval <- cox_summary$coefficients[1, "Pr(>|z|)"]
    cox_text <- paste0("HR = ", round(hr, cox.digits),
                        " (", round(hr_ci[1], cox.digits), ", ", round(hr_ci[2], cox.digits), ")")

     cox_results <- list(
      cox_fit = cox_fit,
      hr = hr,
      hr_ci = hr_ci,
      pval = pval,
      score = cox_score,
      cox_text = cox_text
    )
  }
  ans$cox_results <- cox_results
  logrank_results <- NULL
  surv_obj <- survival::Surv(df_fit[[tte.name]], df_fit[[event.name]])
  group <- df_fit[[treat.name]]
 # survdiff does not support case weights

  lr_rho <- 0
  if(scheme == "fh") lr_rho <- scheme_params$rho

  if (sum(events_by_arm) == 0) {
    na_reasons <- c(na_reasons, "log-rank test (survdiff): no events")
    logrank_fit <- list(chisq = NA_real_)
  } else {
    logrank_fit <- tryCatch({
      if (!is.null(strata.name)) {
        strata_var <- df_fit[[strata.name]]
        survival::survdiff(surv_obj ~ group + strata(strata_var), rho = lr_rho)
      } else {
        survival::survdiff(surv_obj ~ group, rho = lr_rho)
      }
    }, error = function(e) {
      na_reasons <<- c(na_reasons, paste("log-rank test (survdiff) failed:", conditionMessage(e)))
      list(chisq = NA_real_)
    })
  }

  chisq <- logrank_fit$chisq
  pval <- 1 - pchisq(chisq, df = 1)
  logrank_text <- paste0("Logrank p = ", format_pval(pval, eps = lr.eps, digits = lr.digits))
  logrank_results <- list(
    chisq = chisq,
    pval = pval,
    logrank_text = logrank_text
  )

  ans$logrank_results <- logrank_results

  # Implement timefix per survival package
  tfixed <- aeqSurv(Surv(df[[tte.name]],df[[event.name]]))
  time<- tfixed[,"time"]
  delta <- tfixed[,"status"]
  z <- df[[treat.name]]

  strata <- if (!is.null(strata.name)) df[[strata.name]] else rep("All", length(time))
  wgt <- if (!is.null(weight.name)) df[[weight.name]] else rep(1, length(time))

 if (!all(z %in% c(0, 1))) stop("Treatment must be numerical indicator: 0=control, 1=experimental")

  if (is.unsorted(time)) {
    ord <- order(time)
    time <- time[ord]
    delta <- delta[ord]
    z <- z[ord]
    strata <- strata[ord]
    wgt <- wgt[ord]
  }

  ans$time <- time
  ans$delta <- delta
  ans$z <- z
  ans$strata <- strata
  ans$w_hat <- wgt

  at_points_all <- time
  ans$at_points_all <- at_points_all

  # Pooled KM estimates at all timepoints for input to cox estimation
  risk_event <- calculate_risk_event_counts(time, delta, wgt, at_points_all)
  ans$survP_all <- KM_estimates(ybar = risk_event$ybar, nbar = risk_event$nbar, sig2w_multiplier = risk_event$sig2w_multiplier)$S_KM
  ans$ybar_all <- risk_event$ybar
  ans$nbar_all <- risk_event$nbar
  rm("risk_event")
  # Also the censoring distribution
  risk_event <- calculate_risk_event_counts(time, 1-delta, wgt, at_points_all)
  ans$survG_all <- KM_estimates(ybar = risk_event$ybar, nbar = risk_event$nbar, sig2w_multiplier = risk_event$sig2w_multiplier)$S_KM
  rm("risk_event")



  stratum <- unique(strata)
  risk.points <- sort(unique(c(seq(time.zero.label, max(time), by = ifelse(is.null(by.risk), 1, by.risk)), risk.add)))
  ans$risk.points <- risk.points
  ans$risk.points.label <- as.character(c(time.zero.label, risk.points[-1]))

  at_points <- sort(unique(c(time, time.zero, tpoints.add, risk.points)))


  # Treatment arm
  group_data <- extract_group_data(time, delta, wgt, z, group = 1)
  risk_event <- calculate_risk_event_counts(group_data$U, group_data$D, group_data$W, at_points, draws, seedstart)
  cens_ev <- get_censoring_and_events(time, delta, z, 1, censoring_allmarks, at_points)
  riskpoints1 <- get_riskpoints(risk_event$ybar, risk.points, at_points)
  temp <- KM_estimates(ybar = risk_event$ybar, nbar = risk_event$nbar, sig2w_multiplier = risk_event$sig2w_multiplier)
  surv1 <- temp$S_KM
  sig2_surv1 <- temp$sig2_KM

  nbar1 <- risk_event$nbar
  ybar1 <- risk_event$ybar

  # Store in ans
  ans$idx1 <- cens_ev$idx_cens
  ans$idv1 <- cens_ev$idx_ev_full
  idv1.check <- cens_ev$idx_ev
  ans$idv1.check <- idv1.check
  ans$ev1 <- cens_ev$ev
  ans$cens1 <- cens_ev$cens
  ans$riskpoints1 <- riskpoints1

  # Control arm
  group_data <- extract_group_data(time, delta, wgt, z, group = 0)
  risk_event <- calculate_risk_event_counts(group_data$U, group_data$D, group_data$W, at_points, draws, seedstart)
  cens_ev <- get_censoring_and_events(time, delta, z, 0, censoring_allmarks, at_points)
  riskpoints0 <- get_riskpoints(risk_event$ybar, risk.points, at_points)
  temp <- KM_estimates(ybar = risk_event$ybar, nbar = risk_event$nbar, sig2w_multiplier = risk_event$sig2w_multiplier)
  surv0 <- temp$S_KM
  sig2_surv0 <- temp$sig2_KM

  nbar0 <- risk_event$nbar
  ybar0 <- risk_event$ybar

  # Store in ans
  ans$idx0 <- cens_ev$idx_cens
  ans$idv0 <- cens_ev$idx_ev_full
  idv0.check <- cens_ev$idx_ev
  ans$idv0.check <- idv0.check
  ans$ev0 <- cens_ev$ev
  ans$cens0 <- cens_ev$cens
  ans$riskpoints0 <- riskpoints0

  # Quantiles: the time at which S drops to qprob (survfit's quantile(probs = 1 - qprob)).
  # Each arm's curve is read up to that arm's last observed time (positive-weight rows), so that
  # a plateau at the level that runs to the end of follow-up is handled as survfit does.
  arm_q <- function(a, S, sig2) {
    last_a <- max(time[z == a & wgt > 0])
    sel <- at_points <= last_a
    km_quantile(at_points[sel], S[sel], sqrt(sig2[sel]), qprob = qprob, type = "midpoint", conf_level = conf_level)
  }
  kmq1 <- arm_q(1, surv1, sig2_surv1)
  kmq0 <- arm_q(0, surv0, sig2_surv0)
  get_kmq <- rbind(data.frame(group = arms[1], quantile = kmq1$qhat, lower = kmq1$lower, upper = kmq1$upper),
                   data.frame(group = arms[2], quantile = kmq0$qhat, lower = kmq0$lower, upper = kmq0$upper))
  ans$quantile_results <- get_kmq


  # Pooled KM estimates
  risk_event <- calculate_risk_event_counts(time, delta, wgt, at_points)
  temp <- KM_estimates(ybar = risk_event$ybar, nbar = risk_event$nbar, sig2w_multiplier = risk_event$sig2w_multiplier)
  survP <- temp$S_KM
  sig2_survP <- temp$sig2_KM

  # Pooled Censoring estimates
  risk_event <- calculate_risk_event_counts(time, 1- delta, wgt, at_points)
  temp <- KM_estimates(ybar = risk_event$ybar, nbar = risk_event$nbar, sig2w_multiplier = risk_event$sig2w_multiplier)
  survG <- temp$S_KM
  ans$survG <- survG

  # data for weighting
  df_weights <- data.frame(at_points, nbar0, nbar1, ybar0, ybar1, surv1, surv0, survP, survG)

  get_score <- wlr_cumulative(df_weights, scheme, scheme_params, return_cumulative = FALSE)

  z.score <- get_score$z.score
  if (!is.finite(z.score)) {
    na_reasons <- c(na_reasons, "weighted log-rank z: its variance is 0 (no event with both arms at risk)")
    z.score <- NA_real_
  }
  ans$z.score <- z.score
  ans$get_score <- get_score
  pval <- 1 - pnorm(z.score)
  ans$zlogrank_text <- paste0("logrank (1-sided) p = ", format_pval(pval, eps = lr.eps, digits = lr.digits))

  # Return lr and sig2_lr directly (duplicating "score" but for historical naming purposes also return as "lr")
  ans$lr <- get_score$score
  ans$sig2_lr <- get_score$sig2.score

  # KM curve checks
  if (check.KM) {
    # survfit on the positive-weight rows (zero-weight rows contribute nothing);
    # control arm will be first stratum
    cd <- data.frame(time = time, delta = delta, z = z, wgt = wgt)[wgt > 0, , drop = FALSE]
    check_fit <- survfit(Surv(time, delta) ~ z, data = cd, weights = wgt, conf.int = conf_level, conf.type = "log")
    # survfit per arm, at that arm's event times (an arm without events has no rows)
    arm_check <- function(a) {
      s <- summary(survfit(Surv(time, delta) ~ 1, data = cd[cd$z == a, , drop = FALSE], weights = wgt))
      data.frame(time = as.numeric(s$time), surv = as.numeric(s$surv), se = as.numeric(s$std.err))
    }
    df0_check <- arm_check(0)
    df1_check <- arm_check(1)

    # survfit's quantile at the same convention (probs = 1 - qprob) and the same confidence level
    qtab <- quantile(check_fit, probs = 1 - qprob)
    if(qprob != 0.50){
      strata_names <- rownames(qtab$quantile)
      quantile_table <- data.frame(
        stratum = rep(strata_names, each = length(qprob)),
        quantile = rep(qprob, times = length(strata_names)),
        time = as.vector(qtab$quantile),
        lower = as.vector(qtab$lower),
        upper = as.vector(qtab$upper)
      )
    ans$quantile_check <- quantile_table
    }
    if(qprob == 0.50){
      quantile_table <- data.frame(
        time = as.vector(qtab$quantile),
        lower = as.vector(qtab$lower),
        upper = as.vector(qtab$upper),
        row.names = rownames(qtab$quantile)
      )
      ans$quantile_check <- quantile_table
    }
    # A real discrepancy is reported whatever verbose is (an error under stop.onerror).
    # Point estimates are always compared; the limits only without case weights, because the
    # weighted limits use the package's own weighted variance.
    q_differs <- function(a, b) any(is.na(a) != is.na(b)) || any(abs(a - b) > 1e-6, na.rm = TRUE)
    cols_cmp <- if (is.null(weight.name)) 1:3 else 1
    # First row is control here
    aa <- c(unlist(quantile_table[1, c("time","lower","upper")]))
    bb <- c(unlist(get_kmq[2, c("quantile","lower","upper")]))
    if (q_differs(aa[cols_cmp], bb[cols_cmp])) {
      msg <- paste0(arms[2]," : ", "Control: discrepancy in quantile calculations")
      if (stop.onerror) stop(msg) else warning(msg)
    }
    aa <- c(unlist(quantile_table[2, c("time","lower","upper")]))
    bb <- c(unlist(get_kmq[1, c("quantile","lower","upper")]))
    if (q_differs(aa[cols_cmp], bb[cols_cmp])) {
      msg <- paste0(arms[1]," : ", "Treatment: discrepancy in quantile calculations")
      if (stop.onerror) stop(msg) else warning(msg)
    }
    # internal check against survfit (distinct from the exported check_km_curve())
    check_km_fit <- function(time, S.KM, se.KM, df_check, group_name = "Group", check.seKM = TRUE) {
        if (any(S.KM < 0 | S.KM > 1)) {
        msg <- paste0(group_name, " : ","KM curve has values outside [0,1].")
        if (stop.onerror) stop(msg) else warning(msg)
      }
      if (any(diff(S.KM) > 0)) {
        msg <- paste0(group_name, " : ", " KM curve is not non-increasing.")
        if (stop.onerror) stop(msg) else warning(msg)
      }

      if(nrow(df_check) > 0 && round(max(abs(S.KM-df_check$surv)),8)) {
        msg <- paste0(group_name," : ", "Discrepancy in KM curve fit.")
        if (stop.onerror) stop(msg) else warning(msg)
      }
      if(check.seKM && nrow(df_check) > 0){
      if(round(max(abs(se.KM-df_check$se)),8)) {
          if (verbose) {
          msg <- paste0(group_name, " : Discrepancy in se(KM) curve fit.")
          # Check for zero or near-zero standard errors to avoid division by zero
          if (any(abs(with(df_check, se)) < .Machine$double.eps)) {
            warning(paste0(group_name, " : Zero or near-zero standard error encountered in df_check$se. Ratio calculation may be invalid."))
            ratio <- NA
          } else {
            ratio <- se.KM / with(df_check, se)
            print(summary(ratio))
          }

          if (stop.onerror) {
            stop(msg)
          } else {
            warning(msg)
          }
        }
      } # checking seKM

              }
      }
    # the package curve is read at survfit's event times
    idc0 <- match(df0_check$time, at_points)
    idc1 <- match(df1_check$time, at_points)
    check_km_fit(df0_check$time, surv0[idc0], sqrt(sig2_surv0[idc0]), df0_check, "control", check.seKM = check.seKM)
    check_km_fit(df1_check$time, surv1[idc1], sqrt(sig2_surv1[idc1]), df1_check, "treat", check.seKM = check.seKM)
     }

  plot_km_comparison <- function(time,surv,se,dfcheck,group_name = "control"){
  if (length(time) == 0 || nrow(dfcheck) == 0) {
    graphics::plot.new(); title(main = paste(group_name, "(no events)")); return(invisible(NULL))
  }
  yymax <- max(c(se, dfcheck$se))
  plot(time,se, type="s", lty=1, col="lightgrey", lwd=4, ylim=c(0,yymax), xlab="time", ylab="SE(KM)")
  with(dfcheck, lines(time, se, type="s", lty=2, lwd=1, col="red"))
   legend("topleft",c("Mine","Survfit"), lty=c(1,2),col=c("lightgrey","red"), lwd=c(4,1),bty="n", cex=0.8)
   title(main=group_name)
  }

  if(check.seKM){
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  par(mfrow=c(1,2))
  plot_km_comparison(time = at_points[idv0.check], surv = surv0[idv0.check], se = sqrt(sig2_surv0[idv0.check]), dfcheck = df0_check, group_name = "control")
  plot_km_comparison(time = at_points[idv1.check], surv = surv1[idv1.check], se = sqrt(sig2_surv1[idv1.check]), dfcheck = df1_check, group_name = "treat")
  }

  ans$at_points <- at_points
  ans$ybar0 <- ybar0
  ans$nbar0 <- nbar0
  ans$surv0 <- surv0
  ans$sig2_surv0 <- sig2_surv0
  ans$ybar1 <- ybar1
  ans$nbar1 <- nbar1
  ans$surv1 <- surv1
  ans$sig2_surv1 <- sig2_surv1
  ans$survP <- survP
  ans$sig2_survP <- sig2_survP

  if (length(stratum) > 1) {
    ybar0_mat <- matrix(NA, nrow = length(at_points), ncol = length(stratum))
    nbar0_mat <- matrix(NA, nrow = length(at_points), ncol = length(stratum))
    ybar1_mat <- matrix(NA, nrow = length(at_points), ncol = length(stratum))
    nbar1_mat <- matrix(NA, nrow = length(at_points), ncol = length(stratum))
    surv0_mat <- matrix(NA, nrow = length(at_points), ncol = length(stratum))
    surv1_mat <- matrix(NA, nrow = length(at_points), ncol = length(stratum))
    sig2_surv0_mat <- matrix(NA, nrow = length(at_points), ncol = length(stratum))
    sig2_surv1_mat <- matrix(NA, nrow = length(at_points), ncol = length(stratum))
    survP_mat <- matrix(NA, nrow = length(at_points), ncol = length(stratum))
    sig2_survP_mat <- matrix(NA, nrow = length(at_points), ncol = length(stratum))
    lr_stratified <- 0.0
    sig2_lr_stratified <- 0.0
    score_stratified <- 0.0
    sig2_score_stratified <- 0.0
    for (ss in seq_along(stratum)) {
      this_stratum <- stratum[ss]
      U0_s <- time[z == 0 & strata == this_stratum]
      D0_s <- delta[z == 0 & strata == this_stratum]
      W0_s <- wgt[z == 0 & strata == this_stratum]

      U1_s <- time[z == 1 & strata == this_stratum]
      D1_s <- delta[z == 1 & strata == this_stratum]
      W1_s <- wgt[z == 1 & strata == this_stratum]

      risk_event <- calculate_risk_event_counts(U1_s, D1_s, W1_s, at_points, draws, seedstart)
      temp <- KM_estimates(ybar = risk_event$ybar, nbar = risk_event$nbar, sig2w_multiplier = risk_event$sig2w_multiplier)
      surv1_mat[, ss] <- temp$S_KM
      sig2_surv1_mat[, ss] <- temp$sig2_KM

      nbar1_s <- risk_event$nbar
      ybar1_s <- risk_event$ybar

      risk_event <- calculate_risk_event_counts(U0_s, D0_s, W0_s, at_points, draws, seedstart)
      temp <- KM_estimates(ybar = risk_event$ybar, nbar = risk_event$nbar, sig2w_multiplier = risk_event$sig2w_multiplier)
      surv0_mat[, ss] <- temp$S_KM
      sig2_surv0_mat[, ss] <- temp$sig2_KM

      nbar0_s <- risk_event$nbar
      ybar0_s <- risk_event$ybar

      # An arm absent from this stratum: its stratum curve does not exist (NA), its counts are 0,
      # so the stratum contributes nothing to the stratified statistics
      for (a in c(0, 1)) {
        if (length(if (a == 0) U0_s else U1_s) == 0) {
          if (a == 0) { surv0_mat[, ss] <- NA; sig2_surv0_mat[, ss] <- NA }
          else { surv1_mat[, ss] <- NA; sig2_surv1_mat[, ss] <- NA }
          na_reasons <- c(na_reasons, sprintf("stratum '%s' has no %s observations: its %s stratum curve is NA and the stratum contributes nothing",
                                              this_stratum, if (a == 0) "control-arm" else "treatment-arm", if (a == 0) "control" else "treatment"))
        }
      }

      U_s <- time[strata == this_stratum]
      D_s <- delta[strata == this_stratum]
      W_s <- wgt[strata == this_stratum]

      risk_event <- calculate_risk_event_counts(U_s, D_s, W_s, at_points, draws, seedstart)
      temp <- KM_estimates(ybar = risk_event$ybar, nbar = risk_event$nbar, sig2w_multiplier = risk_event$sig2w_multiplier)
      survP_mat[, ss] <- temp$S_KM
      sig2_survP_mat[, ss] <- temp$sig2_KM
      S_pool <- temp$S_KM
      # Censoring distribution
      risk_event <- calculate_risk_event_counts(U_s, 1 - D_s, W_s, at_points, draws, seedstart)
      temp <- KM_estimates(ybar = risk_event$ybar, nbar = risk_event$nbar, sig2w_multiplier = risk_event$sig2w_multiplier)
      G_pool <- temp$S_KM

      df_weights <- data.frame(at_points, nbar0 = nbar0_s, nbar1 = nbar1_s, ybar0 = ybar0_s, ybar1 = ybar1_s, survP = S_pool, survG = G_pool)
      get_score <- wlr_cumulative(df_weights, scheme, scheme_params, return_cumulative = FALSE)

      score_stratified <- score_stratified + get_score$score
      sig2_score_stratified <- sig2_score_stratified + get_score$sig2.score

      lr_stratified <- lr_stratified + get_score$score
      sig2_lr_stratified <- sig2_lr_stratified + get_score$sig2.score

      ybar0_mat[, ss] <- ybar0_s
      nbar0_mat[, ss] <- nbar0_s
      ybar1_mat[, ss] <- ybar1_s
      nbar1_mat[, ss] <- nbar1_s
 }

    # For each observation, get the index of its time in at_points
    id_time <- match(time, at_points)  # 'time' is the vector of observed times
    # For each observation, get the index of its stratum in stratum
    id_stratum <- vapply(strata, function(x) which(stratum == x), integer(1))
    # Now extract for each observation
    ans$ybar0_stratum      <- ybar0_mat[cbind(id_time, id_stratum)]
    ans$nbar0_stratum      <- nbar0_mat[cbind(id_time, id_stratum)]
    ans$ybar1_stratum      <- ybar1_mat[cbind(id_time, id_stratum)]
    ans$nbar1_stratum      <- nbar1_mat[cbind(id_time, id_stratum)]
    ans$surv0_stratum      <- surv0_mat[cbind(id_time, id_stratum)]
    ans$sig2_surv0_stratum <- sig2_surv0_mat[cbind(id_time, id_stratum)]
    ans$surv1_stratum      <- surv1_mat[cbind(id_time, id_stratum)]
    ans$sig2_surv1_stratum <- sig2_surv1_mat[cbind(id_time, id_stratum)]
    ans$survP_stratum      <- survP_mat[cbind(id_time, id_stratum)]
    ans$sig2_survP_stratum <- sig2_survP_mat[cbind(id_time, id_stratum)]

    ans$lr_stratified <- lr_stratified
    ans$sig2_lr_stratified <- sig2_lr_stratified
    ans$z.score_stratified <- score_stratified/sqrt(sig2_score_stratified)

    }

  # Compare with survdiff for gamma=0 (survdiff only handles gamma=0)
  # Helper function for robust comparison and messaging
  check_logrank_discrepancy <- function(z, chisq, context = "") {
    tol <- 1e-6
    if (!is.finite(z) || !is.finite(chisq)) return(invisible(NULL))  # NA already reported
    diff <- z^2 - chisq
    if (abs(diff) > tol) {
      warning(sprintf("Discrepancy with log-rank and survdiff%s", context))
      message(sprintf("Discrepancy with survdiff%s: z^2 = %.8f, chisq = %.8f, diff = %.8g",
                      context, z^2, chisq, diff))
    }
  }

  # Compare with survdiff for gamma=0 (survdiff only handles gamma=0)
  if (is.null(weight.name) && is.null(strata.name) && scheme == "fh" && scheme_params$gamma == 0) {
    if (!is.null(get_score$z.score) && !is.null(logrank_results$chisq)) {
      check_logrank_discrepancy(get_score$z.score, logrank_results$chisq)
    } else {
      warning("Required elements missing in get_score or logrank_results (unstratified).")
    }
  }

  if (is.null(weight.name) && !is.null(strata.name) && scheme == "fh" && scheme_params$gamma == 0) {
    if (exists("lr_stratified") && exists("sig2_lr_stratified") && !is.null(logrank_results$chisq)) {
      z_lr <- lr_stratified / sqrt(sig2_lr_stratified)
      check_logrank_discrepancy(z_lr, logrank_results$chisq, " (stratified)")
    } else {
      warning("Required elements missing for stratified log-rank comparison.")
    }
  }
  if (length(na_reasons)) {
    warning("df_counting(): some quantities could not be computed and are NA: ",
            paste(unique(na_reasons), collapse = "; "), call. = FALSE)
  }
  ans
}


