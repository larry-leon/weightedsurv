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
#' @param weight_variance How the weights are treated in standard errors, intervals, bands and tests.
#'   \code{"known"} (default): as fixed, known numbers. \code{"estimated"}: as estimated by the logistic
#'   propensity model given in \code{weight_model}, whose estimation is carried into the inference. Estimates
#'   never change. See the section "Estimated weights".
#' @param weight_model With \code{weight_variance = "estimated"}, the object returned by
#'   \code{\link{ps_weights}} for the rows of \code{df}, in the same order; its weights must be the column
#'   \code{weight.name}. Ignored with \code{"known"}.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{cox_results}{List with the results of \code{survival::coxph()}: \code{cox_fit}, \code{hr},
#'     \code{hr_ci}, \code{pval}, \code{score}, \code{cox_text} and \code{se}, the standard error of the
#'     log hazard ratio behind \code{hr_ci} and \code{pval} (the robust one of \code{coxph(robust = TRUE)};
#'     \code{NA} when the model could not be fitted). \code{score} is the score test of the model: with no
#'     weights, or with all positive weights equal, \code{summary(cox_fit)$sctest}; with subject weights
#'     that are not all equal, the robust score test \code{summary(cox_fit)$robscore} (the model-based
#'     one reads the weights as counts).}
#'   \item{logrank_results}{List with \code{chisq}, \code{pval}, \code{logrank_text} and \code{method}.
#'     With no weights, or with all positive weights equal, the test of \code{survival::survdiff()}
#'     (\code{method = "survdiff"}; \code{rho} from the scheme when it is \code{"fh"}). With subject weights
#'     that are not all equal \code{survdiff()} has no counterpart, and the result is the robust score
#'     test at 0 of the weighted Cox model (\code{coxph(weights = , init = 0, iter.max = 0, robust = TRUE,
#'     ties = "breslow")} on the positive-weight rows, strata included; \code{method =
#'     "coxph_robust_score"}): the \code{rho = 0} test whatever the scheme. Its \code{chisq} equals
#'     \code{z.score^2} (\code{z.score_stratified^2}) for Fleming-Harrington (0, 0) weights, so
#'     \code{\link{check_results}} agrees on weighted data.}
#'   \item{z.score}{Standardized weighted log-rank test statistic (with subject weights that are not
#'     all equal, standardized by the score-residual variance; see the section on subject weights)}
#'   \item{at_points}{Vector of all time points used in calculations}
#'   \item{surv0, surv1}{Kaplan-Meier survival estimates for control and treatment groups}
#'   \item{sig2_surv0, sig2_surv1}{Variance estimates for survival curves}
#'   \item{survP}{Pooled survival estimates}
#'   \item{survG}{Censoring distribution estimates}
#'   \item{quantile_results}{Data frame with median survival and confidence intervals by group}
#'   \item{lr, sig2_lr}{Weighted log-rank statistic and its variance (\code{get_score} holds the same as
#'     returned by \code{\link{wlr_cumulative}}; \code{zlogrank_text} is the one-sided p-value of \code{z.score})}
#'   \item{riskpoints0, riskpoints1}{Risk set counts at specified time points}
#'   \item{z.score_stratified}{Stratified z-score (if stratified analysis), with \code{lr_stratified} and
#'     \code{sig2_lr_stratified}}
#'   \item{ess}{Data frame with one row per arm, in the order of \code{quantile_results}
#'     (\code{arms[1]}, the treated arm, then \code{arms[2]}): \code{group}; \code{n}, the rows of the arm
#'     with positive weight; \code{ess}, the effective sample size
#'     \eqn{(\sum w)^2 / \sum w^2}{(sum w)^2 / sum(w^2)} over those rows; and \code{ratio = ess / n}.
#'     Without weights, and with constant weights, \code{ess = n}. A ratio near 1 means the weights are
#'     even; a low ratio means that a few subjects carry most of the weight, and standard errors,
#'     intervals and tests can then be too optimistic, whether the weights are treated as known or
#'     as estimated. No warning or message depends on it.}
#'   \item{weight_variance}{\code{"known"} or \code{"estimated"}: how the weights were treated.}
#'   \item{ew_parts}{Only with \code{weight_variance = "estimated"}: list with \code{IF} and \code{dw}
#'     of \code{weight_model}, rows in the order of \code{time}; read by \code{\link{cox_rhogamma}},
#'     \code{\link{wlr_cumulative}} and \code{\link{wlr_dhat_estimates}}.}
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
#' fitted on the positive-weight rows. Zero-weight rows do not affect which variance formula runs
#' either: whether the weights are "all equal" is decided on the positive weights, so a 0/1 weight
#' column gives the results of the positive-weight rows analysed without weights. A stratum in which one arm is empty contributes nothing to the
#' stratified statistics (that arm's stratum curve is \code{NA}). On degenerate data (an arm without events, no events at all, one subject
#' per arm) a result is still returned; quantities that cannot be computed are \code{NA}, with one
#' warning that says why.
#' For stratified analyses, stratum-specific estimates are computed and combined appropriately.
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
#' @section Variance of the log-rank and Cox score with subject weights:
#' With no subject weights, or with all subject weights equal, the weighted log-rank test and the
#' weighted Cox estimate are computed as before, and their results are unchanged. When the subject
#' weights are not all equal, the variance of the weighted score is the sum of squared per-subject
#' score contributions,
#' \deqn{\widehat{Var}\{U(\beta)\} = \sum_i \psi_i(\beta)^2,}{Var(U(b)) = sum_i psi_i(b)^2,}
#' \deqn{\psi_i(\beta) = a_i \Big( D_i K(U_i) \{Z_i - \bar z(U_i, \beta)\} - e^{\beta Z_i} \sum_{s \le U_i} K(s) \{Z_i - \bar z(s, \beta)\} d\hat\Lambda(s, \beta) \Big),}{psi_i(b) = a_i ( D_i K(U_i) (Z_i - zbar(U_i, b)) - exp(b Z_i) sum_{s <= U_i} K(s) (Z_i - zbar(s, b)) dL(s, b) ),}
#' where \eqn{U_i}, \eqn{D_i}, \eqn{Z_i} and \eqn{a_i} are subject \eqn{i}'s time, event indicator,
#' arm (0/1) and subject weight, \eqn{K(s)} is the time weight of the scheme,
#' \eqn{S_0(s, \beta) = \sum_j a_j 1(U_j \ge s) e^{\beta Z_j}}{S0(s, b) = sum_j a_j 1(U_j >= s) exp(b Z_j)},
#' \eqn{\bar z(s, \beta) = \sum_j a_j 1(U_j \ge s) Z_j e^{\beta Z_j} / S_0(s, \beta)}{zbar(s, b) = sum_j a_j 1(U_j >= s) Z_j exp(b Z_j) / S0(s, b)} and
#' \eqn{d\hat\Lambda(s, \beta) = \sum_j a_j D_j 1(U_j = s) / S_0(s, \beta)}{dL(s, b) = sum_j a_j D_j 1(U_j = s) / S0(s, b)}.
#' Sums over \eqn{s} run over distinct times (events at a tied time are added up), and a term whose
#' \eqn{S_0} is 0 is 0. The score itself, \eqn{-\sum_i \psi_i(\beta)}{-sum_i psi_i(b)}, is the one
#' computed before and does not change, and neither does any point estimate. The log-rank z-score
#' is the score at \eqn{\beta = 0} over \eqn{\sqrt{\sum_i \psi_i(0)^2}}{sqrt(sum_i psi_i(0)^2)}; the
#' variance of the log hazard ratio estimate is
#' \eqn{\sum_i \psi_i(\hat\beta)^2 / i(\hat\beta)^2}{sum_i psi_i(bhat)^2 / i(bhat)^2} with
#' \eqn{i(\hat\beta)}{i(bhat)} the information (unchanged); in a stratified analysis the
#' contributions are computed within each stratum and the score and its variance are sums over
#' strata. For Fleming-Harrington (0, 0) weights this is the robust variance of
#' \code{survival::coxph(weights = w, robust = TRUE, ties = "breslow")}: the same robust score test
#' and the same robust standard error. Results no longer depend on how the weights are scaled
#' (weights \code{w}, \code{10 * w} and \code{w / sum(w)} give the same z-score and standard error).
#' By default the weights are treated as fixed and known. For weights estimated by a logistic
#' propensity model, \code{weight_variance = "estimated"} carries their estimation into the test and
#' the intervals (section "Estimated weights"). When the weights are extreme (low effective sample
#' size) the intervals can be too short under either option.
#'
#' @section Estimated weights:
#' Propensity-score weights are estimated from the data. \code{weight_variance} says how standard
#' errors, intervals, bands and tests treat them.
#' \describe{
#'   \item{\code{"known"} (default)}{The weights are fixed, known numbers: the formulas of the two
#'     sections above. \code{weight_model} is ignored.}
#'   \item{\code{"estimated"}}{The weights were estimated by a logistic regression for the treatment,
#'     given as \code{weight_model} (the object returned by \code{\link{ps_weights}}), and that
#'     estimation is carried into the inference.}
#' }
#' Estimates never change: survival curves, their difference, RMST, the log-rank score, hazard
#' ratios, counts and risk tables are identical under the two options. Only standard errors,
#' intervals, bands and tests differ.
#'
#' Notation. There are \eqn{n} subjects, the rows analysed. Subject \eqn{i} has arm \eqn{A_i},
#' row \eqn{x_i} of the model matrix of the propensity model (\eqn{p} columns), fitted probability
#' \eqn{e_i} and weight \eqn{w_i} (written \eqn{a_i} in the score formulas above). The influence of
#' subject \eqn{i} on the fitted coefficients and the derivative of its weight in the coefficients are
#' \deqn{IF_i = \Big\{\sum_j x_j x_j' e_j (1 - e_j)\Big\}^{-1} x_i (A_i - e_i), \qquad dw_i = -\frac{1 - e_i}{e_i} x_i \ (A_i = 1), \quad dw_i = \frac{e_i}{1 - e_i} x_i \ (A_i = 0)}{IF_i = { sum_j x_j x_j' e_j (1 - e_j) }^(-1) x_i (A_i - e_i);  dw_i = -((1 - e_i) / e_i) x_i (treated),  dw_i = (e_i / (1 - e_i)) x_i (control)}
#' (for stabilised weights \eqn{w_i} and \eqn{dw_i} are multiplied by the arm proportion).
#'
#' Kaplan-Meier estimate of a group \eqn{g} (an arm; an arm within a stratum; all subjects). With
#' \eqn{\psi_i(t)}{psi_i(t)} the influence value of the section "Variance of the Kaplan-Meier
#' estimate", computed from the group's own weighted counts,
#' \deqn{D_g(t) = \sum_{k \in g} \psi_k(t)\, dw_k, \qquad \Phi_i(t) = 1(i \in g)\, w_i \psi_i(t) + IF_i' D_g(t) \quad (i = 1, \ldots, n),}{D_g(t) = sum_{k in g} psi_k(t) dw_k,   Phi_i(t) = 1(i in g) w_i psi_i(t) + IF_i' D_g(t)   (i = 1, ..., n),}
#' \deqn{\widehat{Var}[\hat S_g(t)] = \hat S_g(t)^2 \sum_{i=1}^n \Phi_i(t)^2.}{Var[S_g(t)] = S_g(t)^2 * sum_{i=1}^n Phi_i(t)^2.}
#' Every subject contributes, also those outside the group, through the model. The variance of the
#' difference of the two arms is
#' \deqn{\widehat{Var}[\hat S_1(t) - \hat S_0(t)] = \sum_{i=1}^n \{\hat S_1(t) \Phi^1_i(t) - \hat S_0(t) \Phi^0_i(t)\}^2,}{Var[S_1(t) - S_0(t)] = sum_{i=1}^n ( S_1(t) Phi1_i(t) - S_0(t) Phi0_i(t) )^2,}
#' which is \strong{not} the sum of the two arm variances: the arms share the model's coefficients.
#' A resampled curve of arm \eqn{a} is
#' \eqn{-\hat S_a(t) \sum_{i=1}^n G_i \Phi^a_i(t)}{-S_a(t) sum_{i=1}^n G_i Phi^a_i(t)} with the
#' standard-normal multipliers \eqn{G_i} that \code{\link{KM_diff}} generates in any case; the
#' resampled difference curves, the simultaneous band and the RMST variance and band of
#' \code{\link{cumulative_rmst_bands}} follow from them unchanged.
#'
#' Score of the log-rank test and of the weighted Cox estimate. With
#' \eqn{\psi_i(\beta)}{psi_i(b)} the contribution of the section "Variance of the log-rank and Cox
#' score with subject weights", computed within each stratum,
#' \deqn{B(\beta) = \sum_{j=1}^n \frac{\psi_j(\beta)}{a_j}\, dw_j, \qquad \psi^{ew}_i(\beta) = \psi_i(\beta) + IF_i' B(\beta), \qquad \widehat{Var}\{U(\beta)\} = \sum_{i=1}^n \psi^{ew}_i(\beta)^2.}{B(b) = sum_{j=1}^n (psi_j(b) / a_j) dw_j,   psi_ew_i(b) = psi_i(b) + IF_i' B(b),   Var(U(b)) = sum_{i=1}^n psi_ew_i(b)^2.}
#' The log-rank z-score, the variance of the log hazard ratio estimate (the sum at
#' \eqn{\hat\beta}{bhat} over the squared information) and a resampled score
#' (\eqn{-\sum_i G_i \psi^{ew}_i(\beta)}{-sum_i G_i psi_ew_i(b)}) are formed as before with
#' \eqn{\psi^{ew}}{psi_ew} in place of \eqn{\psi}{psi}. In a stratified analysis the variance is this
#' single sum over all subjects, not a sum of per-stratum variances. The cumulative variance of
#' \code{\link{wlr_cumulative}} is \eqn{\sum_i \{c_i(t) + IF_i' B(t)\}^2}{sum_i (c_i(t) + IF_i' B(t))^2}
#' with \eqn{c_i(t)} the contribution of subject \eqn{i} to the score up to \eqn{t} and
#' \eqn{B(t) = \sum_j \{c_j(t) / a_j\} dw_j}{B(t) = sum_j (c_j(t) / a_j) dw_j}; in
#' \code{\link{wlr_dhat_estimates}} the contribution to the survival difference is
#' \eqn{d_i + IF_i' (\hat S_0 D_0 - \hat S_1 D_1)}{d_i + IF_i' (S_0 D_0 - S_1 D_1)} at the grid point
#' read for \code{tzero}, and the covariance uses both adjusted contributions. The outputs of
#' \code{df_counting()} that come from the \code{survival} package (\code{logrank_results},
#' \code{cox_results$score}, \code{$se}, \code{$hr_ci}, \code{$pval}) carry the same adjustment,
#' applied to \code{coxph()}'s weighted score and dfbeta residuals \eqn{r}:
#' \eqn{r_i + IF_i' \sum_j (r_j / w_j) dw_j}{r_i + IF_i' sum_j (r_j / w_j) dw_j}.
#'
#' How to use it. Fit the propensity model on exactly the rows to be analysed with
#' \code{\link{ps_weights}}, put its weights in a column of \code{df}, and pass
#' \code{weight_variance = "estimated"} and \code{weight_model} to \code{df_counting()},
#' \code{\link{KM_diff}}, \code{\link{cumulative_rmst_bands}} or \code{\link{plotKM.band_subgroups}}.
#' \code{\link{cox_rhogamma}}, \code{\link{wlr_cumulative}} and \code{\link{wlr_dhat_estimates}} take the
#' option from the \code{df_counting()} result they are given (its elements \code{weight_variance} and
#' \code{ew_parts}). The rows of \code{df}, the treatment column and the weight column must be those of
#' \code{weight_model}; anything else is an error.
#'
#' Scope and limits of this version.
#' \itemize{
#'   \item A logistic \code{glm} for a binary treatment, fitted on exactly the rows analysed. For
#'     weights from any other model use \code{"known"}.
#'   \item A subgroup analysis that reuses weights from a model fitted on the full sample must use
#'     \code{"known"}: the model's rows are then not the rows analysed. (The subgroup curves of
#'     \code{plotKM.band_subgroups()} are point estimates and are not affected.)
#'   \item \code{draws > 0} (the resampling-based pointwise variance) is not available with
#'     \code{"estimated"}; \code{draws.band} and the resampling of \code{cox_rhogamma()} are.
#'   \item The time weights of the scheme (Fleming-Harrington and the others) depend on the weighted
#'     pooled curve; that dependence on the estimated weights is ignored.
#'   \item With stabilised weights the arm proportion is treated as a fixed number. The Kaplan-Meier
#'     quantities do not depend on that: they are unchanged when an arm's weights are multiplied by a
#'     constant.
#'   \item With \code{"estimated"} the influence-based formulas always run; whether the weights are all
#'     equal is not tested.
#'   \item Carrying the estimation through does not repair extreme weights: when the effective sample
#'     size (\code{df_counting()$ess}) is a small fraction of an arm, standard errors, intervals and
#'     tests can be too optimistic under either option.
#' }
#' The device of carrying the estimation of a nuisance quantity into the inference by perturbing its
#' estimate is that of Leon, Cai and Wei (2009); it is used here in linearised form.
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
#' Leon, L. F., Cai, T. and Wei, L. J. (2009). Robust inferences for covariate effects on survival
#' time with censored linear regression models. Statistics in Biosciences, 1, 50-64.
#' \doi{10.1007/s12561-009-9002-5}
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
                        stop.onerror=FALSE,censoring_allmarks=TRUE,
                        weight_variance = c("known", "estimated"), weight_model = NULL) {

  weight_variance <- match.arg(weight_variance)

  # Columns present (strata too); no NA; times finite and >= 0; event 0/1, logical or 1/2;
  # weights >= 0; both arms non-empty. Errors name the column and the problem.
  vin <- ws_validate_survdata(df, tte.name, event.name, treat.name, weight.name, strata.name)
  # 1/2 event coding: use the 0/1 recoding everywhere, so that subsets cannot be misread
  if (any(vin$event != df[[event.name]])) df[[event.name]] <- vin$event
  # "estimated": the pieces of the propensity model (NULL for "known"), rows in the order of df
  ew <- ws_ew_setup(weight_variance, weight_model, df, treat.name, weight.name, draws, "df_counting")
  ew_on <- !is.null(ew)
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
                 cox_text = "HR = NA", se = NA_real_)
  # Subject weights whose positive values are not all equal (Rb): the survival-based outputs are the
  # robust ones of the weighted Cox model (score test of cox_fit; robust score test at 0 as log-rank).
  rb_unequal <- (!is.null(weight.name) && ws_weights_unequal(vin$weight)) || ew_on
  # Estimated weights: a vector r of weighted residuals of a coxph fit (rows of df) gets the term through
  # the propensity model, r + IF (sum_j (r_j / w_j) dw_j)
  ew_df <- ew                                   # rows in the order of df (before sorting by time)
  ew_adj <- function(r) {
    r <- as.numeric(r)
    as.numeric(r + ew_df$IF %*% colSums((r / vin$weight) * ew_df$dw))
  }
  treat.name.strata <- treat.name
  if (!is.null(strata.name)) {
    treat.name.strata <- paste(treat.name,"+",paste("strata(",eval(strata.name),")"))
  }
  cox_formula <- as.formula(paste0("survival::Surv(", tte.name, ",", event.name, ") ~ ", treat.name.strata))
  if (get.cox && any(events_by_arm == 0)) {
    na_reasons <- c(na_reasons, sprintf("Cox model: no events in the %s; the hazard ratio is not identified",
                                        paste(c("control arm", "treatment arm")[events_by_arm == 0], collapse = " and ")))
    cox_results <- cox_na
  } else if (get.cox) {
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
    cox_score <- if (rb_unequal) cox_summary$robscore[1] else cox_summary$sctest[1]
    hr <- exp(cox_fit$coef)
    hr_ci <- exp(confint(cox_fit))
    pval <- cox_summary$coefficients[1, "Pr(>|z|)"]
    cox_se <- sqrt(cox_fit$var[1, 1])
    if (ew_on) {
      # Estimated weights: the same adjustment through survival's residuals. Standard error from the
      # dfbeta residuals of cox_fit; score test from the score residuals of the fit at 0 (cox_fit's ties).
      ew_cox <- tryCatch({
        se <- sqrt(sum(ew_adj(stats::residuals(cox_fit, type = "dfbeta", weighted = TRUE))^2))
        fit0 <- survival::coxph(cox_formula, data = df_fit, weights = df_fit[[weight.name]], init = 0,
                                control = survival::coxph.control(iter.max = 0), robust = TRUE)
        r0 <- as.numeric(stats::residuals(fit0, type = "score", weighted = TRUE))
        list(se = se, score = sum(r0)^2 / sum(ew_adj(r0)^2))
      }, error = function(e) {
        na_reasons <<- c(na_reasons, paste("Cox model with estimated weights failed:", conditionMessage(e)))
        list(se = NA_real_, score = NA_real_)
      })
      cox_se <- ew_cox$se
      cox_score[1] <- ew_cox$score
      bh <- unname(cox_fit$coef[1])
      hr_ci[] <- exp(bh + c(-1, 1) * stats::qnorm(0.975) * cox_se)
      pval <- 2 * stats::pnorm(-abs(bh / cox_se))
    }
    cox_text <- paste0("HR = ", round(hr, cox.digits),
                        " (", round(hr_ci[1], cox.digits), ", ", round(hr_ci[2], cox.digits), ")")

     cox_results <- list(
      cox_fit = cox_fit,
      hr = hr,
      hr_ci = hr_ci,
      pval = pval,
      score = cox_score,
      cox_text = cox_text,
      se = cox_se
    )
  }
  ans$cox_results <- cox_results
  logrank_results <- NULL
  surv_obj <- survival::Surv(df_fit[[tte.name]], df_fit[[event.name]])
  group <- df_fit[[treat.name]]
 # survdiff does not support case weights

  lr_rho <- 0
  if(scheme == "fh") lr_rho <- scheme_params$rho

  logrank_method <- if (rb_unequal) "coxph_robust_score" else "survdiff"
  if (sum(events_by_arm) == 0) {
    na_reasons <- c(na_reasons, "log-rank test (survdiff): no events")
    logrank_fit <- list(chisq = NA_real_)
  } else if (rb_unequal) {
    # survdiff() takes no case weights: the robust score test at 0 of the weighted Cox model
    # (positive-weight rows, the formula of cox_fit, Breslow ties) -- the rho = 0 test whatever the scheme
    logrank_fit <- tryCatch({
      lr_fit <- survival::coxph(cox_formula, data = df_fit, weights = df_fit[[weight.name]], init = 0,
                                control = survival::coxph.control(iter.max = 0), robust = TRUE, ties = "breslow")
      if (ew_on) {
        r <- as.numeric(stats::residuals(lr_fit, type = "score", weighted = TRUE))
        list(chisq = sum(r)^2 / sum(ew_adj(r)^2))
      } else list(chisq = as.numeric(lr_fit$rscore))
    }, error = function(e) {
      na_reasons <<- c(na_reasons, paste("log-rank test (robust score test of coxph) failed:", conditionMessage(e)))
      list(chisq = NA_real_)
    })
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
    logrank_text = logrank_text,
    method = logrank_method
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
    ew <- ws_ew_reorder(ew, ord)
  }
  # Estimated weights: variance of the Kaplan-Meier curve of a group (rows `rows` of the sorted data)
  ew_sig2 <- function(rows, risk_event, S) {
    ws_km_ew_sig2(time[rows], delta[rows], wgt[rows], rows, at_points, risk_event$ybar, risk_event$nbar, S, ew)
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
  sig2_surv1 <- if (ew_on) ew_sig2(which(z == 1), risk_event, temp$S_KM) else temp$sig2_KM

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
  sig2_surv0 <- if (ew_on) ew_sig2(which(z == 0), risk_event, temp$S_KM) else temp$sig2_KM

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
  sig2_survP <- if (ew_on) ew_sig2(seq_along(time), risk_event, temp$S_KM) else temp$sig2_KM

  # Pooled Censoring estimates
  risk_event <- calculate_risk_event_counts(time, 1- delta, wgt, at_points)
  temp <- KM_estimates(ybar = risk_event$ybar, nbar = risk_event$nbar, sig2w_multiplier = risk_event$sig2w_multiplier)
  survG <- temp$S_KM
  ans$survG <- survG

  # data for weighting
  df_weights <- data.frame(at_points, nbar0, nbar1, ybar0, ybar1, surv1, surv0, survP, survG)
  # Subject weights not all equal: hand the subject-level data to wlr_cumulative(), which then uses
  # the score-residual variance (in every stratum too). Otherwise nothing is attached.
  # With weight_variance = "estimated" the influence-based formulas always run, with the term through the
  # propensity model (the pieces of the model go along with the subject-level data).
  lr_unequal <- ws_weights_unequal(wgt) || ew_on
  if (lr_unequal) attr(df_weights, "ws_subject") <- list(time = time, delta = delta, z = z, w_hat = wgt, unequal = TRUE,
                                                         ew = if (ew_on) list(IF = ew$IF, dw = ew$dw) else NULL)

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
    # with estimated weights the standard error differs from survfit's by design: no discrepancy warning for it
    check_km_fit(df0_check$time, surv0[idc0], sqrt(sig2_surv0[idc0]), df0_check, "control", check.seKM = check.seKM && !ew_on)
    check_km_fit(df1_check$time, surv1[idc1], sqrt(sig2_surv1[idc1]), df1_check, "treat", check.seKM = check.seKM && !ew_on)
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
    psi_strata <- numeric(length(time))      # estimated weights: the score contributions, computed within stratum
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
      sig2_surv1_mat[, ss] <- if (ew_on) ew_sig2(which(z == 1 & strata == this_stratum), risk_event, temp$S_KM) else temp$sig2_KM

      nbar1_s <- risk_event$nbar
      ybar1_s <- risk_event$ybar

      risk_event <- calculate_risk_event_counts(U0_s, D0_s, W0_s, at_points, draws, seedstart)
      temp <- KM_estimates(ybar = risk_event$ybar, nbar = risk_event$nbar, sig2w_multiplier = risk_event$sig2w_multiplier)
      surv0_mat[, ss] <- temp$S_KM
      sig2_surv0_mat[, ss] <- if (ew_on) ew_sig2(which(z == 0 & strata == this_stratum), risk_event, temp$S_KM) else temp$sig2_KM

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
      sig2_survP_mat[, ss] <- if (ew_on) ew_sig2(which(strata == this_stratum), risk_event, temp$S_KM) else temp$sig2_KM
      S_pool <- temp$S_KM
      # Censoring distribution
      risk_event <- calculate_risk_event_counts(U_s, 1 - D_s, W_s, at_points, draws, seedstart)
      temp <- KM_estimates(ybar = risk_event$ybar, nbar = risk_event$nbar, sig2w_multiplier = risk_event$sig2w_multiplier)
      G_pool <- temp$S_KM

      df_weights <- data.frame(at_points, nbar0 = nbar0_s, nbar1 = nbar1_s, ybar0 = ybar0_s, ybar1 = ybar1_s, survP = S_pool, survG = G_pool)
      if (lr_unequal) {
        in_s <- strata == this_stratum
        attr(df_weights, "ws_subject") <- list(time = time[in_s], delta = delta[in_s], z = z[in_s], w_hat = wgt[in_s], unequal = TRUE)
      }
      get_score <- wlr_cumulative(df_weights, scheme, scheme_params, return_cumulative = FALSE)
      if (ew_on) psi_strata[in_s] <- ws_lr_psi(df_weights, scheme, scheme_params)

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

    # Estimated weights: one sum over all subjects of (psi_i + IF_i' B)^2, B summed over all strata
    # (not a sum of per-stratum variances: the model's coefficients are shared)
    if (ew_on) sig2_lr_stratified <- sig2_score_stratified <- sum(ws_score_ew(psi_strata, wgt, ew$IF, ew$dw)^2)
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
  # Effective sample size by arm, (sum w)^2 / sum(w^2) over the rows with positive weight (Ra);
  # rows in the order of quantile_results. No warning or message depends on it.
  ess_arm <- function(a) {
    w <- wgt[z == a & wgt > 0]
    c(n = length(w), ess = sum(w)^2 / sum(w^2))
  }
  e1 <- ess_arm(1); e0 <- ess_arm(0)
  ans$ess <- data.frame(group = arms[1:2], n = unname(c(e1["n"], e0["n"])), ess = unname(c(e1["ess"], e0["ess"])),
                        stringsAsFactors = FALSE)
  ans$ess$ratio <- ans$ess$ess / ans$ess$n
  # How the weights were treated; with "estimated", the pieces of the propensity model in the order of `time`
  ans$weight_variance <- weight_variance
  if (ew_on) ans$ew_parts <- list(IF = ew$IF, dw = ew$dw)

  if (length(na_reasons)) {
    warning("df_counting(): some quantities could not be computed and are NA: ",
            paste(unique(na_reasons), collapse = "; "), call. = FALSE)
  }
  ans
}


