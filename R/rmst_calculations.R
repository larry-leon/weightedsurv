#' Cumulative RMST bands for survival curves
#'
#' Calculates cumulative Restricted Mean Survival Time (RMST) and confidence bands for survival curves using resampling.
#' Optionally plots the cumulative RMST curve, pointwise confidence intervals, and simultaneous confidence bands.
#'
#' @param df Data frame containing survival data.
#' @param fit Survival fit object (output from KM_diff).
#' @param tte.name Name of time-to-event variable in \code{df}.
#' @param event.name Name of event indicator variable in \code{df}.
#' @param treat.name Name of treatment group variable in \code{df}.
#' @param weight.name Optional name of weights variable in \code{df}.
#' @param draws_sb Number of resampling draws for simultaneous bands (default: 1000).
#' @param xlab Label for x-axis (default: "months").
#' @param ylim_pad Padding for y-axis limits (default: 0.5).
#' @param rmst_max_legend Position for RMST legend (default: "left").
#' @param rmst_max_cex Text size for RMST legend (default: 0.7).
#' @param plot Logical; if TRUE, plot the results. Default is TRUE.
#' @param mean_rule Character; how the survival difference is integrated over the time grid of
#'   \code{fit} (\code{fit$at_points}). \code{"trapezoid"} (default) averages the difference at
#'   consecutive grid points: \eqn{R(t_j) = \sum_{k \le j} \{d(t_{k-1}) + d(t_k)\}/2 \, (t_k - t_{k-1})}.
#'   \code{"rectangle"} is the left-rectangle sum
#'   \eqn{R(t_j) = \sum_{k \le j} d(t_{k-1}) (t_k - t_{k-1})}, the exact area under the step
#'   function. In both, \eqn{R = 0} at the first grid point. \code{"rectangle"} equals the
#'   difference of \code{survfit()}'s restricted means (\code{summary(fit, rmean = tau)}) when the
#'   time grid holds every event time, which it does for the default \code{KM_diff()} grid; it does
#'   not when the fit was made with \code{modify_tau = TRUE}. The rule applies to the estimate and
#'   to both sets of resampled curves; nothing else changes.
#'
#' @details The pointwise variance of the cumulative RMST difference and the simultaneous band come
#' from the resampled survival-difference curves of \code{fit} (\code{fit$dhat_star}) and of a
#' second \code{\link{KM_diff}} call with \code{seedstart = 99999}. With weights that are not all
#' equal, those resampled curves are built from the per-subject contributions of the robust
#' variance (see \code{\link{resampling_survival}}), so the RMST variance and bands are the robust
#' ones; with no weights or all weights equal they are unchanged. The weights are treated as fixed
#' and known: with estimated propensity-score weights the intervals are conservative, and with
#' extreme weights (low effective sample size) they are too narrow.
#'
#' @return A list with elements:
#'   \item{at_points}{Time points used for RMST calculation}
#'   \item{rmst_time}{Cumulative RMST estimates}
#'   \item{sig2_rmst_time}{Variance of RMST estimates}
#'   \item{rmst_time_lower}{Pointwise lower confidence interval}
#'   \item{rmst_time_upper}{Pointwise upper confidence interval}
#'   \item{rmst_maxtau_ci}{RMST and CI at maximum time}
#'   \item{rmst_text}{Text summary for legend}
#'   \item{c_alpha_band}{Critical value for simultaneous band}
#'   \item{rmst_time_sb_lower}{Simultaneous band lower bound}
#'   \item{rmst_time_sb_upper}{Simultaneous band upper bound}
#'   \item{mean_rule}{The integration rule used}
#'
#' @importFrom stats quantile rnorm
#' @importFrom utils head tail
#' @importFrom graphics polygon lines legend abline plot
#' @export

cumulative_rmst_bands <- function(df, fit, tte.name, event.name, treat.name, weight.name = NULL, draws_sb = 1000, xlab="months", ylim_pad = 0.5,
                                  rmst_max_legend = "left", rmst_max_cex = 0.7, plot = TRUE,
                                  mean_rule = c("trapezoid", "rectangle")) {
  mean_rule <- match.arg(mean_rule)
  # Input validation
  if (!is.data.frame(df)) stop("df must be a data frame.")
  required_names <- c(tte.name, event.name, treat.name, weight.name)
  missing <- setdiff(required_names, names(df))
  if (length(missing) > 0) stop(sprintf("Missing columns in df: %s", paste(missing, collapse = ", ")))
  if (!is.list(fit)) stop("fit must be a list (output from KM_diff).")
  if (is.null(fit[["dhat_star"]])) {
    stop("cumulative_rmst_bands() needs a KM_diff() fit made with draws.band > 0 (fit$dhat_star is missing).")
  }

ans <- list()
at_points <- fit$at_points
dhat <- fit$dhat
dhat_star_mat <- fit$dhat_star
risk.points <- fit$risk.points

ans$at_points <- at_points

dt <- diff(at_points)

rmst_time <- ws_cumulative_mean(dhat, dt, mean_rule)

ans$rmst_time <- rmst_time

# For multiple draws (matrix: rows=time, cols=draws)
cum_rmst_stars <- apply(dhat_star_mat, 2, ws_cumulative_mean, dt = dt, rule = mean_rule)
sig2_rmst_time <- apply(cum_rmst_stars, 1, var, na.rm=TRUE)

ans$sig2_rmst_time <- sig2_rmst_time

# Pointwise CIs
rmst_time_lower <- rmst_time - 1.96*sqrt(sig2_rmst_time)
rmst_time_upper <- rmst_time + 1.96*sqrt(sig2_rmst_time)

ans$rmst_time_lower <- rmst_time_lower
ans$rmst_time_upper <- rmst_time_upper

# Extract RMST at max time
at_maxtau <- length(rmst_time)
rmst_maxtau <- rmst_time[at_maxtau]
rmst_maxtau_lower <- rmst_time_lower[at_maxtau]
rmst_maxtau_upper <- rmst_time_upper[at_maxtau]
rmst_maxtau_ci <- c(rmst_maxtau, rmst_maxtau_lower, rmst_maxtau_upper)

ans$rmst_maxtau_ci <- rmst_maxtau_ci

rmst_text <- paste0("RMST(tau*) = ", round(rmst_maxtau, 1),
                   " (", round(rmst_maxtau_lower, 1), ", ", round(rmst_maxtau_upper, 1), ")")

ans$rmst_text <- rmst_text

# Centered resamples
fit_draws <- KM_diff(
  df = df, tte.name = tte.name, event.name = event.name , treat.name = treat.name, weight.name = weight.name,
  at_points = at_points, alpha = 0.05, risk.points = risk.points,
  modify_tau = FALSE,
  draws.band = draws_sb, seedstart = 99999, show_resamples = FALSE
)

dhat_star_mat2 <- fit_draws$dhat_star

cum_rmst_stars2 <- apply(dhat_star_mat2, 2, ws_cumulative_mean, dt = dt, rule = mean_rule)

# Standardized resamples
rmst_time_star <- cum_rmst_stars2 / sqrt(sig2_rmst_time)

# simultaneous band
sups <- apply(abs(rmst_time_star), 2, max, na.rm = TRUE)
c_alpha_band <- quantile(sups,c(0.95))

ans$c_alpha_band <- c_alpha_band

# Simultaneous bands
rmst_time_sb_lower <- rmst_time - c_alpha_band * sqrt(sig2_rmst_time)
rmst_time_sb_upper <- rmst_time + c_alpha_band * sqrt(sig2_rmst_time)

ans$rmst_time_sb_lower <- rmst_time_sb_lower
ans$rmst_time_sb_upper <- rmst_time_sb_upper
ans$mean_rule <- mean_rule

if(plot){
x <- at_points
mean.value <- rmst_time
l0_pw <- rmst_time_lower
u0_pw <- rmst_time_upper
l0_sb <- rmst_time_sb_lower
u0_sb <- rmst_time_sb_upper

time.zero.label <- 0.0

ymin <- min(l0_sb)
ymax <- max(u0_sb) + ylim_pad

plot(
  x[order(x)], mean.value[order(x)], type = "n", xlab = xlab, lty = 1,
  ylab = "Cumulative RMST", ylim = c(ymin,ymax), cex.lab = 1
)
polygon(
  c(x[order(x)], rev(x[order(x)])),
  c(l0_sb[order(x)], rev(u0_sb[order(x)])),
  col = "lightgrey", border = FALSE
)
lines(x[order(x)], l0_pw[order(x)], lty=2, type="s")
lines(x[order(x)], u0_pw[order(x)], lty=2, type="s")

lines(x[order(x)], mean.value[order(x)], lty = 1, lwd = 1, type = "s")
abline(h = time.zero.label, lty = 1, col = "blue", lwd = 0.5)

legend(rmst_max_legend, legend = rmst_text, cex = rmst_max_cex, bty = "n")
}

return(invisible(ans))
}

