#' Propensity-score weights with the pieces needed to carry their estimation into the inference
#'
#' Fits (or takes) a logistic regression for a binary treatment and returns inverse-probability
#' weights together with the influence of the fitted coefficients and the derivative of each weight
#' in those coefficients. Passed as `weight_model` to [df_counting()], [KM_diff()],
#' [cumulative_rmst_bands()] or [plotKM.band_subgroups()] with `weight_variance = "estimated"`,
#' the object lets standard errors, intervals, bands and tests account for the fact that the
#' weights were estimated (see the section "Estimated weights" of [df_counting()]).
#'
#' @param formula A formula for the treatment model, treatment (coded 0/1) on the left-hand side.
#'   Used with `data`; fitted by `stats::glm(formula, family = stats::binomial(), data = data)`
#'   with the default control, so that weights computed by hand from the same call agree.
#' @param data Data frame holding the variables of `formula`: exactly the rows that will be analysed.
#' @param type `"ipw"` (default): `1 / e` for treated subjects and `1 / (1 - e)` for controls, with
#'   `e` the fitted probability of treatment. `"stabilized"`: those weights multiplied by the
#'   proportion treated (treated subjects) and the proportion not treated (controls).
#' @param model An existing logistic `glm` fit, instead of `formula` and `data`. Give either
#'   `formula` and `data`, or `model`, not both.
#'
#' @details
#' With \eqn{x_i} the row of the model matrix of subject \eqn{i}, \eqn{A_i} the treatment and
#' \eqn{e_i} the fitted probability,
#' \deqn{IF_i = \Big\{\sum_j x_j x_j' e_j (1 - e_j)\Big\}^{-1} x_i (A_i - e_i)}{IF_i = { sum_j x_j x_j' e_j (1 - e_j) }^(-1) x_i (A_i - e_i)}
#' is the influence of subject \eqn{i} on the fitted coefficients (the estimate minus its limit is
#' approximately \eqn{\sum_i IF_i}{sum_i IF_i}), and the derivative of the weight of subject \eqn{i}
#' in the coefficients is
#' \eqn{dw_i = -\{(1 - e_i) / e_i\} x_i}{dw_i = -((1 - e_i) / e_i) x_i} for a treated subject and
#' \eqn{dw_i = \{e_i / (1 - e_i)\} x_i}{dw_i = (e_i / (1 - e_i)) x_i} for a control. With
#' `type = "stabilized"` the weights and `dw` are those of `"ipw"` multiplied by the arm proportion,
#' which is treated as a fixed number; `IF` is the same.
#'
#' The following are errors, each naming the problem: a `model` that is not a `glm`; a family other
#' than binomial or a link other than logit; prior weights other than 1 or an offset; a fit that did
#' not converge; a fitted probability within `1e-10` of 0 or 1; rows dropped for missing values; a
#' response that does not take both values 0 and 1.
#'
#' The weights are never truncated or rescaled. The rows of `weights`, `treat`, `ps`, `IF` and `dw`
#' are in the order of the rows of the data the model was fitted on.
#'
#' @return An object of class `"ps_weights"`, a list with
#' \describe{
#'   \item{weights}{The weights, one per row.}
#'   \item{type}{`"ipw"` or `"stabilized"`.}
#'   \item{treat}{The treatment (0/1), one per row.}
#'   \item{ps}{The fitted probabilities of treatment.}
#'   \item{IF}{Matrix, one row per subject and one column per estimated coefficient (aliased columns
#'     dropped): the influence of each subject on the coefficients.}
#'   \item{dw}{Matrix of the same shape: the derivative of each weight in the coefficients.}
#'   \item{model}{The `glm` fit.}
#'   \item{call}{The call.}
#' }
#' `print()` shows the type, the number of rows and of coefficients, and by arm the size, the range
#' of the weights and the effective sample size \eqn{(\sum w)^2 / \sum w^2}{(sum w)^2 / sum(w^2)}.
#'
#' @examples
#' library(survival)
#' gbsg$tte <- gbsg$rfstime / 30.4375
#' pw <- ps_weights(hormon ~ age + meno + size + grade + nodes + pgr + er, data = gbsg,
#'                  type = "stabilized")
#' pw
#' gbsg$sw <- pw$weights
#' range(gbsg$sw)
#'
#' @seealso [df_counting()] (section "Estimated weights"), [KM_diff()], [cox_rhogamma()]
#' @family survival_analysis
#' @importFrom stats glm binomial model.matrix fitted
#' @export
ps_weights <- function(formula = NULL, data = NULL, type = c("ipw", "stabilized"), model = NULL) {
  type <- match.arg(type)
  cl <- match.call()
  if (!is.null(model) && (!is.null(formula) || !is.null(data))) {
    stop("ps_weights(): give either 'formula' and 'data', or 'model', not both.", call. = FALSE)
  }
  if (is.null(model)) {
    if (is.null(formula) || is.null(data)) {
      stop("ps_weights(): 'formula' and 'data' are both needed (or give a fitted 'model').", call. = FALSE)
    }
    if (!inherits(formula, "formula")) stop("ps_weights(): 'formula' must be a formula.", call. = FALSE)
    if (!is.data.frame(data)) stop("ps_weights(): 'data' must be a data frame.", call. = FALSE)
    # glm.fit's own warnings (no convergence, fitted probabilities 0 or 1) are reported below as errors
    model <- withCallingHandlers(
      stats::glm(formula, family = stats::binomial(), data = data),
      warning = function(w) if (grepl("^glm\\.fit:", conditionMessage(w))) invokeRestart("muffleWarning"))
  }
  if (!inherits(model, "glm")) {
    stop("ps_weights(): 'model' must be a glm fit (a logistic regression for the treatment); it has class '",
         class(model)[1], "'.", call. = FALSE)
  }
  fam <- model$family
  if (!identical(fam$family, "binomial") || !identical(fam$link, "logit")) {
    stop(sprintf("ps_weights(): the model must be a logistic regression (family binomial, link logit); it has family '%s' and link '%s'.",
                 fam$family, fam$link), call. = FALSE)
  }
  if (!is.null(model$na.action) && length(model$na.action) > 0) {
    stop(sprintf("ps_weights(): %d row(s) were dropped by the model for missing values; fit the model on exactly the rows to be analysed (remove or impute the missing values first).",
                 length(model$na.action)), call. = FALSE)
  }
  pw <- model$prior.weights
  if (!is.null(pw) && any(pw != 1)) {
    stop("ps_weights(): the model has prior weights other than 1; fit it without 'weights' and with a 0/1 response.", call. = FALSE)
  }
  if (!is.null(model$offset) && any(model$offset != 0)) {
    stop("ps_weights(): the model has an offset; this is not supported.", call. = FALSE)
  }
  A <- model$y
  if (is.null(A)) stop("ps_weights(): the model does not hold its response; fit it with y = TRUE (the default).", call. = FALSE)
  A <- as.numeric(A)
  if (!all(A %in% c(0, 1)) || !any(A == 0) || !any(A == 1)) {
    stop("ps_weights(): the response of the model (the treatment) must take both values 0 and 1, and no other.", call. = FALSE)
  }
  if (!isTRUE(model$converged)) {
    stop("ps_weights(): the logistic regression did not converge.", call. = FALSE)
  }
  e <- as.numeric(stats::fitted(model))
  if (any(e < 1e-10 | e > 1 - 1e-10)) {
    stop(sprintf("ps_weights(): %d fitted probabilit%s numerically 0 or 1 (the covariates separate the arms); the weights are not defined.",
                 sum(e < 1e-10 | e > 1 - 1e-10), if (sum(e < 1e-10 | e > 1 - 1e-10) == 1) "y is" else "ies are"), call. = FALSE)
  }
  xm <- stats::model.matrix(model)
  est <- !is.na(stats::coef(model))
  xm <- xm[, est, drop = FALSE]
  attr(xm, "assign") <- NULL; attr(xm, "contrasts") <- NULL
  IF <- (xm * (A - e)) %*% solve(crossprod(xm * sqrt(e * (1 - e))))
  w <- ifelse(A == 1, 1 / e, 1 / (1 - e))
  dw <- ifelse(A == 1, -(1 - e) / e, e / (1 - e)) * xm
  if (type == "stabilized") {
    p1 <- mean(A)
    sc <- ifelse(A == 1, p1, 1 - p1)
    w <- sc * w
    dw <- sc * dw
  }
  dimnames(IF) <- dimnames(dw) <- list(NULL, colnames(xm))
  structure(list(weights = as.numeric(w), type = type, treat = A, ps = e, IF = IF, dw = dw, model = model, call = cl),
            class = "ps_weights")
}

#' @rdname ps_weights
#' @param x A `"ps_weights"` object.
#' @param digits Number of significant digits to print.
#' @param ... Not used.
#' @export
print.ps_weights <- function(x, digits = 4, ...) {
  n <- length(x$weights)
  cat(sprintf("Propensity-score weights (%s): n = %d, %d coefficient%s in the logistic model\n",
              x$type, n, ncol(x$IF), if (ncol(x$IF) == 1) "" else "s"))
  tab <- do.call(rbind, lapply(c(1, 0), function(a) {
    w <- x$weights[x$treat == a]
    data.frame(arm = if (a == 1) "treat = 1" else "treat = 0", n = length(w), min = min(w), max = max(w),
               ess = sum(w)^2 / sum(w^2), ratio = (sum(w)^2 / sum(w^2)) / length(w))
  }))
  names(tab) <- c("arm", "n", "weight min", "weight max", "ess", "ess / n")
  print(format(tab, digits = digits), row.names = FALSE)
  invisible(x)
}
