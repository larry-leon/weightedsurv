# B2: survival values in [0, 1]; exactly 0 (variance 0) once everything at risk has failed;
# no NaN in the log-transform limits.

last_event_sim <- function(seed) {
  old <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) get(".Random.seed", envir = globalenv()) else NULL
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else assign(".Random.seed", old, envir = globalenv()))
  set.seed(seed)
  n <- 40; tt <- sort(stats::rexp(n)); st <- stats::rbinom(n, 1, 0.7); st[n] <- 1
  list(time = tt, status = st, w = stats::runif(n, 0.2, 4))
}

test_that("KM_estimates(): weighted curve ending on an event ends at exactly 0, values in [0, 1]", {
  for (s in 1:60) {
    d <- last_event_sim(s)
    re <- calculate_risk_event_counts(d$time, d$status, d$w, d$time)
    km <- KM_estimates(re$ybar, re$nbar, re$sig2w_multiplier)
    expect_true(all(km$S_KM >= 0 & km$S_KM <= 1))
    expect_identical(utils::tail(km$S_KM, 1), 0)
    expect_identical(utils::tail(km$sig2_KM, 1), 0)
  }
})

test_that("df_counting(): survP, survG and the arm curves need no clipping and raise no false warnings", {
  d <- ws_sim(160, 21, weighted = TRUE)
  for (g in 0:1) { i <- which(d$treat == g); d$status[i[which.max(d$time[i])]] <- 1 }
  ws <- character(0)
  r <- withCallingHandlers(df_counting(d, "time", "status", "treat", weight.name = "w"),
                           warning = function(w) { ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning") })
  expect_false(any(grepl("outside \\[0,1\\]|NaNs produced", ws)))
  for (v in list(r$surv0, r$surv1, r$survP, r$survG)) expect_true(all(v >= 0 & v <= 1))
  expect_identical(min(r$surv0), 0)
  expect_no_error(check_km_curve(r$surv0))
})

test_that("km_quantile(): no NaN warning where S = 0", {
  d <- last_event_sim(3)
  re <- calculate_risk_event_counts(d$time, d$status, d$w, d$time)
  km <- KM_estimates(re$ybar, re$nbar, re$sig2w_multiplier)
  expect_no_warning(km_quantile(d$time, km$S_KM, sqrt(km$sig2_KM)))
})

test_that("plot_km_confint_polygon(): upper capped at 1, band closes onto the curve where S = 0, no NaN", {
  got <- NULL
  local_mocked_bindings(polygon = function(x, y, ...) got <<- y)
  x <- c(1, 2, 3, 4); s <- c(0.99, 0.6, 0.3, 0); se <- c(0.05, 0.1, 0.1, 0)
  plot_km_confint_polygon(x, s, se, 0.95, "grey")
  expect_false(anyNA(got))
  expect_true(all(got <= 1))
  expect_equal(got[c(4, 5)], c(0, 0))   # last lower and first reversed upper at S = 0
})
