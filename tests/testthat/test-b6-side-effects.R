# B6: exported functions leave the caller's RNG state alone and do not print.

rng_state <- function() if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) get(".Random.seed", envir = globalenv()) else NULL

calls_with_draws <- function(d) list(
  KM_diff = function() KM_diff(d, "time", "status", "treat", draws.band = 30, show_resamples = FALSE),
  KM_diff_draws = function() KM_diff(d, "time", "status", "treat", weight.name = "w", draws = 30),
  calculate_risk_event_counts = function() calculate_risk_event_counts(d$time, d$status, d$w, sort(d$time), draws = 30),
  df_counting = function() df_counting(d, "time", "status", "treat", weight.name = "w", draws = 30, get.cox = FALSE),
  cox_rhogamma = function() cox_rhogamma(df_counting(d, "time", "status", "treat", get.cox = FALSE), draws = 30)
)

test_that("B6: caller's .Random.seed is unchanged after calls that draw", {
  d <- ws_sim(80, seed = 3, weighted = TRUE)
  for (nm in names(calls_with_draws(d))) {
    f <- calls_with_draws(d)[[nm]]
    set.seed(42)
    before <- rng_state()
    suppressWarnings(f())
    expect_identical(rng_state(), before, info = nm)
  }
})

test_that("B6: .Random.seed is not created when it did not exist", {
  d <- ws_sim(80, seed = 3, weighted = TRUE)
  old <- rng_state()
  on.exit(if (!is.null(old)) assign(".Random.seed", old, envir = globalenv()), add = TRUE)
  for (nm in names(calls_with_draws(d))) {
    f <- calls_with_draws(d)[[nm]]
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv())
    suppressWarnings(f())
    expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE), info = nm)
  }
})

test_that("B6: results for a given seedstart do not depend on the caller's seed", {
  d <- ws_sim(80, seed = 3)
  set.seed(1); a <- KM_diff(d, "time", "status", "treat", draws.band = 30, show_resamples = FALSE)
  set.seed(2); b <- KM_diff(d, "time", "status", "treat", draws.band = 30, show_resamples = FALSE)
  expect_identical(a$dhat_star, b$dhat_star)
})

test_that("B6: get_validated_weights() puts the missing mb_tstar in the error, prints nothing", {
  dfw <- data.frame(at_points = 1:5, S.pool = seq(1, 0.6, length.out = 5))
  out <- utils::capture.output(err <- tryCatch(get_validated_weights(dfw, scheme = "MB", scheme_params = list(rho = 0)),
                                               error = function(e) conditionMessage(e)))
  expect_length(out, 0)
  expect_match(err, "mb_tstar")
  expect_match(err, "rho")
})
